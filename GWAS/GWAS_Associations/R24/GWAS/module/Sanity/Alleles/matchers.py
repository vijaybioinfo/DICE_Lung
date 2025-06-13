from typing import Protocol
import pandas as pd

from os.path import exists
from ...utils.columns import Columns
from ...utils.snpdb import ReferenceSNP

#################################################################################
########################### Panel Matching Protocol #############################
#################################################################################
class PANELNOTFOUND(Exception):
    pass

class NOTREADY(Exception):
    pass

class PANEL_MATCHER(Protocol):
    PANEL: str

    def is_ready(self) -> bool: return not hasattr(self, "REFERENCE_PANEL") and bool(self.PANEL)

    def add_panel(self, panel: str):
        if not exists(panel): raise PANELNOTFOUND(panel)  
        self.PANEL = panel

    def _match(self, data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]: ...

    def match(self, data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        if not self.is_ready():
            raise NOTREADY("Allele Matcher is not ready, Reference Panel is not set or is invalid. Provide by using add_panel methodl")
        
        return self._match(data)
    
#################################################################################
########################### PANEL DATABASE MATCH ################################
#################################################################################

######################### Helper Functions ######################################
def _match_query(data: pd.DataFrame, query: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    isna = data[Columns.CHR].isna()
    na, data = data[isna], data[~isna]                                                             
    
    data[Columns.CHR] = data[Columns.CHR].astype(str)
    data[Columns.POS] = data[Columns.POS].astype(int)
    na["Found"] = "left_only"

    merged = pd.merge(data, query, 
                            how="left", 
                            left_on=[Columns.CHR, Columns.POS], 
                            right_on=[Columns.CHR, Columns.POS], 
                            indicator="Found"
                            )
            
    merged = pd.concat([merged, na])
    found = merged.loc[merged.Found == "both", merged.columns != "Found"]                                               #type: ignore
    not_found = merged.loc[merged.Found == "left_only", data.columns]                                                   #type: ignore

    return found, not_found

######################### Implementation ##########################################
class PANEL_DATABASE_MATCH(PANEL_MATCHER):
    def _select_from_db(self, data: pd.DataFrame):
        database = ReferenceSNP(self.PANEL)
        chromosome = data[Columns.CHR].unique()[0]
        position = data[Columns.POS].tolist()

        return database.search_chromosome_position(chromosomes = chromosome, 
                                                   positions = position, 
                                                   columns = [Columns.CHR, Columns.POS, Columns.REF, Columns.ALT]
                                                  )

    def _match(self, data: pd.DataFrame):
        def split(data):                                                                     
            for chromosome in data[Columns.CHR].unique():
                chunk = data[data[Columns.CHR] == chromosome]
                yield chunk
    
        match = [_match_query(chunk, self._select_from_db(chunk)) for chunk in split(data)]

        return pd.concat([i[0] for i in match]), pd.concat([k[1] for k in match]) 