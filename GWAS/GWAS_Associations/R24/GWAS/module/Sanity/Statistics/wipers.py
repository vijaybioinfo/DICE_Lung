import pandas as pd

from ...utils.columns import Columns
from ..base import FLOAT_TYPE, ColumnWiper, WIPERNOTFOUND, ColumnWiperFactory


###############################################################
###################### Column Wipers ##########################
###############################################################

class STATS_WIPER(ColumnWiper):                                                                           #type: ignore
    def _main(self, data: pd.DataFrame) -> pd.DataFrame:
        data = self.PARSER.parse(data, self.NAME)

        ### If empty recast to float ###
        if data[self.NAME].isna().all() and self.PARSER.TYPE == float:
            data[self.NAME] = data[self.NAME].astype(self.PARSER.TYPE)                                    #type: ignore
        return data
    
class STATS_WIPER_WITH_FILTER(ColumnWiper):                                                               #type: ignore
    def _main(self, data: pd.DataFrame) -> pd.DataFrame:
        data = self.PARSER.parse(data, self.NAME)

        ### If empty recast to float ###
        if data[self.NAME].isna().all() and self.PARSER.TYPE == float:
            data[self.NAME] = data[self.NAME].astype(self.PARSER.TYPE)                                    #type: ignore

        ######## Remove everything that is <0 or  >1 ###############
        data = data[(0 <= data[self.NAME]) & (data[self.NAME] <= 1) | (data[self.NAME].isna())]

        return data
    

class WiperFactory(ColumnWiperFactory):
    WIPERS = (Columns.PVAL, Columns.BETA, Columns.OR, Columns.SE, Columns.Z, Columns.AF)

    PARSERS = {
                Columns.PVAL: FLOAT_TYPE(),
                Columns.BETA: FLOAT_TYPE(signed=True),
                Columns.OR: FLOAT_TYPE(signed=True),
                Columns.SE: FLOAT_TYPE(),
                Columns.Z: FLOAT_TYPE(signed=True),
                Columns.AF: FLOAT_TYPE()
            }
    
    TYPES = {
            Columns.PVAL: STATS_WIPER_WITH_FILTER,
            Columns.BETA: STATS_WIPER,
            Columns.OR: STATS_WIPER,
            Columns.SE: STATS_WIPER,
            Columns.Z: STATS_WIPER,
            Columns.AF: STATS_WIPER_WITH_FILTER
            }     

    @classmethod
    def new(cls, name: str) -> ColumnWiper:
        if not name in cls.WIPERS:
            raise WIPERNOTFOUND(f"Wiper for column {name} not available...")
        
        NEW = cls.TYPES[name]()

        NEW.NAME = name
        NEW.PARSER = cls.PARSERS[name]

        return NEW
