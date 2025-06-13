import numpy as np
import pandas as pd

from typing import Protocol
from .matchers import PANEL_MATCHER
from ...utils.columns import Columns

####################################################
################### CONSTANTS ######################
####################################################
CATEGORIES = {
              "2-1": (1, 0, 0),
              "3-4": (0, 1, 0),
              "1-0": (0, 0, 1),
              "0-2": (0, 0, 1),
              "4-3": (1, 1, 0),
              "0-1": (1, 0, 1),
              "2-0": (1, 0, 1),
              "3-0": (0, 1, 1),
              "0-4": (0, 1, 1),
              "0-3": (1, 1, 1),
              "4-0": (1, 1, 1)
              }

NOTES = {
         "1-2": "Found",
         "2-1": "Flip",
         "3-4": "Complement",
         "1-0": "Fill",
         "0-2": "Fill",
         "4-3": "Flip-Complement",
         "0-1": "Flip-Fill",
         "2-0": "Flip-Fill",
         "3-0": "Complement-Fill",
         "0-4": "Complement-Fill",
         "0-3": "Flip-Complement-Fill",
         "4-0": "Flip-Complement-Fill"
         }

VALID = tuple(NOTES.keys())

ORDER = ("Found", "Flip", "Complement", "Fill", "Flip-Complement", "Flip-Fill", "Complement-Fill", "Flip-Complement-Fill", "Not Found")

####################################################
################ DNA SEQUENCE MAP ##################
####################################################
_dna_mapper = ('ACTGNactgnYRWSKMDVHBXyrwskmdvhbx', 'TGACNtgacnRYWSMKHBDVXrywsmkhbdvx')
_invalid_characters_string = ''.join(set( chr(x) for x in range(256) if chr(x) not in _dna_mapper[0]))
DNA_SEQ_MAP = (str.maketrans(_dna_mapper[0], _dna_mapper[1], _invalid_characters_string), )

####################################################
################ Allele Modifiers ##################
####################################################
def flip(d: pd.DataFrame) -> pd.DataFrame: 
    d["A1"], d["A2"] = d["A2"].to_list(), d["A1"].to_list()
    return d

def complement(d: pd.DataFrame) -> pd.DataFrame: 
    d["A2"] = d["A2"].astype(object).str.translate(*DNA_SEQ_MAP)                  #type: ignore
    d["A1"] = d["A1"].astype(object).str.translate(*DNA_SEQ_MAP)                  #type: ignore
    return d

def fill(d: pd.DataFrame) -> pd.DataFrame: 
    wherex, wherey = d["A1"].isna(), d["A2"].isna()

    d.loc[wherex, "A1"] = d.loc[wherex, Columns.ALT]
    d.loc[wherey, "A2"] = d.loc[wherey, Columns.REF]

    return d

#####################################################################
######################## CORE FUNCTIONS #############################
#####################################################################

def _match_alleles(target: "pd.Series[str]", opt1: "pd.Series[str]", opt2: "pd.Series[str]") -> "pd.Series[str]":
    ctarget = target.astype(object).str.translate(*DNA_SEQ_MAP)              

    a = target.isna().replace({True: 0, False: 5})
    b = (target == opt1).replace({True: 1, False: 5})
    c = (target == opt2).replace({True: 2, False: 5})
    d = (ctarget == opt1).replace({True: 3, False: 5})
    e = (ctarget == opt2).replace({True: 4, False: 5})

    match = pd.concat([a, b, c, d, e], axis=1).min(axis=1)
    match = match.astype(int).astype(str)
    return match

def _classify(data: pd.DataFrame, categories = CATEGORIES):
    x = _match_alleles(data[Columns.A2], data[Columns.REF], data[Columns.ALT])
    y = _match_alleles(data[Columns.A1], data[Columns.REF], data[Columns.ALT])
    
    data["category"] = x.str[:] + "-" + y.str[:]
    data["todo"] = data.category.map(categories)
    
    where = data["todo"].isna()
    nona, na = data[~where], data[where]
    na["todo"] = [(0, 0, 0)] * na.shape[0]

    return pd.concat([nona, na])

def _clean(data: pd.DataFrame) -> pd.DataFrame:
    cleaned, functions = [], [flip, complement, fill]
    
    for case in data["todo"].unique():
        tmp = data[data["todo"] == case]
        
        methods = [f for i, f in zip(case, functions) if i]
        
        for method in methods:
            tmp = method(tmp)                                                           

        cleaned.append(tmp)

    return pd.concat(cleaned) if cleaned else pd.DataFrame(columns=data.columns) 

def _choose_case(data: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    def _fun(tmp):
        tmp["GroupID"] = tmp.groupby(cols, dropna=False).ngroup()

        size = tmp.groupby("GroupID").size()
        mins = tmp.groupby("GroupID")["Notes"].min()
        tmp.sort_values(["GroupID", "Notes"], inplace=True)

        tmp["Mins"] = list(np.repeat(mins, size))
        tmp = tmp[tmp["Notes"] == tmp["Mins"]]

        w = tmp["Notes"] == "Not Found"
        tmp = pd.concat([tmp[~w], tmp[w].drop_duplicates("GroupID")])
        tmp.drop(columns=["GroupID", "Mins"], inplace=True)

        return tmp
    
    if data.empty: return data
    
    return pd.concat([_fun(data[data[Columns.CHR] == chrom]) for chrom in data[Columns.CHR].unique()])

################################################################
##################### ALLELE HARMONIZER ########################
################################################################
class NOTREADY(Exception):
    pass

class ALLELE_HARMONIZER(Protocol):

    def is_ready(self) -> bool: ...

    def _harmonize(self, data: pd.DataFrame) -> pd.DataFrame: ...

    def harmonize(self, data: pd.DataFrame) -> pd.DataFrame:
        if not self.is_ready():
            raise NOTREADY("Allele Harmonizer is not ready. Use proper builder class to initialize this method")
        
        data = self._harmonize(data)
        return data

################### BiAllelic Implementation #####################
class BIALLELIC_HARMONIZER(ALLELE_HARMONIZER):
    def is_ready(self) -> bool: return True

    def _harmonize(self, data: pd.DataFrame) -> pd.DataFrame:
        data = data.drop_duplicates().reset_index(drop=True)
        tmp_cols = [c for c in data.columns if not c == Columns.REF and not c == Columns.ALT]
        
        ################## Classify each allele case match ###################
        new = _classify(data)
        new["Notes"] = new["category"].map(NOTES).fillna("Not Found")
        new["Notes"] = new["Notes"].astype("category").cat.set_categories(ORDER)
        new["Notes"] = new["Notes"].cat.as_ordered()
        
        ################## Check Duplicates (SNPs in different cases) #######################
        where = new.duplicated(tmp_cols, keep=False)                                                                       
        ok, dups = new[~where], new[where]
        
        ################## Choose which case to keep #####################
        dups = _choose_case(dups, tmp_cols)                                                                            
        new = pd.concat([ok, dups])

        ################## Solve each case ######################
        new = _clean(new)
        new["found"] = new["category"].isin(VALID)
        new.drop(columns=["category", "todo"], inplace=True)

        ################## Flag not found cases ######################
        new.loc[new["Notes"] == 'Not Found', [Columns.REF, Columns.ALT]] = np.nan
        
        return new.drop_duplicates()


################### MultiAllelic Implemantation ###############################
class MULTI_ALLELE_HARMONIZER(ALLELE_HARMONIZER):
    def is_ready(self) -> bool: return True 

    def _harmonize(self, data: pd.DataFrame):
        data = data.drop_duplicates().reset_index(drop=True)
        tmp_cols = [c for c in data.columns if not c == Columns.REF and not c == Columns.ALT]

        ################## Classify each allele case match ###################
        tmp = _classify(data)

        ################## Special Handling for both Alleles Missing ###################
        where = tmp["category"] == "0-0"
        multi, new = tmp[where], tmp[~where]

        multi["A2"] = multi[Columns.REF]
        multi["A1"] = multi[Columns.ALT]

        multi["Notes"] = "Fill"
        multi["found"] = True

        ################# Standard Handling ######################
        new["Notes"] = new["category"].map(NOTES).fillna("Not Found")                                                  
        new["Notes"] = new["Notes"].astype("category").cat.set_categories(ORDER)                                       
        new["Notes"] = new["Notes"].cat.as_ordered()                                                                        
        
        ################## Check Duplicates (SNPs in different cases) #######################
        where = new.duplicated(tmp_cols, keep=False)                                                                        
        ok, dups = new[~where], new[where]
        
        ################## Choose which case to keep #####################
        dups = _choose_case(dups, tmp_cols)                                                                            
        new = pd.concat([ok, dups])

        ################## Solve each case ######################
        new = _clean(new)                                                                                            
        new["found"] = new["category"].isin(VALID)

        ################## Flag not found cases ######################
        combined = pd.concat([new, multi])
        combined.drop(columns=["category", "todo"], inplace=True)   
        combined.loc[combined["Notes"] == 'Not Found', [Columns.REF, Columns.ALT]] = np.nan 

        return combined.drop_duplicates()