from typing import Union
from ..utils.columns import Columns

#############################################
################## PROGRAMS #################
#############################################
ALLELES_PROGRAM = [Columns.EA, Columns.NEA]
STATS_PROGRAM = [Columns.PVAL, Columns.BETA, Columns.OR, Columns.SE, Columns.Z, Columns.AF]
COORDINATES_PROGRAM = [Columns.CHR, Columns.POS]

#####################################################################################################
###################################### CLEANING PROGRAM #############################################
#####################################################################################################
class EMPTY_CLEANING_PROGRAM(Exception):
    pass

class PROGRAM:
    def __init__(self) -> None:
        self.values: list[str] = []

    def add(self, items: Union[list[str], str]):
        if isinstance(items, str):
            items = [items]
            
        for i in items:
            if not i in self.values:
                self.values.append(i)

    def reset(self):
        self.values = []

    def get(self):
        return self.values
    
    def empty(self): return not bool(self.values)