import pandas as pd

from numpy import nan
from ...utils.columns import Columns
from ..base import ColumnWiper, NUCLEOTIDE_TYPE, WIPERNOTFOUND, ColumnWiperFactory

###############################################################
##################### Column Wipers ###########################
###############################################################

class ALLELE_WIPER(ColumnWiper):                                                            #type: ignore
    PARSER = NUCLEOTIDE_TYPE()

    def _main(self, data: pd.DataFrame) -> pd.DataFrame:
        data[self.NAME].replace({"?": nan}, inplace=True)
        return self.PARSER.parse(data, self.NAME)

class AlleleWiperFactory(ColumnWiperFactory):
    WIPERS = (Columns.NEA, Columns.EA, Columns.REF, Columns.ALT)
    TYPES = {W: ALLELE_WIPER for W in WIPERS}

    @classmethod
    def new(cls, name: str) -> ColumnWiper:
        if not name in cls.WIPERS:
            raise WIPERNOTFOUND(f"Wiper for column {name} not available...")
        
        NEW = cls.TYPES[name]()
        NEW.NAME = name

        return NEW