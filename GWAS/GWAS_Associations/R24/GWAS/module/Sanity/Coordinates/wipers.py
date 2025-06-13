import re
import pandas as pd

from ...utils.columns import Columns
from ..base import WIPERNOTFOUND, ColumnWiper, CHROMOSOME_TYPE, INT_TYPE, _split_missing, ColumnWiperFactory

###############################################################
##################### Column Wipers ###########################
###############################################################
CHOROMOSOME_VALUES = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

class Chromosome(ColumnWiper):
    NAME: str =  Columns.CHR
    PARSER = CHROMOSOME_TYPE()

    def _main(self, data: pd.DataFrame) -> pd.DataFrame:
        nans, data = _split_missing(data, self.NAME)

        if not data.empty:
            ######## Remove non-digits at begining of string ###############
            regex = lambda x: re.sub(r"^[^\d+]+(?=([0-9]{1,2}|X|Y|MT))", "", x)
            data[self.NAME] = list(map(regex, data[self.NAME].astype(str)))

            ######## Remove '.0' at end of string caused by nan values in column ################
            data[self.NAME] = data[self.NAME].str.replace("\.0$","")

            ######## Convert 23 and 24 to X and Y respectively ##################
            data[self.NAME].replace({"23": "X", "24": "Y"}, inplace=True)
            
            ######## Remove everything that doesnt look like chromosome ###############
            data = self.PARSER.parse(data, self.NAME)

        data = pd.concat([data, nans])
        data[self.NAME] = pd.Categorical(data[self.NAME], categories=CHOROMOSOME_VALUES, ordered=True)

        return data

   
class Position(ColumnWiper):
    NAME: str = Columns.POS
    PARSER = INT_TYPE()

    def __init__(self) -> None:
        self.PARSER.TYPE = "UInt32"
        super().__init__()

    def _main(self, data: pd.DataFrame) -> pd.DataFrame:
        nans, data = _split_missing(data, self.NAME)
        ######## Remove '.0' at end of string caused by nan values in column ################
        if not data.empty:
            data[self.NAME] = data[self.NAME].astype(str).str.replace("\.0$","")

        data = pd.concat([data, nans])
        return self.PARSER.parse(data, self.NAME)
    


class WiperFactory(ColumnWiperFactory):
    WIPERS = (Columns.CHR, Columns.POS)
    TYPES = {Columns.CHR: Chromosome, Columns.POS: Position}

    @classmethod
    def new(cls, name: str) -> ColumnWiper:
        if not name in cls.WIPERS:
            raise WIPERNOTFOUND(f"Wiper for column {name} not available...")

        return cls.TYPES[name]()