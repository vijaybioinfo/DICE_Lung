import re
import pandas as pd

from abc import ABC, abstractmethod
from typing import Union, Type, Protocol
from .programs import PROGRAM, EMPTY_CLEANING_PROGRAM

###############################################################################################
#################################### DataType Parsers #########################################
###############################################################################################
def _split_missing(data: pd.DataFrame, NAME: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    where = data[NAME].isna()
    nans, data = data[where], data[~where]
    return nans, data

class DataTypeParserProtocol(Protocol):
    TYPE: Union[Type[str], Type[int], Type[float], str]
    REGEX: str

    @property
    def _engine(self) -> bool: ...

    def parse(self, data: pd.DataFrame, column: str) -> pd.DataFrame: ...


class DataTypeParser(DataTypeParserProtocol):
    TYPE: Union[Type[str], Type[int], Type[float], str]
    REGEX: str

    @property
    def _engine(self): return lambda x: bool(re.match(self.REGEX, x))

    def parse(self, data: pd.DataFrame, column: str) -> pd.DataFrame:
        nans, data = _split_missing(data, column)

        if not data.empty:
            data[column] = data[column].astype(str)
            data = data[data[column].apply(self._engine)]
            data[column] = data[column].astype(self.TYPE)                                        #type: ignore
        
        data = pd.concat([data, nans])
        
        return data


##################### Types Implementation ##########################
class NUMERIC(DataTypeParser):
    TYPE: Union[Type[int], Type[float], str]
    SIGNED: bool

    def __init__(self, signed: bool = False) -> None:
        self.SIGNED = signed
        super().__init__()


class INT_TYPE(NUMERIC):
    TYPE = int                                 #type: ignore

    @property
    def REGEX(self) -> str: return "^(-)?\d+$" if self.SIGNED else "^\d+$"


class FLOAT_TYPE(NUMERIC):
    TYPE = float                               #type: ignore

    @property
    def REGEX(self) -> str: return "^(-)?\d+(\.\d+)?(e-\d+)?$" if self.SIGNED else "^\d+(\.\d+)?(e-\d+)?$"


class NUCLEOTIDE_TYPE(DataTypeParser):
    TYPE = str                                  #type: ignore
    REGEX = "^[A|T|C|G]+$"                      #type: ignore


class CHROMOSOME_TYPE(DataTypeParser):
    TYPE = "category"                           #type: ignore
    REGEX = '^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$'    #type: ignore

#---------------------------------------------------------------------------#

####################################################################################################
########################################## WIPERS ##################################################
####################################################################################################
class WIPERNOTFOUND(Exception):
    pass

class ColumnWiper(Protocol):
    NAME: str
    PARSER: DataTypeParser
    
    def _main(self, data: pd.DataFrame) -> pd.DataFrame: ...

    def wipe(self, data: pd.DataFrame, dropna: bool=True) -> pd.DataFrame:
        if not self.NAME in data.columns:
            raise ValueError(f"Column '{self.NAME}' not in columns")
        
        cleaned = self._main(data)
        return cleaned if not dropna else cleaned[~cleaned[self.NAME].isna()]


class ColumnWiperFactory(Protocol):
    WIPERS: tuple[str, ...]
    TYPES: dict[str, Type[ColumnWiper]]

    @classmethod
    def new(cls, name: str) -> ColumnWiper: ...


#####################################################################################################
####################################### RECOVERY METHOD ############################################# 
#####################################################################################################
class RECOVER_METHOD_NOT_READY(Exception):
    pass

class RECOVERY_METHOD(Protocol):
    REQUIREMENTS = {}

    def check(self): ...
            
    def recover(self, data: pd.DataFrame) -> pd.DataFrame: ...

#####################################################################################################
######################################### Sanitizers ################################################
#####################################################################################################

#class SANITIZER(Protocol):
class SANITIZER:
    FACTORY: ColumnWiperFactory
    RECOVERY_METHOD: RECOVERY_METHOD

    def __init__(self) -> None:
        self.PROGRAM = PROGRAM()

    def check(self):
        if not hasattr(self, "FACTORY") or not hasattr(self, "RECOVERY_METHOD"):
            raise Exception("SANITIZER is not properly initialize, please use builder class to properly create this object...")
        
        #if not isinstance(self.FACTORY, ColumnWiperFactory):
        #    raise ValueError(f"FACTORY attribute must be of tyep ColumnWiperFactory...")
        
        #if not isinstance(self.RECOVERY_METHOD, type(RECOVERY_METHOD)):
        #    raise ValueError(f"RECOVERY_METHOD attribute must be of tyep RECOVERY_METHOD...")
        
        self.RECOVERY_METHOD.check()

    def _build_wipers(self):
        self.WIPERS = {W: self.FACTORY.new(W) for W in self.PROGRAM.get()}

    def _first_wipe(self, data: pd.DataFrame) -> pd.DataFrame:
        for C in self.PROGRAM.get():
            data = self.WIPERS[C].wipe(data, dropna=False).reset_index(drop=True)
        return data

    def _final_scrub(self, data: pd.DataFrame, drop_empty: list[str] = False) -> pd.DataFrame:                 #type: ignore
        return data[~data[drop_empty].isna().any(axis=1)] if drop_empty else data

    def clean(self, data: pd.DataFrame, dropna: Union[list[str], bool] = False) -> pd.DataFrame:
        print("[INFO] - Checking Cleaning Program...")
        if self.PROGRAM.empty():
            raise EMPTY_CLEANING_PROGRAM("Cleaning Program is empty...")
        else:
            print("[INFO] - \tAll good..")
        
        print("[INFO] - Building Column Wipers...")
        self._build_wipers()

        print("[INFO] - Perfoming First Wipe...")
        data = self._first_wipe(data).reset_index(drop=True)

        print("[INFO] - Perfoming Recovering when possible...")
        if not data.empty:
            data = self.RECOVERY_METHOD.recover(data).reset_index(drop=True)

        print("[INFO] - Final data scrubbing...")
        if dropna == True:
            dropna = self.PROGRAM.get()
        
        ############### RESET ###############
        self.WIPERS = {}
        self.PROGRAM.reset()

        return self._final_scrub(data, dropna)                                                                  #type: ignore

######## BUILDER #######
class SANITIZER_BUILDER(ABC):
    def __init__(self) -> None:
        self.reset()

    def reset(self): 
        self.product = SANITIZER() 
    
    @abstractmethod
    def set_wiper_factory(self): ...

    @abstractmethod
    def set_recovery_method(self, *args, **kwargs): ...

    def get_sanitizer(self) -> SANITIZER:
        self.product.check()
        return self.product
