import pandas as pd

from .wipers import WiperFactory
from ...utils.columns import Columns
from ...utils.snpdb import ReferenceSNP
from ..base import RECOVER_METHOD_NOT_READY, RECOVERY_METHOD, SANITIZER_BUILDER
    
################################################################
####################### GWAS CATALOG ########################### 
################################################################
def _merge(x: pd.DataFrame, y, del_cols: list[str]) -> pd.DataFrame:
    return pd.merge(x.drop(columns=del_cols), y, how="left")

def _find_where_missing(data: pd.DataFrame) -> tuple[pd.DataFrame, ...]:
    where = data[[Columns.CHR, Columns.POS]].isna()
    
    both = where.all(axis=1)
    only_chr = where.any(axis=1) & ~where[Columns.POS]                                                             #type: ignore
    only_pos = where.any(axis=1) & ~where[Columns.CHR]                                                             #type: ignore
    no_missing = ~where.any(axis=1)

    return data[both], data[only_chr], data[only_pos], data[no_missing]                                        #type: ignore  

def retrieve_coordinates_using_rsid(miss: pd.DataFrame, db: str):
    database = ReferenceSNP(db)
    results = database.search_rsid(miss[Columns.RSID].tolist(), [Columns.CHR, Columns.POS, Columns.RSID])

    if results.empty: return miss

    both, chromosome, position, ok = _find_where_missing(miss)

    both = _merge(both, results, [Columns.CHR, Columns.POS])
    chromosome = _merge(chromosome, results, [Columns.CHR])
    position = _merge(position, results, [Columns.POS])

    return pd.concat([both, chromosome, position, ok]).drop_duplicates().reset_index(drop=True)

##################################################################################
########################## RECOVERY IMPLEMENTATION ###############################
##################################################################################
class CATALOG_RECOVERY_METHOD(RECOVERY_METHOD):
    FROM: str

    def check(self):
        if not hasattr(self, "FROM"):
            raise RECOVER_METHOD_NOT_READY("Method not ready, use builder class to properly iniitialize this method")

        if not isinstance(self.FROM, str):
            raise ValueError(f"FROM attribute must be of type {type(str)}, but {type(self.FROM)} found instead")

    def recover(self, data: pd.DataFrame) -> pd.DataFrame:
        where = data[Columns.RSID].isna()
        wrsid, norsid = data[~where], data[where]
        if not wrsid.empty:
            wrsid = retrieve_coordinates_using_rsid(wrsid, self.FROM)
        return pd.concat([wrsid, norsid]).reset_index(drop=True)
    
class NULL_RECOVERY_METHOD(RECOVERY_METHOD):
    def check(self): pass

    def recover(self, data: pd.DataFrame) -> pd.DataFrame: return data

####################################################################################
################################## BUILDER #########################################
####################################################################################

class NULL_COORDINATES_BUILDER(SANITIZER_BUILDER):
    def set_wiper_factory(self):
        self.product.FACTORY = WiperFactory()

    def set_recovery_method(self):                                                           
        self.product.RECOVERY_METHOD = NULL_RECOVERY_METHOD()

class CATALOG_COORDINATES_BUILDER(SANITIZER_BUILDER):
    def set_wiper_factory(self):
        self.product.FACTORY = WiperFactory()

    def set_recovery_method(self, _from: str):                                                           
        recovery = CATALOG_RECOVERY_METHOD()
        recovery.FROM = _from

        recovery.check()
        self.product.RECOVERY_METHOD = recovery


def make_coordinates_sanitizer(panel: str = None):                                         #type: ignore
    if panel:
        builder = CATALOG_COORDINATES_BUILDER()
        builder.set_recovery_method(panel)
    else:
        builder = NULL_COORDINATES_BUILDER()
        builder.set_recovery_method()

    builder.set_wiper_factory()
    return builder.get_sanitizer()