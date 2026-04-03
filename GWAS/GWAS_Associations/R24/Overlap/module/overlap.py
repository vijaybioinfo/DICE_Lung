import pandas as pd
from .utils import Columns

def _overlap(x: pd.DataFrame, y: pd.DataFrame, on: str = Columns.SNPID) -> pd.DataFrame:
    return pd.merge(x, y, on=on)

def overlap(x: pd.DataFrame, y: pd.DataFrame, x_prefix: str = "", y_prefix: str = ""):
    collision = set(x.columns) & set(y.columns)
    collision.remove(Columns.SNPID)

    if x_prefix:
        x.rename(columns={c: x_prefix+"_"+c for c in collision}, inplace=True)
    
    if y_prefix:
        y.rename(columns={c: y_prefix+"_"+c for c in collision}, inplace=True)

    merge = _overlap(x, y)

    if merge.empty: return merge
    
    tmp = merge[Columns.SNPID].str.split(":", expand=True)                                                                                                           #type: ignore
    tmp.columns = [Columns.CHR, Columns.POS, Columns.REF, Columns.ALT]

    return pd.merge(merge, tmp, right_index=True, left_index=True)
    


class Overlap:
    def __call__(self, x: pd.DataFrame, y: pd.DataFrame, x_prefix: str = "", y_prefix: str = ""):
        merge = _overlap(x, y)

        if merge.empty: return merge
        
        tmp = merge[Columns.SNPID].str.split(":", expand=True)                                                                                                           #type: ignore
        tmp.columns = [Columns.CHR, Columns.POS, Columns.REF, Columns.ALT]

        return pd.merge(merge, tmp, right_index=True, left_index=True)


