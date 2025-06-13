import pandas as pd
from typing import Any

from datatable import fread
from .chain_map import ChainFiles
from .utils import TEMPDIR, Columns, get_temp_file, get_filename, read_file, exists, shell

################# Functions ####################
def check_liftover_exe():
    if not shell("which liftOver").returncode == 0:
        raise Exception("LiftOver executable doesn't exists...")

def _write_data(data: pd.DataFrame, path: str) -> str:
    data = data[[Columns.CHR, Columns.POS]]

    data[Columns.CHR] = "chr" + data[Columns.CHR].astype(str)
    data["End"] = data[Columns.POS] + 1
    data["ID"] = range(1, data.shape[0] + 1)
    
    temp = get_temp_file(path)
    data.to_csv(temp, index=False, sep="\t", header=False)

    return temp

def _split(x: pd.DataFrame, y: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:        #type: ignore
    x[Columns.CHR] = x[Columns.CHR].astype("category")
    y[Columns.CHR] = y[Columns.CHR].astype("category")
    
    results = []
    for chromosome in y[Columns.CHR].unique():
        results.append((x[x[Columns.CHR] == chromosome], y[y[Columns.CHR] == chromosome]))
    
    return results                                                                         #type: ignore

def read_liftover_file(file: str) -> pd.DataFrame:
    data = fread(file, header=False, columns={0: str, 1: int, 2: int, 3: int}).to_pandas()
    data.columns = [Columns.CHR, Columns.POS, "END", "ID"]
    data[Columns.CHR] = data[Columns.CHR].str.replace("chr", "")
    return data

def _map_results(original: pd.DataFrame, results: pd.DataFrame) -> pd.DataFrame:
    def _merge(tup):
        x, y = tup 
        merge = pd.merge(x, y[[Columns.POS, "ID"]], on="ID")
        merge[f"{Columns.POS}_x"] = merge[f"{Columns.POS}_y"]
        merge = merge.rename(columns={f"{Columns.POS}_x": Columns.POS})
        return merge.drop(columns=[f"{Columns.POS}_y", "ID"])


    original["ID"] = range(1, original.shape[0] + 1)

    if results.empty: return pd.DataFrame(columns=original.columns)

    return pd.concat(map(_merge, _split(original, results)))                                 #type: ignore


################ LiftOver ######################
class LiftOver:
    def __init__(self, chains: ChainFiles, temp_dir = "."):
        self.check()

        self.chains = chains
        self.workdir = TEMPDIR(f"{temp_dir}/LiftOver")

    def check(self):
        return check_liftover_exe()

    def _create_plan(self, version: str, target: str) -> list[tuple[Any, Any]]:
        path = self.chains.find_path(version, target)
        return [(path[i], path[i+1]) for i in range(len(path) - 1)]

    def _liftover(self, file: str, chain: str, out: str) -> str:
        shell(f"liftOver {file} {chain} {out}.mapped {out}.unmapped")
        
        if not exists(f"{out}.mapped", b=True):
            raise Exception("Output file not found...")

        return f"{out}.mapped"

    def _execute_plan(self, data: pd.DataFrame, plan: list[tuple[Any, Any]]) -> pd.DataFrame:
        if not plan: raise Exception("Execution plan is empty...")
        self.workdir.create()
        
        for s, t in plan:
            if not data.empty:
                chain_file = self.chains.get(s, t)

                tempfile = _write_data(data, self.workdir.path)
                file = self._liftover(tempfile, chain_file, out=f"{self.workdir.path}/{get_filename(tempfile)}.out.{s}.{t}")
                data = read_liftover_file(file)

        return data

    def __call__(self, data: pd.DataFrame, version: str, target: str) -> pd.DataFrame:
        if version == target:
            print("Original and Target versions are the same...")
            return data

        na = data[[Columns.CHR, Columns.POS]].isna().any(axis=1)
        ok, nans = data[~na], data[na]                  

        ok[Columns.CHR] = ok[Columns.CHR].astype(str)                                         #type: ignore                                       
        ok[Columns.POS] = ok[Columns.POS].astype(int)                                         #type: ignore
        
        output = self._execute_plan(ok, self._create_plan(version, target))                   #type: ignore
        results = _map_results(ok, output)                                                    #type: ignore

        return pd.concat([results, nans])

