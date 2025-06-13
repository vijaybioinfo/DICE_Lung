import os
import random
import string
import numpy as np
import pandas as pd
import datatable as dt

from subprocess import run
from functools import partial, singledispatch

class Columns:
    CHR: str = "Chromosome"
    POS: str = "Position"

class TEMPDIR:
    def __init__(self, path: str):
        self.path = path

    @property
    def exists(self):
        return exists(self.path, b=True)

    def create(self):
        if not self.exists:
            mkdir(self.path)
    
    def rm(self):
        rm(self.path, r=True)

### Helper functions ###
shell = partial(run, shell=True)

def mkdir(path) -> None:
    os.makedirs(path, exist_ok=True)

@singledispatch
def rm(path, r=False):
    shell(f"rm -r {path}" if r else f"rm {path}")

@rm.register(list)
def _(path, r=False):
    files = " ".join(path)
    rm(files, r)

def trim_file_extension(path: str) -> str:
    return os.path.splitext(path)[0]

def exists(path, b=False):
    _exists = os.path.exists(path)
    if b: return _exists
    if not _exists: raise Exception(f"Path not found: {path}")

def get_filename(path: str) -> str:
    return os.path.basename(path)

def get_temp_file(path: str) -> str:
    filename = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    temp = trim_file_extension(get_filename(filename))
    return f"{path}/{temp}.temp"

def read_file(file: str):
    data = dt.fread(file).to_pandas()
    
    data[Columns.CHR] = data[Columns.CHR].astype(str).str.replace("chr", "")

    data[Columns.CHR].replace({"True": "1", 
                                True: "1", 
                                "23": "X", 
                                23: "X", 
                                "24": "Y", 
                                24: "Y"},
                                inplace=True
                            )

    data[Columns.CHR] = data[Columns.CHR].str.replace("\.0", "")
    data[Columns.CHR] = data[Columns.CHR].replace("", np.nan)
    data[Columns.POS] = pd.to_numeric(data[Columns.POS])

    for col, dtype in data.dtypes.items():
        if dtype == object:
            data[col] = data[col].replace("b''", np.nan)
            data[col] = data[col].where(data[col].apply(type) != bytes, np.nan)

    return data