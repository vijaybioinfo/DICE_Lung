import pandas as pd
import sqlite3 as sql

from numpy import int64
from typing import Any, List, Union

class DataBase:
    def __init__(self, database: str):
        self.file = database

    def __enter__(self):
        self.conn = sql.connect(self.file)
        return self.conn

    def __exit__(self, type, value, traceback):
        self.conn.close()


class ReferenceSNP(DataBase):
    def _search(self, where: str, 
                    value: Union[Any, List[Any]],
                    columns: Union[str, List[str]] = "*") -> pd.DataFrame:

        columns = ", ".join(columns) if isinstance(columns, list) else columns
        base_command = f"SELECT {columns} FROM Variants WHERE {where} " 
        
        ### Check that input is not empty
        if len(value) == 0:
            raise Exception("Empty query !!!")
        
        if isinstance(value, list) and len(value) == 1:
            value = value[0]

        if isinstance(value, str):
            base_command = base_command + f"= '{value}' ;"

        elif len(value) == 1:
            base_command = base_command + f"= {value} ;"

        else:
            base_command = base_command + f"in {str(tuple(value))}"

        with self as connection:
            results = pd.read_sql(base_command, connection)
        
        return results

    def search_rsid(self, value: Union[str, List[str]], columns: Union[str, List[str]] = "*") -> pd.DataFrame:
        return self._search(where="RSID", value=value, columns=columns)
        
    def search_snpid(self, value: Union[str, List[str]], columns: Union[str, List[str]] = "*") -> pd.DataFrame:
        return self._search(where="SNPID", value=value, columns=columns)

    def search_chromosome(self, value: Union[int, List[int]], columns: Union[str, List[str]] = "*") -> pd.DataFrame:
        return self._search(where="Chromosome", value=value, columns=columns)

    def search_chromosome_position(self, chromosomes: Union[str, List[str]],
                                         positions: Union[int, List[int]],
                                         columns: Union[str, List[str]] = "*"
                                  ) -> pd.DataFrame:

        if isinstance(chromosomes, list) and len(chromosomes) == 1:
            chromosomes = chromosomes[0]

        if isinstance(positions, list) and len(positions) == 1:
            positions = int(positions[0])
        
        chr_cmd = f"Chromosome = '{chromosomes}'" if isinstance(chromosomes, str) else f"Chromosome in {str(tuple(chromosomes))}"
        pos_cmd = f"Position = {int(positions)}" if isinstance(positions, (int, int64, float)) else f"Position in {str(tuple(positions))}"                       #type: ignore

        full_cmd = '''SELECT {columns} FROM Variants WHERE {chr_cmd} and {pos_cmd}'''.format(columns=", ".join(columns), chr_cmd=chr_cmd, pos_cmd=pos_cmd)                                                       
        
        with self as connection:
            query = pd.read_sql(full_cmd, connection)

        return query