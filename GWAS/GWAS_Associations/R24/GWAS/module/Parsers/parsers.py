import re
import numpy as np
import pandas as pd
import datatable as dt

from abc import ABC, abstractmethod
from pandas.errors import IntCastingNaNError
from ..utils.columns import Columns, ColumnTypes

####################################################################
######################## COLUMN TYPES ##############################
####################################################################
TYPES = {n: ColumnTypes().__getattribute__(r) for r, n in vars(Columns).items() if not r.startswith("__")}

####################################################################
######################### Parser Protocol ##########################
####################################################################
class Parser(ABC):
    @abstractmethod
    def _parse(self, file: str) -> pd.DataFrame: ...

    def parse(self, file: str) -> pd.DataFrame:
        data = self._parse(file)

        for col in data.columns:
            if col in TYPES and not data[col].isna().any():
                try:
                    data[col] = data[col].astype(TYPES[col])
                except (ValueError, IntCastingNaNError):
                    pass

        return data

    
#####################################################################
#################### Summary Statistics Parser ######################
#####################################################################
class SummaryStatisticsParser(Parser):
    def _parse(self, file: str):
        data = dt.fread(file).to_pandas()

        if Columns.CHR in data.columns:
            where = data[Columns.CHR].isna()
            
            ok, nans = data[~where], data[where]
            ok[Columns.CHR] = ok[Columns.CHR].astype(str).str.replace("\.0", "")

            data = pd.concat([ok, nans])
            data[Columns.CHR].replace({"True": "1", True: "1"}, inplace=True)

        for col, dtype in data.dtypes.items():
            if dtype == object:
                data[col] = data[col].replace("b''", np.nan)
                data[col] = data[col].where(data[col].apply(type) != bytes, np.nan)                 #type: ignore

        return data
    
#####################################################################
########################## GWAS Catalog Parser ######################
#####################################################################

######################## Helper Functions ###########################
def _fill(l: list[str], length: int) -> list[str]:
    return l + [np.nan] * (length - len(l))          #type: ignore

def _match(pattern, items: list[str]):
    func = lambda x: re.search(pattern, x)
    return list(map(func, items))

def _match_risk_allele(snps: list[str]):
    matches = _match(re.compile("[ACTGN]+$|\?$"), snps)
    alleles = [m.group() if m else np.nan for m in matches]
    return [a if a != "?" else np.nan for a in alleles]

def _match_snp_id(ids: list[str]) -> list[str]:
    matches = _match(re.compile(".+(?=-[ACTGN\?]+$)|.+"), ids)
    return [match.group() if match else i for i, match in zip(ids, matches)]

def _match_chr_pos(ids: list[str]):
    pattern = re.compile('(?<![^A-z])(1[0-9]|2[0-4]|[1-9]|X|Y):\d+$')
    return [m.group() if bool(m) else None for m in _match(pattern, ids)]

def _is_chr_pos_id(ids: list[str]) -> list[bool]:
    return [bool(m) for m in _match_chr_pos(ids)]

def _handle_special_ids(data: pd.DataFrame) -> pd.DataFrame:
    where = _is_chr_pos_id(data[Columns.RSID].tolist())

    chr_pattern = re.compile('(1[0-9]|2[0-4]|[1-9]|X|Y)(?=:)')
    pos_pattern = re.compile('(?<=:)\d+')

    chr_pos = _match_chr_pos(data.loc[where, Columns.RSID].tolist())
    chr_match = [m.group() for m in _match(chr_pattern, chr_pos)]        #type: ignore
    pos_match = [m.group() for m in _match(pos_pattern, chr_pos)]        #type: ignore

    data.loc[where, Columns.CHR] = chr_match                                                            #type: ignore                                                       
    data.loc[where, Columns.POS] = pos_match                                                            #type: ignore
    data.loc[where, Columns.RSID] = np.nan
    
    return data

def _make_table(fields: list[list[str]]) -> pd.DataFrame:
    def _get_column_values(val: list[str], n: int):
        if len(val) < 2: return val[0]
        return _fill(fields[3], n)
    
    frame = pd.DataFrame(fields[:3]).T
    frame.columns = [Columns.CHR, Columns.POS, Columns.RSID]
    frame[Columns.EA] = _match_risk_allele(frame[Columns.RSID].tolist())
    frame[Columns.RSID] = _match_snp_id(frame[Columns.RSID].tolist())

    size = frame.shape[0]

    frame[Columns.NEA] = np.nan
    frame[Columns.PVAL] = _get_column_values(fields[3], size)
    frame[Columns.BETA] = _get_column_values(fields[4], size)
    frame[Columns.AF] = _get_column_values(fields[5], size)

    frame = frame.fillna(np.nan).replace("", np.nan)
    frame[Columns.AF].replace("NR", np.nan, inplace=True)
    frame[Columns.EA].replace("?", np.nan, inplace=True)

    frame = _handle_special_ids(frame)

    return frame

def _get_header_indexes(head: list[str], cols: list[str]) -> dict[str, int]:
    return {col: head.index(col) for col in cols if col in head}

def _parse_row(headers: dict[str, int], row: str):
    fields = row.strip("\n").split('\t')
    fields_kept = []

    for index in headers.values():
        value = re.sub("\s[x;]\s|;\s|\s", ";", fields[index])
        fields_kept.append(value.split(";"))
        
    return _make_table(fields_kept)

############################ Parser ################################

class CatalogGWASParser(Parser):   
    def _parse(self, file: str):
        with open(file, "r") as stream:
            headers = _get_header_indexes(
                                        head = stream.readline().strip("\n").split("\t"),
                                        cols = ["CHR_ID", "CHR_POS", "STRONGEST SNP-RISK ALLELE", "P-VALUE", "OR or BETA", "RISK ALLELE FREQUENCY"]
                                        )
            
            data = pd.concat([_parse_row(headers, row) for row in stream.readlines()])
        
        return data.reset_index(drop=True)
    
######################################################################
######################### Parser Constructor #########################
######################################################################
class PARSER_CONSTRUCTOR:
    @staticmethod
    def construct_sumstats_parser() -> Parser:
        return SummaryStatisticsParser()
    
    @staticmethod
    def construct_catalog_parser() -> Parser:
        return CatalogGWASParser()