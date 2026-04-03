import pandas as pd

from typing import Protocol
from .utils import Columns


#########################
###### QTL Columns ######
#########################
QTL_COLUMNS = {
               "eQTL": [Columns.SNPID, Columns.GENEID, Columns.GENENAME, Columns.PVAL, Columns.BETA],
               "chromQTL": [Columns.SNPID, Columns.PEAK, Columns.PVAL, Columns.BETA],
               "tfQTL": [Columns.SNPID, Columns.TF, Columns.PEAK, Columns.PVAL, Columns.BETA]
               }

QTL_TYPES = tuple(QTL_COLUMNS.keys())



class NOTAVAILABLEFORMATTER(Exception):
    pass


class Processing(Protocol):
    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        ...

class QTL(Processing):
    def __init__(self, keep: list[str]):
        self.keep = keep

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        if Columns.RSID in data.columns:
            cols.append(Columns.RSID)
        return data[self.keep].copy(deep=True)

class QTL_FORMATTER:
    AVAILABLE = QTL_TYPES

    COLUMNS = QTL_COLUMNS
    
    def get(self, which: str = "eQTL"):
        if which not in self.AVAILABLE:
            raise NOTAVAILABLEFORMATTER(f"Formatter {which} is not available!!!")
        
        return QTL(self.COLUMNS[which])


#class eQTLs(Processing):
#    def process(self, data: pd.DataFrame) -> pd.DataFrame:
#        cols = [Columns.SNPID, Columns.GENEID, Columns.GENENAME, Columns.PVAL, Columns.BETA]
#        if Columns.RSID in data.columns:
#            cols.append(Columns.RSID)
#        return data[cols].copy(deep=True)

#class chromQTLs(Processing):
#    def process(self, data: pd.DataFrame) -> pd.DataFrame:
#        cols = [Columns.SNPID, Columns.PEAK, Columns.PVAL, Columns.BETA]
#        if Columns.RSID in data.columns:
#            cols.append(Columns.RSID)
#        return data[cols].copy(deep=True)


class GWAS(Processing):
    def __init__(self, ld_engine = None, threshold: float = 5*10**-8) -> None:
        self.ld_engine = ld_engine
        self.threshold = threshold
    
    def _run_ld(self, data: pd.DataFrame):
        snps = data.loc[data[Columns.PVAL] <= self.threshold, Columns.SNPID].dropna().tolist()

        if self.ld_engine is None:
            return pd.DataFrame([snps, snps], columns = ["Lead", "LD"])

        if not len(snps):
            columns = ["Lead", "LD"] + list(self.ld_engine._params["populations"].keys())
            return pd.DataFrame(columns = columns)
        
        return self.ld_engine(snps)
    
    def _classify_ld_snps(self, data: pd.DataFrame) -> pd.DataFrame:
        lead = data[data["LD"] == data["Lead"]]
        lead["Category"] = "Lead"

        ld = data[~data["LD"].isin(lead["Lead"].unique().tolist())]
        ld.reset_index(drop=True, inplace=True)

        if not ld.empty:
            ld = ld.loc[ld.groupby("LD").apply(lambda x: (x[x.columns[2:]].max(axis=1)).idxmax())]                    #type: ignore
        
        ld["Category"] = "LD"

        return pd.concat([lead, ld])
    
    def _add_info(self, ld: pd.DataFrame, gwas: pd.DataFrame):
        ld.rename(columns={"LD": Columns.SNPID}, inplace=True)                                                                                                       
        ld = pd.merge(ld, gwas, on=Columns.SNPID, how="left")
        ld = ld.rename(columns={Columns.PVAL: f"GWAS_{Columns.PVAL}", Columns.BETA: f"GWAS_{Columns.BETA}"})                                                                      
        
        ld = pd.merge(ld, gwas.rename(columns={Columns.SNPID: "Lead"}), on="Lead", how="left")
        ld = ld.rename(columns={Columns.PVAL: f"GWAS_Lead_{Columns.PVAL}", Columns.BETA: f"GWAS_Lead_{Columns.BETA}"})                                                
        return ld.drop_duplicates()
    
    def process(self, data: pd.DataFrame):
        if not self.ld_engine is None:
            ld = self._run_ld(data)
            return self._add_info(self._classify_ld_snps(ld), data)
        
        ### No LD engine ###
        cols = [Columns.SNPID, Columns.PVAL, Columns.BETA]
        data = data.loc[data[Columns.PVAL] <= self.threshold, cols]
        data = data.rename(columns={Columns.PVAL: f"GWAS_Lead_{Columns.PVAL}", Columns.BETA: f"GWAS_Lead_{Columns.BETA}"}) 

        return data
