import os
import pandas as pd
import pyranges as pr

from typing import Union
from collections.abc import Iterable
from R24.GWAS.module.utils.columns import Columns



def read_haploblocks(file: str):
    return pd.read_table(file)

def read_summary_stats(file: str, chunksize: int = 10_000_000):
    for chunk in pd.read_table(file, chunksize=chunksize):
        yield chunk

class PartitionMethod:
    def _method(self, stats: pd.DataFrame, haploblocks: pd.DataFrame) -> pd.DataFrame: ...

    def __call__(self, stats: pd.DataFrame, haploblocks: pd.DataFrame) -> pd.DataFrame: ...

class StandardPartition:
    def _method(self, stats: pd.DataFrame, haploblocks: pd.DataFrame) -> pd.DataFrame:
        print("Start")
        cols = list(stats.columns)

        S = stats#.copy(deep=True)
        S["Start"] = S[Columns.POS]
        S["End"] = S[Columns.POS]
        
        S = pr.PyRanges(S)
        H = pr.PyRanges(haploblocks)
        P = H.join(S, slack=1).df

        if P.empty:
            return pd.DataFrame(columns=cols + ["Loci"])

        #P["Loci"] = P[Columns.CHR].astype(str).str[:] + ":" + P["Start"].astype(str).str[:] + "-" + P["End"].astype(str).str[:]
        P["Loci"] = [f"{c}:{s}-{e}" for c, s, e in zip(P[Columns.CHR], P["Start"], P["End"])]
        
        P = P.drop(columns=["Start", "End", "Start_b", "End_b"])
        print("End")
        return P
    
    def __call__(self, stats: Union[pd.DataFrame, list], haploblocks: pd.DataFrame) -> pd.DataFrame:
        if isinstance(stats, Iterable) and not isinstance(stats, pd.DataFrame):
            return pd.concat([self._method(chunk, haploblocks) for chunk in stats])                            #type: ignore
        return self._method(stats, haploblocks)

class eQTL_Partition(StandardPartition):
    def __call__(self, stats: pd.DataFrame, haploblocks: pd.DataFrame) -> pd.DataFrame:
        if isinstance(stats, Iterable) and not isinstance(stats, pd.DataFrame):
            tmp = pd.concat([self._method(chunk, haploblocks) for chunk in stats])                            #type: ignore
        else:
            tmp = self._method(stats, haploblocks)

        tmp["Loci"] = tmp[Columns.GENEID].astype(str).str[:] + "_" + tmp["Loci"].str[:]
        return tmp
    
METHODS = {"GWAS": StandardPartition(),
           "eQTL": eQTL_Partition()
           }


def main(args):
    print(args.S)
    haploblocks = read_haploblocks(args.H)
    stats = read_summary_stats(args.S)
    method = METHODS[args.method]
    
    PARTITION = method(stats, haploblocks)
    
    if args.filter:
        filtered = [loci for _, loci in PARTITION.groupby("Loci") if any(loci[args.metric] <= args.threshold)]
        if filtered:
            PARTITION = pd.concat(filtered)
        else:
            PARTITION = pd.DataFrame(columns=list(PARTITION.columns))
    
    folder = os.path.abspath(os.path.join(args.out, os.pardir))
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    #PARTITION.to_csv(args.out + ".tsv", index=False, sep="\t")
    PARTITION.to_csv(args.out, index=False, sep="\t")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--H", required=True, help="Path to Haploblock file. Format: Chromosome\tStart\tEnd")
    parser.add_argument("--S", required=True, help="Path to summary statistics file")
    parser.add_argument("--method", choices=METHODS.keys(), default="GWAS", help=f"One of {set(METHODS.keys())}")
    parser.add_argument("--filter", action="store_true", help="Use this flag if you which to filter your data")
    parser.add_argument("--metric", default=Columns.PVAL, help="Column to filter output loci data")
    parser.add_argument("--threshold", default=5*10**-8, type=float)
    parser.add_argument("--out", required=True, help="Output prefix")

    args = parser.parse_args()

    main(args)







