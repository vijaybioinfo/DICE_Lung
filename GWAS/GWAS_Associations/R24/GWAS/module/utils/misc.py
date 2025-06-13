import os
import pandas as pd
import pyranges as pr

from subprocess import run
from typing import Callable
from functools import partial

from .columns import Columns


shell: Callable = partial(run, shell=True)

def mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def parent_directory(path: str) -> str:
    return os.path.abspath(os.path.join(path, os.pardir))

def get_filename(path: str, ext: bool = True) -> str:
    if ext:
        return os.path.basename(path)
    return os.path.splitext(os.path.basename(path))[0]
    

def trim_file_extension(path: str) -> str:
    return os.path.splitext(path)[0]


############# Locus summaries ##############
def _locus_summary(locus: pd.DataFrame, threshold: float = 5*10**-8):
    n_snps = locus[Columns.SNPID].nunique()
    n_sig = locus.loc[locus[Columns.PVAL] <= threshold, Columns.SNPID].nunique()

    locus = locus.loc[locus[Columns.PVAL].idxmin()]

    locus["n_sig_snps"] = n_sig
    locus["n_snps"] = n_snps

    return locus

def summarize_haploblocks(loci: pd.DataFrame, threshold: float = 5*10**-8):
    loci = loci[["Haploblock", Columns.SNPID, Columns.RSID, Columns.PVAL]]
    loci = loci.groupby(Columns.LOCI).apply(_locus_summary)
    return loci.reset_index(drop=True)

def summarize_loci(loci: pd.DataFrame, threshold: float = 5*10**-8):
    n_haploblocks = loci.groupby(Columns.LOCI).apply(lambda x: x["Haploblock"].nunique()).reset_index().rename(columns={0: "n_haploblocks"})

    loci = loci[[Columns.LOCI, Columns.SNPID, Columns.RSID, Columns.PVAL]]
    loci = loci.groupby(Columns.LOCI, sort=False).apply(_locus_summary).reset_index(drop=True)
    loci = pd.merge(loci, n_haploblocks)
    
    return loci.reset_index(drop=True)


def haploblocks_to_loci(haploblocks: pd.DataFrame, width: int = 5*10**5, summarize: bool = False):
    haploblocks = haploblocks.copy(deep=True)
    haploblocks[Columns.CHR] = haploblocks[Columns.HAPLOBLOCK].str.split(":").str[0]

    tmp = haploblocks[Columns.HAPLOBLOCK].str.split(":").str[1].str.split("-")
    haploblocks["Start"] = tmp.str[0]
    haploblocks["End"] = tmp.str[1]

    ranges = pr.PyRanges(haploblocks)
    merged = ranges.merge(slack=width)
    joined = merged.join(ranges, slack=1).df

    joined[Columns.LOCI] = joined[Columns.CHR].astype(str) + ":" + joined.Start.astype(str).str[:] + "-" + joined.End.astype(str).str[:]
    joined.drop(columns=["Start", "End", "Start_b", "End_b"], inplace=True)

    return joined


