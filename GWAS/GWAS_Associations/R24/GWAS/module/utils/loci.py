import os
import pandas as pd

from R24.GWAS.module.utils.columns import Columns
from R24.GWAS.module.utils.misc import parent_directory

def select_significat_snps(data: pd.DataFrame, threshold: float = 5*10**-8):
    return data[data[Columns.PVAL] <= threshold]

def select_top_signals(data: pd.DataFrame, window=10**5):
    signals, distance = [], round(window/2)

    while (not data.empty):
        top = data.loc[data[Columns.PVAL].idxmin()]
        lower, upper = max(top[Columns.POS] - distance, 0), top[Columns.POS] + distance

        data = data[(data[Columns.POS] <= lower) | (data[Columns.POS] >= upper)]
        signals.append(top)

    if not signals:
        return pd.DataFrame(columns=data.columns.tolist())
    
    return pd.concat(signals, axis=1).T.reset_index(drop=True)

def define_loci(data: pd.DataFrame, threshold: float = 5*10**-8, window: int = 5*10**5, shrink: bool = False):
    sig = select_significat_snps(data, threshold)
    top = select_top_signals(sig, window)
    
    loci = []
    for ix, position in enumerate(top[Columns.POS]):
        locus = data[(data[Columns.POS] >= position - round(window/2)) & (data[Columns.POS] <= position + round(window/2))]

        if shrink:
            sh = select_significat_snps(locus)
            mn, mx = sh[Columns.POS].min(), sh[Columns.POS].max()
            locus = locus[(locus[Columns.POS] >= mn) & (locus[Columns.POS] <= mx)]

        locus[Columns.LOCI] = ix + 1
        locus["Start"] = locus[Columns.POS].min()
        locus["End"] = locus[Columns.POS].max()
        loci.append(locus)

    if not loci:
        return pd.DataFrame(columns=data.columns.tolist() + [Columns.LOCI, "Start", "End"])

    return pd.concat(loci)

def merge_loci(data: pd.DataFrame):
    if data.empty:
        return pd.DataFrame(columns=data.columns.tolist() + [Columns.LOCI, "Start", "End"])

    loci = list(data.sort_values(["Start", "End"])[Columns.LOCI].unique())
    merge, new = [[loci.pop(0)]], []

    while loci:
        end = data.loc[data[Columns.LOCI] == merge[-1][-1], Columns.POS].max()

        next = loci.pop(0)
        start = data.loc[data[Columns.LOCI] == next, Columns.POS].min()

        if end >= start:
            merge[-1].append(next)
        else:
            merge.append([next])

    for ix, m in enumerate(merge):
        tmp = data[data[Columns.LOCI].isin(m)]
        tmp[Columns.LOCI] = ix + 1
        tmp["Start"] = tmp[Columns.POS].min()
        tmp["End"] = tmp[Columns.POS].max()
        new.append(tmp)
    
    return pd.concat(new).reset_index(drop=True)

def identify_loci(data: pd.DataFrame, threshold: float = 5*10**-8, window: int = 10**5, shrink: bool = False):
    merge, add = [], 0

    for chromosome in data[Columns.CHR].unique():
        chr_data = data[data[Columns.CHR] == chromosome]

        ############ Define Loci ############
        loci = define_loci(chr_data, threshold, window, shrink)

        ############ Merge Overlapping Loci ###############
        loci = merge_loci(loci)

        if loci.empty: continue

        loci[Columns.LOCI] = loci[Columns.LOCI] + add
        add += loci[Columns.LOCI].nunique()
        
        merge.append(loci)

    return pd.concat(merge).reset_index(drop=True)


def main(args):
    ############ Load Data ##############
    data = pd.read_table(args.file)

    ############ Define Loci ############
    loci = identify_loci(data, args.threshold, args.window, args.shrink)

    os.makedirs(parent_directory(args.out))
    loci.to_csv(f"{args.out}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="Summary statistics file to estimate loci")
    parser.add_argument("--threshold", type=float, default=5*10**-8, help="Threshold to select significant snps, default = 5e-8 pvalue")
    parser.add_argument("--window", type=int, default=10**5, help="Window around the top snp to define locus")
    parser.add_argument("--shrink", action="store_true", help="Shrink boundaries to significant snps")
    parser.add_argument("--out", type=str, required=True, help="Output filename")

    args = parser.parse_args()

    main(args)




