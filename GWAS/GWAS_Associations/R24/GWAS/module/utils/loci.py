import os
import pandas as pd

from R24.GWAS.module.utils.columns import Columns
from R24.GWAS.module.utils.misc import parent_directory

################################################################################
############################## Signal Selection ################################ 
################################################################################
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


################################################################################
############################## Loci Definition ################################# 
################################################################################
class ShrinkMethods:
    def significant(locus, threshold):
        sh = select_significat_snps(locus, threshold)
        mn, mx = sh[Columns.POS].min(), sh[Columns.POS].max()
        locus = locus[(locus[Columns.POS] >= mn) & (locus[Columns.POS] <= mx)]

        locus["Start"] = locus[Columns.POS].min()
        locus["End"] = locus[Columns.POS].max()
        return locus

    def edges(locus):
        locus["Start"] = locus[Columns.POS].min()
        locus["End"] = locus[Columns.POS].max()
        return locus
    
    def none(locus, left_edge, right_edge):
        locus["Start"] = left_edge
        locus["End"] = right_edge
        return locus


shrink_methods = {
    "significant": ShrinkMethods.significant,
    "edges": ShrinkMethods.edges,
    "none": ShrinkMethods.none
}

def define_loci(data: pd.DataFrame, threshold: float = 5*10**-8, window: int = 5*10**5, shrink: str = None):
    sig = select_significat_snps(data, threshold)
    top = select_top_signals(sig, window)
    
    loci = []
    for ix, position in enumerate(top[Columns.POS]):
        left_edge, right_edge =  position - round(window/2), position + round(window/2)
        locus = data[(data[Columns.POS] >= left_edge) & (data[Columns.POS] <= right_edge)]
        locus[Columns.LOCI] = ix + 1

        if shrink == "significant":
            locus = ShrinkMethods.significant(locus, threshold)
        elif shrink == "edges":
            locus = ShrinkMethods.edges(locus)
        else:
            locus = ShrinkMethods.none(locus, left_edge, right_edge)

        loci.append(locus)

    if not loci:
        return pd.DataFrame(columns=data.columns.tolist() + [Columns.LOCI, "Start", "End"])

    return pd.concat(loci)

###########################################################################
############################## Loci Merge ################################# 
###########################################################################
def merge_loci(data: pd.DataFrame):
    if data.empty:
        return pd.DataFrame(columns=data.columns.tolist() + [Columns.LOCI, "Start", "End"])

    loci = list(data.sort_values(["Start", "End"])[Columns.LOCI].unique())
    merge, new = [[loci.pop(0)]], []

    while loci:
        end = data.loc[data[Columns.LOCI] == merge[-1][-1], "End"].values[0]

        next = loci.pop(0)
        start = data.loc[data[Columns.LOCI] == next, "Start"].values[0]
        
        if end >= start:
            merge[-1].append(next)
        else:
            merge.append([next])

    for ix, m in enumerate(merge):
        tmp = data[data[Columns.LOCI].isin(m)]
        tmp[Columns.LOCI] = ix + 1
        tmp["Start"] = tmp["Start"].min()
        tmp["End"] = tmp["End"].max()
        new.append(tmp)
    
    return pd.concat(new).reset_index(drop=True)


################################################################################
######################## Loci Identifiction Methods ############################ 
################################################################################
class Identify_Loci:
    @staticmethod
    def standard(data: pd.DataFrame, threshold: float = 5*10**-8, window: int = 5*10**5, shrink: str = None):
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

        if not merge:
            return pd.DataFrame(columns=[Columns.LOCI, "Start", "End"])

        return pd.concat(merge).reset_index(drop=True)

    @staticmethod
    def pooled(data: list[pd.DataFrame], threshold: float = 5*10**-8, window: int = 10**5, shrink: str = None):
        data = pd.concat(data).drop_duplicates().reset_index(drop=True)
        
        return Identify_Loci.standard(data, threshold, window, shrink)


### For compatibility with older scripts ###
identify_loci = Identify_Loci.standard


def main(args):
    ############ Load Data ##############
    data = pd.read_table(args.file)

    ############ Define Loci ############
    loci = Identify_Loci.standard(data, args.threshold, args.window, args.shrink)

    os.makedirs(parent_directory(args.out), exist_ok=True)
    loci.to_csv(f"{args.out}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="Summary statistics file to estimate loci")
    parser.add_argument("--threshold", type=float, default=5*10**-8, help="Threshold to select significant snps, default = 5e-8 pvalue")
    parser.add_argument("--window", type=int, default=5*10**5, help="Window around the top snp to define locus")
    parser.add_argument("--shrink", default=None, choices=list(shrink_methods.keys()), help="Shrink boundaries to significant snps")
    parser.add_argument("--out", type=str, required=True, help="Output filename")

    args = parser.parse_args()

    main(args)




