import numpy as np
import pandas as pd

from R24.Utils.LD import SNPNOTFOUND
from R24.Utils.Files import read_file
from R24.Utils.plots.Locus import Locus
from R24.Utils.Headers import Columns as Col
from R24.Utils.plots.Utilities import pdf_merge
from R24.Utils.Utilities import rm, deepcopy, shell, warning, exists


class OutputIDGenerator:
    def __init__(self, function):
        self.calls = 0
        self.function = function
    
    def __call__(self, *args, **kwargs):
        self.calls += 1
        lab, nolab = self.function(*args, **kwargs)

        if "out" in lab.keys() and exists(lab["out"]+".pdf", True):
            lab["out"] = f"{lab['out']}.{self.calls}"
            nolab["out"] = f"{nolab['out']}.{self.calls}"

        return lab, nolab


def get_snp(data, column, how="min"):
    if how == "min":
        snp = data.loc[data[column].idxmin(), Col.SNPID]
    elif how == "max":
        snp = data.loc[data[column].idxmax(), Col.SNPID]
        print(data.loc[data[column].idxmax()])
    else:
        raise ValueError(f"{how} not supported, possible values {'min', 'max'}")
    return snp

def merge_locus_traits(locus1, locus2):
    merge = pd.merge(locus1, locus2, on=Col.SNPID, how="outer")
    merge.rename(columns={f"{Col.POS}_x": Col.POS, f"{Col.PVAL}_x":"Feature1", f"{Col.PVAL}_y":"Feature2"}, inplace=True)

    merge.loc[:, "Feature1"] = -np.log10(merge["Feature1"])
    merge.loc[:, "Feature2"] = -np.log10(merge["Feature2"])

    return merge[[Col.SNPID, Col.POS, "Feature1", "Feature2"]]

@OutputIDGenerator
def update_parameters(labels, nolabels, xcol=None, ycol=None, xlab=None, ylab=None, out=None):
    labels = deepcopy(labels)
    nolabels = deepcopy(nolabels)

    if not xcol is None:
        labels["plot"]["x_column"] = xcol
        nolabels["plot"]["x_column"] = xcol

    if not ycol is None:
        labels["plot"]["y_column"] = ycol
        nolabels["plot"]["y_column"] = ycol

    if not xlab is None:
        labels["plot"]["xlab"] = xlab

    if not ylab is None:
        labels["plot"]["ylab"] = ylab

    if not out is None:
        labels["out"] = out
        nolabels["out"] = f"{out}_nolabels"

    return labels, nolabels

def write_empty(out):
     shell(f"touch {out}.pdf")

def main(args):
    #get feature names
    ft1 = args.feature1
    ft2 = args.feature2

    hits = read_file(args.coloc)

    if hits.empty:
        warning("No coloc hits, writing empty plot")
        write_empty(args.out)
        return

    #excluding chromosome because there not LD information for them
    hits = hits[~hits[Col.CHR].isin(["X", "Y"])]
    hits = hits.drop_duplicates(["Feature1", "Feature2"])[["Feature1", "Feature2", Col.SNPID]].values

    trait1 = read_file(args.trait1)
    trait1 = trait1[trait1[ft1].isin([pair[0] for pair in hits])]

    trait2 = read_file(args.trait2)
    trait2 = trait2[trait2[ft2].isin([pair[1] for pair in hits])]

    labels = {"plot": {"height": 10, "width": 10, "snp_column": Col.SNPID, "labels": True, "no_labels": False}}
    nolabels = {"plot": {"height": 10, "width": 10, "no_labels": True, "labels": False}}
    
    temporary = []

    for f1, f2, snp in hits:
        #Extract and merge locus from each trait
        locus1 = trait1[trait1[ft1] == f1]
        locus2 = trait2[trait2[ft2] == f2]
        merge = merge_locus_traits(locus1, locus2)

        #get most significant snp
        #snp = get_snp(merge, "Feature1", how='max')
        labels["plot"]["snp"] = snp

        # Plot configuration parameters #
        try:
            plot = Locus(merge, snp, args.out)
        except SNPNOTFOUND:
            continue

        # Temporal fix to rename axis from GENEID to GENE NAMES
        lab1 = locus1[Col.GENENAME].unique()[0] if ft1 == Col.GENEID else f1
        lab2 = locus2[Col.GENENAME].unique()[0] if ft2 == Col.GENEID else f2

        p1, p1_nolabels = update_parameters(labels, nolabels, ycol="Feature1", ylab=f"-log10(P) {lab1}", out=f"{args.out}_{f1}")
        p2, p2_nolabels = update_parameters(labels, nolabels, ycol="Feature2", ylab=f"-log10(P) {lab2}", out=f"{args.out}_{f2}")
        p3, p3_nolabels = update_parameters(labels, nolabels, xcol="Feature2", xlab=f"-log10(P) {lab2}", ycol="Feature1", ylab=f"-log10(P) {lab1}", out=f"{args.out}_{f1}_{f2}")
        configurations = [p1, p1_nolabels, p2, p2_nolabels, p3, p3_nolabels]
        
        for conf in configurations:
            plot(conf)
            temporary.append(f"{conf['out']}.pdf")

    #merge pdf files
    if not temporary:
        write_empty(args.out)
        return 

    pdf_merge(temporary, args.out)
    rm(temporary)



if __name__ == "__main__":
    if "snakemake" in globals():
        from collections import namedtuple

        args = namedtuple("snakemake", ["trait1", "trait2", "feature1", "feature2", "out"])
        args.trait1 = snakemake.input.trait1                                                #type: ignore                                               
        args.trait2 = snakemake.input.trait2                                                #type: ignore
        args.coloc = snakemake.input.coloc                                                  #type: ignore
        args.feature1 = snakemake.params.feature1                                           #type: ignore
        args.feature2 = snakemake.params.feature2                                           #type: ignore
        args.out = snakemake.params.prefix                                                  #type: ignore
    else:
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--trait1", required=True, help="Summary statistics for trait1")
        parser.add_argument("--trait2", required=True, help="Summary statistics for trait2")
        parser.add_argument("--coloc", required=True, help="Colocalization results")
        parser.add_argument("--feature1", required=True, help="Feature for trait1")
        parser.add_argument("--feature2", required=True, help="Feature for trait2")
        parser.add_argument("--out", required=True, help="Output prefix")

        args = parser.parse_args()

    main(args)



