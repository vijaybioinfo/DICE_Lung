import os
import numpy as np
import pandas as pd

def combine_data(files: list[str]):
    data = []

    for file in files:
        tmp = pd.read_table(file)
        tmp["dataset"] = os.path.splitext(os.path.basename(file))[0]
        data.append(tmp)

    return pd.concat(data).reset_index(drop=True)

def make_matrix(data: pd.DataFrame, column: str):
    return pd.pivot_table(data, index="Category", columns="dataset", values=column).T

def main(args):
    data = combine_data(args.files)
    enrichment = make_matrix(data, "Enrichment").reset_index()
    pvalue = -np.log10(make_matrix(data, "Enrichment_p"))
    pvalue = pvalue.reset_index()

    data.to_csv(f"{args.out}.tsv", index=False, sep="\t")
    enrichment.to_csv(f"{args.out}.enrichment.mtx", index=False, sep="\t")
    pvalue.to_csv(f"{args.out}.pvalue.mtx", index=False, sep="\t")


if __name__ == "__main__":
    from collections import namedtuple

    args = namedtuple("args", "files out")

    args.files = snakemake.input
    args.out = snakemake.params.prefix

    main(args)

