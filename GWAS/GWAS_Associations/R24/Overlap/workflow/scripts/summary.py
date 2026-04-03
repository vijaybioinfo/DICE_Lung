import numpy as np
import pandas as pd

from R24.Overlap.module.utils import Columns, make_matrix, merge_files


###### Functions ########
replication = lambda x: pd.DataFrame({"GWAS_datasets": ", ".join(x.dataset.tolist()), "Number_of_GWAS_datasets": x.shape[0]}, index=[0])
    
def create_matrix(data: pd.DataFrame, metadata: pd.DataFrame):
    matrix = make_matrix(data, metadata)
    matrix = -np.log10(matrix).fillna(0)
    
    core = matrix.index.to_frame().reset_index(drop=True)
    reps = data.drop_duplicates([Columns.GENEID, "dataset"])

    if reps.shape[0] == 1:
        reps = pd.DataFrame({Columns.GENEID: reps[Columns.GENEID], "GWAS_datasets": reps["dataset"], "Number_of_GWAS_datasets": 1}, index=[0])
    else:
        reps = reps.groupby(Columns.GENEID).apply(replication).reset_index(level=0)

    joined = pd.merge(core, reps)
    joined = pd.merge(joined, matrix.reset_index())

    return joined

def create_empty_matrix(metadata: pd.DataFrame):
     core = ["GENEID", "GENENAME", "GWAS_datasets", "Number_of_GWAS_datasets"]
     cells = metadata.sort_values("order")["sample_name"].tolist()

     return pd.DataFrame(columns=core + cells)


def main(args):
    combined = merge_files(args.files)
    if not combined.empty:
        matrix = create_matrix(combined, pd.read_csv(args.metadata))
    else:
        matrix = create_empty_matrix(pd.read_csv(args.metadata))

    combined.to_csv(f"{args.out}.tsv", index=False, sep="\t")
    matrix.to_csv(f"{args.out}.mtx", index=False, sep="\t")
    

if __name__ == "__main__":
    from collections import namedtuple

    args = namedtuple("args", ["files", "ld","out"])
    args.files = snakemake.input                                                          #type: ignore
    args.metadata = snakemake.params.metadata                                             #type: ignore  
    args.out = snakemake.params.prefix                                                    #type: ignore

    main(args)
