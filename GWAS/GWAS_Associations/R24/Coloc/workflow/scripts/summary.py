import pandas as pd
from R24.Coloc.module.utils import make_coloc_matrix, merge_coloc_files

replication = lambda x: pd.DataFrame({"Feature2": x["Feature2"].unique(), "datasets": ", ".join(x["dataset1"].tolist()), "Number_of_datasets": x.shape[0]}, index=[0])

def make_empty_matrix(metadata: pd.DataFrame = None):
    if metadata is None or metadata.empty:
        cols = []
    else:
        cols = metadata.sample_name.tolist()

    return pd.DataFrame(columns=["Feature2", "datasets"] + cols)

def create_matrix(data: pd.DataFrame, metadata: pd.DataFrame = None) -> pd.DataFrame:
    if data.empty:
        return make_empty_matrix(metadata)

    matrix = make_coloc_matrix(data, metadata, "Feature2", "dataset2").fillna(0)
    
    core = matrix.index.to_frame().reset_index(drop=True)
    reps = data.drop_duplicates(["Feature2", "dataset1"])
    
    if reps.shape[0] == 1:
        reps = pd.DataFrame({"Feature2": reps["Feature2"], "datasets": reps["dataset1"], "Number_of_datasets": 1}, index=[0])
    else:
        reps = reps.groupby("Feature2", as_index=False).apply(replication).reset_index(drop=True)
    
    joined = pd.merge(core, reps)
    joined = pd.merge(joined, matrix.reset_index())

    return joined

def main(args):
    metadata = pd.read_csv(args.metadata)
    merged = merge_coloc_files(args.files, args.column)
    matrix = create_matrix(merged, metadata)

    merged.to_csv(args.out[0], index=False, sep="\t")
    matrix.to_csv(args.out[1], index=False, sep="\t")


if __name__ == "__main__":
    from collections import namedtuple

    args = namedtuple("args", ["files", "column","out"])
    args.files = snakemake.input                                                          #type: ignore
    args.column = snakemake.params.column                                                 #type: ignore
    args.metadata = snakemake.params.metadata                                             #type: ignore  
    args.out = snakemake.output                                                           #type: ignore

    main(args)


