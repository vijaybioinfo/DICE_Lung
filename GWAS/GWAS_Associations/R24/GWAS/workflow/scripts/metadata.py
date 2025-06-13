import pandas as pd

from R24.GWAS.module.utils.misc import get_filename
from R24.GWAS.module.utils.columns import Columns as Col


######### Helper Functions ############
def _classify(columns: list[str], is_stats: bool, has_n: bool):
    if not Col.PVAL:
        raise Exception("Expected statistics are not present")
    
    OVERLAP = Col.PVAL in columns
    GARFIELD = OVERLAP and is_stats

    TWAS = all(c in columns for c in [Col.PVAL, Col.BETA, Col.SE, Col.Z]) and is_stats
    COLOC = (Col.BETA in columns and Col.SE in columns and is_stats) or (Col.PVAL in columns and Col.AF in columns and has_n and is_stats)

    return OVERLAP, GARFIELD, TWAS, COLOC 

def count_entries(file: str):
    with open(file, "r") as stream:
        stream.readline()
        return sum(1 for _ in stream)
    
def classify_and_QC(original: str, clean: str, is_stats: bool, has_n: bool):
    n_entries = count_entries(original)
    
    with open(clean, "r") as stream:
        header = stream.readline().strip("\n").split("\t")
        
        index = header.index("found")
        rows = tuple(line.strip().split("\t")[index] == "True" for line in stream)
        
    n_variants = len(rows)
    n_in_panel = sum(rows)

    qc = (n_entries, n_variants, n_in_panel, n_in_panel/n_variants if n_variants else 0)
    classification = tuple(i and bool(n_variants) for i in _classify(header, is_stats, has_n))

    return qc, classification



########### MAIN #####################
def main(args):
    analyses = ["OVERLAP", "GARFIELD", "TWAS", "COLOC"]
    qcs = ["Number_of_entries", "Number_of_final_variants", "Number_of_variants_in_panel", "Proportion_of_variants_in_panel"]

    #Read Sample Table
    ST = pd.read_csv(args.sample_table)
    ST.rename(columns={"Sample name": "sample_name"}, inplace=True)

    #Keep current disease and corresponding datasets
    ST = ST[ST["disease"] == args.disease]
    samples = [get_filename(file, ext=False) for file in args.files]
    ST = ST[ST["sample_name"].isin(samples)]
    
    #Annotate genome version and liftover
    genome_versions = pd.DataFrame({"sample_name": args.original_version.keys(), "Genome version": args.original_version.values()})
    ST = pd.merge(genome_versions, ST)
    
    ST = ST.rename(columns={"Genome version": "Original genome version"})
    ST["Genome version"] = args.version
    ST["LiftOver"] = ST["Genome version"] != ST["Original genome version"]

    ### Checkpoints ###
    assert len(args.sources) == len(args.files), "Source files and processed files are not the same"
    assert len(samples) == (ST.shape[0]), "Samples lost while trying to match original genome version"

    ### Get QC metrics ###
    for sample, source, file in zip(samples, args.sources, args.files):
        is_stats = bool(int(ST.loc[ST["sample_name"] == sample , "Sum stats"].values[0]))
        N = float(str(ST.loc[ST["sample_name"] == sample , "sample_size"].values[0]).replace(",",""))

        QCs, classification = classify_and_QC(source, file, is_stats, (N == N) and N > 0)
        ST.loc[ST["sample_name"] == sample , analyses] = classification             
        ST.loc[ST["sample_name"] == sample , qcs] = QCs                 

    #Write results
    ST.to_csv(args.out, index=False)


def parse_snake_args(snake):
    class args:
        sources = snake.input.sources                                       #type: ignore
        files = snake.input.files                                           #type: ignore
        sample_table = snake.input.sample_table                             #type: ignore
        disease = snake.params.disease                                      #type: ignore
        version = snake.params.version                                      #type: ignore
        original_version = snake.params.original_version                     #type: ignore
        out = snake.output[0]                                               #type: ignore  
    
    return args

if __name__ == "__main__":
    main(parse_snake_args(snakemake))
    
















