import os
import re
import numpy as np
import pandas as pd

########################################
#########       READ DATA       ########
########################################
class Read_Data:
    @staticmethod
    def from_path(path: str):
        filename = path.split("/")[-1][:-4]
        disease, category = filename.rsplit("_", 1)

        data = pd.read_table(path)
        data["disease"] = disease
        data["disease_category"] = int(float(category))
        
        data.rename(columns={
                            "dataset1": "sample_name", 
                            "dataset": "sample_name",
                            "Feature2": "Feature",
                            "GENEID": "Feature"
                            }, 
                    inplace=True)

        return data

    @staticmethod
    def from_metadata(meta: pd.DataFrame, path: str):
        files  = []
        mapper = meta[["Disease", "disease"]].drop_duplicates().set_index("disease").to_dict()["Disease"]

        for D in mapper.keys():
            if os.path.exists(f"{path}/{D}"):
                files.extend([f"{path}/{D}/{i}" for i in os.listdir(f"{path}/{D}") if re.match(f"{D}_\d+.tsv", i) or re.match(f"{D}_\d+.0.tsv", i)])

        data = pd.concat(Read_Data.from_list(files))
        data["disease"] = data["disease"].replace(mapper)

        return data

    @staticmethod
    def from_list(l: list[str]):
        return [Read_Data.from_path(file) for file in l]

########################################
######### COUNT NUMBER OF GENES ########
########################################

def add_counts(meta: pd.DataFrame, paths: list[str], name: str):
    collect = []
    data = pd.concat([Read_Data.from_metadata(meta, path) for path in paths])

    for disease in meta.disease.unique():
        ### Read results ###
        subset = meta[meta.disease == disease]
        data = pd.concat([Read_Data.from_metadata(meta, path) for path in paths])

        ### Preprocess results ###
        data["sample_name"] = data["sample_name"].astype("category")
        data["sample_name"] = data["sample_name"].cat.set_categories(subset.sample_name.tolist())

        ### Get number of associations ###
        counts = data.groupby("sample_name").Feature.nunique().reset_index()
        counts.rename(columns={"Feature": name}, inplace=True)

        ### merge counts to metadata ###
        combined = pd.merge(subset[["disease", "sample_name"]], counts)

        assert subset.shape[0] == combined.shape[0], "ERROR: Some datasets are lost after estimaating gene counts per dataset!!!"

        ### collect ###
        collect.append(combined)

    
    count_table = pd.concat(collect)

    return count_table



######################################
######### SUPPLEMENTARY TABLE ########
######################################

### May require some improvements
def supplementary_table(new: pd.DataFrame, original: pd.DataFrame):
    original = original.rename(columns={
                                        "Sample name": "sample_name",
                                        "Genome version": "Original genome version"
                                        }
                                ).drop_duplicates(["Disease_Group", "disease", "sample_name"])

    joined = pd.merge(
                original,
                new[["Disease_Group", "disease", "sample_name"]].drop_duplicates(),
                how="left",
                indicator = "origin"
            )

    full = pd.concat(
                    [joined[joined.origin == "left_only"].drop(columns="origin"), new]
                    ).drop_duplicates()

    full["Link"] = [np.nan if s == 0 else l for s, l in zip(full["Sum stats"], full["Link"])]

    assert full.shape[0] == original.shape[0], "Some rows were lost when processing the data..."

    return full.sort_values(["Disease_Group", "Disease", "disease", "Order_GWAS"])


######################################
#########    MAIN FUNCTION    ########
######################################
def main(args):
    ### Read Data ###
    new = pd.read_csv(args.metadata)
    original = pd.read_csv(args.original)

    ### Get Gene Counts ###
    coloc_counts = add_counts(new[new["COLOC"]], args.coloc_path, "Number_of_associations_COLOC")
    overlap_counts = add_counts(new[new["OVERLAP"]], args.overlap_path, "Number_of_associations_OVERLAP")

    ### Merge Gene Counts ###
    combined = pd.merge(new, overlap_counts, how="left")
    combined = pd.merge(combined, coloc_counts, how="left")

    ### Get supplementary table ###
    supplementary = supplementary_table(combined, original)

    ### Create output file ###
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    supplementary.to_csv(f"{args.out}.csv", index=False)


############################
#########    CLI    ########
############################
def cli():
    import argparse 

    parser = argparse.ArgumentParser()

    parser.add_argument("--metadata", required=True, help="GWAS metadata file")
    parser.add_argument("--original", required=True, help="GWAS original metadata file")
    parser.add_argument("--coloc-path", nargs="+", help="Path to coloc summary results files")
    parser.add_argument("--overlap-path", nargs="+", help="Path to overlap summary results files")
    parser.add_argument("--out", required=True, help="Output prefix")

    args = parser.parse_args()

    return args



if __name__ == "__main__":
    args = cli()

    main(args)

