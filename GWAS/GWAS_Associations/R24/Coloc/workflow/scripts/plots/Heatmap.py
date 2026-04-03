import numpy as np
import pandas as pd
import R24.Utils.Utilities as Utl
import R24.Utils.plots.Heatmap as HM

from pyensembl import EnsemblRelease
from typing import  Dict, Tuple, Any
from R24.Utils.Files import read_file
from R24.Utils.plots.Utilities import PlotPDF
from R24.Utils.Utilities import EnsembleGENEID, geneid_to_genename

####### Special Type Hints #######
Config = Dict[str, Any]
TupleConfig = Tuple[Config, Config, Config]


def create_matrix(coloc_files: str, threshold: float, which: int = 2) -> pd.DataFrame:
    coloc = []
    associations = []
    feature = f"Feature{which}"

    for file in coloc_files:
        name = Utl.trim_file_extension(Utl.get_filename(file))
        
        temp = read_file(file)
        #temporal fix to translate geneid to genename
        #if temp[feature].apply(lambda x: EnsembleGENEID.isID(x)).all():
        #    temp[feature] = temp[Col.GENENAME]
        if temp.empty:
            continue
        
        temp = temp.loc[temp.groupby([feature]).PPH4.idxmax()]
        temp["ID"] = name

        coloc.append(temp)
        associations.append(temp[temp["PPH4"] >= threshold])

    if not coloc:
        return pd.DataFrame()
        
    coloc = pd.concat(coloc).reset_index(drop=True)
    associations = pd.concat(associations)[[feature]].drop_duplicates()

    for name in coloc.ID.unique():
        data = coloc[coloc.ID == name].rename(columns={"PPH4": name})[[feature, name]]      #type: ignore
        data = data.groupby(feature).apply(lambda x: x.loc[x[name].idxmax()]).reset_index(drop=True)
        associations = pd.merge(associations, data, how="left").fillna(0)

    if associations[feature].apply(lambda x: EnsembleGENEID.isID(x)).all():
        associations = geneid_to_genename(associations, feature, feature)

    if not associations.empty:
        associations = associations.groupby(feature).apply(lambda x: x.max())

    return associations.set_index(feature)

def check_missing_columns(data, metadata):
    for col in metadata.sample_name:
        if not col in data.columns:
            data[col] = np.nan
    
    return data

def sort_matrix_by_pp4(data, threshold):
    fake = data.copy(deep=True)
    fake[fake < threshold] = 0
    _sorted = HM.sort.MatrixSorterMethods.standard(fake)

    return data.loc[_sorted.index]

def get_configuration(associations, metadata, args) -> TupleConfig:
    # configuration with PPH4 styler
    config: Config = {"plot": HM.style.HeatmapStyler.pp4()}
    config["plot"]["annotation_colors"] = ",".join(HM.get_colors(associations, metadata))
    config["plot"]["configuration"] = "labels"
    config["out"] = f"{args.out}_labels"

    config_nolabels: Config = Utl.deepcopy(config)
    config_nolabels["plot"]["configuration"] = "nolabels"
    config_nolabels["out"] = f"{args.out}_nolabels" 

    config_nolabels_noheaders: Config = Utl.deepcopy(config)
    config_nolabels_noheaders["plot"]["configuration"] = "nolabels-noheader"
    config_nolabels_noheaders["out"] = f"{args.out}_nolabels_noheader"

    return config, config_nolabels, config_nolabels_noheaders

######## Remove HLA locus ##########
def remove_HLA_locus(data: pd.DataFrame, genome: str = "hg19"):
    ensembl = EnsemblRelease(75 if genome == "hg19" else 108)

    if genome == "hg19":
        locus = (6, 28477897, 33448354)
    else:
        raise Exception("Unknown Locus HLA for genome version hg38")

    genes = [gene.gene_name for gene in ensembl.genes_at_locus(*locus)]
    data = data.loc[~data.index.isin(genes)]

    return data


def main(args):
    metadata = read_file(args.metadata)
    associations = create_matrix(args.coloc_files, args.threshold)
    associations = check_missing_columns(associations, metadata)
    associations = remove_HLA_locus(associations)

    if associations.empty:
        Utl.warning("No colocalization results found, writing empty plot")
        Utl.shell(f"touch {args.out}.pdf")
        return 

    # Reorder columns and sort
    associations = HM.reorder_columns(associations, metadata)
    associations = sort_matrix_by_pp4(associations, args.threshold)

    # configuration with PPH4 styler
    configs = get_configuration(associations, metadata, args)
    
    # Create Heatmap object
    heatmap = HM.Heatmap(data=associations, out=args.out)

    # Create and merge Heatmaps
    with PlotPDF(args.out) as pdf:
        heatmap(configs[0])
        heatmap(configs[1])
        heatmap(configs[2])

        files = [f"{cfg['out']}.pdf" for cfg in configs]
        pdf.extend(files)

        Utl.rm(files)


if __name__ == '__main__':
    if "snakemake" in globals():
        from collections import namedtuple

        args = namedtuple("args", ["coloc_files", "threshold", "feature", "metadata", "out"])           #type: ignore
        args.coloc_files = snakemake.input.coloc_files                                                  #type: ignore
        args.threshold = snakemake.params.threshold                                                     #type: ignore
        args.feature = snakemake.params.feature                                                         #type: ignore
        args.metadata = snakemake.params.metadata                                                       #type: ignore
        args.out = snakemake.params.prefix                                                              #type: ignore

    else:
        import argparse

        parser = argparse.ArgumentParser()

        parser.add_argument("--coloc-files", required=True, help="Comma-separated list of coloc raw files")
        parser.add_argument("--threshold", default=0.8, type=float, help="PPH4 threshold")
        parser.add_argument("--feature", default=2, type=int, choices=[1, 2], help="Feature to show. Choices {1, 2}")
        parser.add_argument("--metadata", required=True, help="Project metadata")
        parser.add_argument("--out", required=True, help="Output prefix" )

        args = parser.parse_args()
        args.coloc_files = args.coloc_files.split(",")

    
    main(args)