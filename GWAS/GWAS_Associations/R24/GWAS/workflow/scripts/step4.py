import os
import pandas as pd

from R24.GWAS.module.utils.columns import Columns

##################################
########### READ FILES ###########
################################## 
def retrieve_data(files: list[str]):
    collector = []

    for file in files:
        for chunk in pd.read_csv(file, sep="\t", chunksize = 10**7, usecols=["SNPID", "RSID", "NEA", "EA", "P", "RA", "NRA"]):
            chunk = chunk[chunk.P < 5*10**-8]
            chunk["dataset"] = os.path.splitext(os.path.basename(file))[0]
            collector.append(chunk)

    table = pd.concat(collector)

    return table


    
#########################################
########### RISK ALLELE TABLE ###########
######################################### 
def make_risk_allele_table_from_files(table: pd.DataFrame):
    def _summarize(data: pd.DataFrame):
        return pd.DataFrame({
                            Columns.SNPID: data[Columns.SNPID].unique()[0],
                            Columns.RSID: data[Columns.RSID].unique()[0],
                            Columns.REF: data[Columns.NEA].unique()[0],
                            Columns.ALT: data[Columns.EA].unique()[0],
                            "REF_risk_counts": data["REF_risk"].sum(),
                            "ALT_risk_counts": data["ALT_risk"].sum(),
                            "Total_number_of_datasets": data.dataset.nunique(),
                            "REF_risk_datasets": ";".join(data.loc[data["REF_risk"], "dataset"].unique().tolist()),
                            "ALT_risk_datasets": ";".join(data.loc[data["ALT_risk"], "dataset"].unique().tolist())
                            }, index=[0])

    ### Who is risk? ###
    table["REF_risk"] = table[Columns.RA] == table[Columns.NEA]
    table["ALT_risk"] = table[Columns.RA] == table[Columns.EA]

    return table.groupby(Columns.SNPID).apply(_summarize).reset_index(drop=True)
                    

####################################
########### MAIN PROGRAM ###########
####################################
def main(args):
    snp_table = retrieve_data(args.files)
    risk_table = make_risk_allele_table_from_files(snp_table)
    risk_table.to_csv(args.out, index=False, sep="\t")


###########################################
########### SNAKEMAKE INTERFACE ###########
###########################################
def parse_snake_args(snake):
    class args:
        files = snake.input
        out = snake.output[0]

    return args


if __name__ == "__main__":
    main(parse_snake_args(snakemake))