def get_heritability_files(wildcards):
    return expand(rules.heritability.output, project=wildcards.project, disease=wildcards.disease, dataset=wildcards.dataset, cell=IP.get_datasets(wildcards.project))

def get_merged_files(wildcards):
    return [OUT + f"/filtered/{wildcards.project}/{wildcards.disease}/{dataset}.tsv" for dataset in GS.get_by_category(wildcards)]

rule merge_by_tissue:
    input:
        get_heritability_files
    output:
        OUT + "/filtered/{project}/{disease}/{dataset}.tsv"
    run:
        import os
        import pandas as pd

        associations = []

        for file in input:
            name = os.path.splitext(os.path.basename(file))[0]
            
            data = pd.read_table(file)
            data = data[data["Category"] == f"{name}L2_0"]
            associations.append(data)

        associations = pd.concat(associations).reset_index(drop=True)
        associations.Category = associations.Category.str[:-4]
        
        associations.to_csv(output[0], index=False, sep="\t")


rule summarize:
    input:
        get_merged_files
    output:
        summary = OUT + "/summary/{project}/{disease}/{disease}_{category}.tsv",
        enrichment = OUT + "/summary/{project}/{disease}/{disease}_{category}.enrichment.mtx",
        pvalue = OUT + "/summary/{project}/{disease}/{disease}_{category}.pvalue.mtx"
    params:
        prefix = OUT + "/summary/{project}/{disease}/{disease}_{category}"
    script:
        "../scripts/summary.py"