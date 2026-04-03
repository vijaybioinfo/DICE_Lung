###########################################################
################# non-summary statistics ##################
###########################################################
def get_h2med_non_sumary(wildcards):
    samples = IP.get_datasets(wildcards.project)
    files = expand(rules.filter_non_summary_stats.output, project=wildcards.project, disease=wildcards.disease, dataset=wildcards.dataset, sample=samples)
    return {"h2med": files}

def merge_h2med_non_summary(wildcards):
    datasets = GS.get_by_category(wildcards)
    files = expand(rules.merge_non_summary.output, project=wildcards.project, disease=wildcards.disease, dataset=datasets)
    return {"h2med": files}


rule merge_non_summary:
    input:
        unpack(get_h2med_non_sumary)
    output:
        h2med = OUT + "/filtered/{project}/non_summary/{disease}/{dataset}.all.h2med"
    params:
        colname = "CELL"
    run:
        from R24.MESC.module.utils import merge_files

        h2med = merge_files(input.h2med, params.colname)
        h2med.to_csv(output.h2med, index=False, sep="\t")


use rule merge_non_summary as merge_non_summary_by_dataset with:
    input:
        unpack(merge_h2med_non_summary)
    output:
        h2med = OUT + "/summary/{project}/{disease}/{disease}_{category}.all.h2med.noss"
    params:
        colname = "dataset"



####################################################################
################# non-summary statistics by group ##################
####################################################################
def get_group_h2med_files(wildcards):
    groups = IP.Groups(True)
    files = expand(rules.filter_group.output, project=wildcards.project, disease=wildcards.disease, dataset=wildcards.dataset, group=groups)
    return {"h2med": files}

def merge_h2med_non_summary_groups(wildcards):
    datasets = GS.get_by_category(wildcards)
    files = expand(rules.merge_non_summary_by_group.output, project=wildcards.project, disease=wildcards.disease, dataset=datasets)
    return {"h2med": files}


use rule merge_non_summary as merge_non_summary_by_group with:
    input:
        unpack(get_group_h2med_files)
    output:
        h2med = OUT + "/filtered/{project}/groups/{disease}/{dataset}.all.h2med"


use rule merge_non_summary as merge_non_summary_groups_by_dataset with:
    input:
        unpack(merge_h2med_non_summary_groups)
    output:
        h2med = OUT + "/summary/{project}/{disease}/{disease}_{category}.all.h2med.groups"
    params:
        colname = "dataset"