###########################################################
################## Summary stats results ##################
###########################################################
def get_ss_h2med_files(wildcards):
    samples = IP.get_datasets(wildcards.project)
    files = expand(rules.filter_summary_stats.output, project=wildcards.project, disease=wildcards.disease, dataset=wildcards.dataset, sample=samples)
    return {"h2med": files}

def merge_ss(wildcards):
    datasets = GS.get_by_category(wildcards)
    files = expand(rules.merge_ss_by_sample.output, project=wildcards.project, disease=wildcards.disease, dataset=datasets)
    return {"h2med": files}


rule merge_ss_by_sample:
    input:
        unpack(get_ss_h2med_files)
    output:
        h2med = OUT + "/filtered/{project}/summary_stats/{disease}/{dataset}.all.h2med"
    params:
        colname = "CELL"
    run:
        from R24.MESC.module.utils import merge_files

        h2med = merge_files(input.h2med, params.colname)
        h2med.to_csv(output.h2med, index=False, sep="\t")


use rule merge_ss_by_sample as merge_ss_by_dataset with:
    input:
        unpack(merge_ss)
    output:
        h2med = OUT + "/summary/{project}/{disease}/{disease}_{category}.all.h2med.ss"
    params:
        colname = "dataset"