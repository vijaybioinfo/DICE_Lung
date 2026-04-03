#################################################################
################# h2med non-summary statistics ##################
#################################################################
def gather_individual_scores(wildcards):
    scores = os.path.dirname(rules.expression_scores_non_summary.output[0])
    path = expand(scores, project=wildcards.project, sample=wildcards.sample)[0]
    return [f"{path}/{wildcards.sample}.{chrom}.lasso" for chrom in range(1,23)]



rule h2med_non_summary_stats:
    input:
        gwas = rules.munge_stats.output,
        scores = gather_individual_scores
    output:
        multiext(OUT + "/heritability/{project}/non_summary/{disease}/{dataset}/{sample}", ".all.h2med", ".categories.h2med")
    params:
        software = CFG.software["MESC"],
        exp = "individual_expression_scores/{project}/{sample}/{sample}",
        out = OUT + "/heritability/{project}/non_summary/{disease}/{dataset}/{sample}"
    conda: "mesc"
    resources:
        mem_gb = 50,
        walltime = 120
    script:
        "../scripts/h2med.py"


rule filter_non_summary_stats:
    input:
        OUT + "/heritability/{project}/non_summary/{disease}/{dataset}/{sample}.all.h2med"
    output:
        OUT + "/filtered/{project}/non_summary/{disease}/{dataset}/{sample}.all.h2med"
    run:
        import pandas as pd 

        data = pd.read_table(input[0])

        #if data.empty:
        #    data.to_csv(output[0], index=False, sep="\t")
        #else:
        #    h2 = data.loc[data.Quantity == "h2", "Estimate"].values[0]
        #    h2med = data.loc[data.Quantity == "h2med", "Estimate"].values[0]
        #
        #    if not((0 <= h2 <= 1) and (0 <= h2med <= 1)):
        #        data = pd.DataFrame(columns=data.columns)
        #        
        #    #data.to_csv(output[0], index=False, sep="\t")
        data.to_csv(output[0], index=False, sep="\t")



#####################################################################
################## h2med for meta analyzed scores ##################
#####################################################################
def gather_meta_scores(wildcards):
    files = expand(rules.meta_scores.output[0], project=wildcards.project, group=wildcards.group, chromosome=list(range(1, 23)))
    return files
    

use rule h2med_non_summary_stats as h2med_by_group with:
    input:
        gwas = rules.munge_stats.output,
        scores = gather_meta_scores
    output:
        multiext(OUT + "/heritability/{project}/groups/{disease}/{dataset}/{group}", ".all.h2med", ".categories.h2med")
    params:
        software = CFG.software["MESC"],
        exp ="meta_expression_scores/{project}/{group}/{group}",
        out = OUT + "/heritability/{project}/groups/{disease}/{dataset}/{group}"
    resources:
        mem_gb = 100,
        walltime = 120


use rule filter_non_summary_stats as filter_group with:
    input:
        rules.h2med_by_group.output 
    output:
         OUT + "/filtered/{project}/groups/{disease}/{dataset}/{group}.all.h2med"
