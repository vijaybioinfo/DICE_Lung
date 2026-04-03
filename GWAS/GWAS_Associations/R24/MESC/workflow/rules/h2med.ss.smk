#############################################################
################# h2med summary statistics ##################
#############################################################
def gather_summary_stats(wildcards):
    folder = checkpoints.preprocess.get(
                                        project = wildcards.project,
                                        sample = wildcards.sample
                                        ).output[0]
    
    suffix = [file.replace(wildcards.sample, "")[1:] for file in os.listdir(folder) if file.endswith(".tsv")]
    chromosomes = set([suff.split(".", 1)[0] for suff in suffix])
    
    return expand(rules.expression_scores_summary_stats.output, project=wildcards.project, sample=wildcards.sample, chromosome=chromosomes)


rule h2med_summary_stats:
    input:
        gwas = rules.munge_stats.output,
        scores = gather_summary_stats
    output:
        multiext(OUT + "/heritability/{project}/summary_stats/{disease}/{dataset}/{sample}", ".all.h2med", ".categories.h2med")
    params:
        software = CFG.software["MESC"],
        exp = "expression_scores_ss/{project}/{sample}/{sample}",
        out = OUT + "/heritability/{project}/summary_stats/{disease}/{dataset}/{sample}"
    conda: "mesc"
    resources:
        mem_gb = 50,
        walltime = 120
    script:
        "../scripts/h2med.py"


rule filter_summary_stats:
    input:
        rules.h2med_summary_stats.output 
    output:
         OUT + "/filtered/{project}/summary_stats/{disease}/{dataset}/{sample}.all.h2med"
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