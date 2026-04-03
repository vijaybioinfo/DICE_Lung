
########################################################################
##################### Individual Expression Scores #####################
########################################################################
rule expression_scores_non_summary:
    input:
        expression = IP.file,
        genotype = multiext(CFG.genotype["COHORT"], ".bed", ".bim", ".fam"),
        panel = multiext(CFG.genotype["PANEL"], ".bed", ".bim", ".fam")
    output:
        "individual_expression_scores/{project}/{sample}/{sample}.{chromosome}.lasso"
    params:
        chromosome = "{chromosome}",
        genotype = CFG.genotype["COHORT"],
        panel = CFG.genotype["PANEL"],
        plink = CFG.software["PLINK"],
        software = CFG.software["MESC"],
        prefix = "individual_expression_scores/{project}/{sample}/{sample}"
    threads: 1
    resources: 
        mem_gb = lambda wildcards, attempt: 80 + 20 * attempt,
        walltime = 600
    conda: "mesc"
    shell:
        "python {params.software} --compute-expscore-indiv --plink-path {params.plink} --expression-matrix {input.expression} --exp-bfile {params.genotype} --geno-bfile {params.panel} --chr {params.chromosome} --out {params.prefix} --tmp {params.prefix}"


#########################################################################################
##################### Meta Analysis of Individual Expression Scores #####################
#########################################################################################
def get_score_files_noss(wildcards):
    samples = IP.get_group_samples(wildcards)
    files = expand(rules.expression_scores_non_summary.output, project=wildcards.project, sample=samples, chromosome=wildcards.chromosome)
    return files

def get_input_prefixes(wildcards):
    files = get_score_files_noss(wildcards)

    input_file = f"meta_expression_scores/input_files/{wildcards.project}/{wildcards.group}.{wildcards.chromosome}.txt"
    os.makedirs(os.path.dirname(input_file), exist_ok=True)
    
    with open(input_file, "w") as stream:
        for path in files:
            print(path.rsplit(f".{wildcards.chromosome}.lasso", 1)[0], file=stream)
    
    return input_file


rule meta_scores:
    input:
        scores = get_score_files_noss,
        prefixes = get_input_prefixes,
        panel = multiext(CFG.genotype["PANEL"], ".bed", ".bim", ".fam"),
    output:
        multiext("meta_expression_scores/{project}/{group}/{group}.{chromosome}", ".gannot.gz", ".G", ".ave_h2cis")
    params:
        chromosome = "{chromosome}",
        panel = CFG.genotype["PANEL"],
        software = CFG.software["META_ANALYSIS"],
        out = "meta_expression_scores/{project}/{group}/{group}"
    threads: 1
    resources: 
        mem_gb = 150,
        walltime = 600
    conda: "mesc"
    shell:
        "{params.software} --input-prefixes {input.prefixes} --bfile {params.panel} --chr {params.chromosome} --out {params.out}"
