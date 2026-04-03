checkpoint preprocess:
    input:
        file = IP.file,
        annot = CFG.annotation,
        db = CFG.DB
    output:
        directory("preprocess/{project}/{sample}")
    params:
        out = "preprocess/{project}/{sample}/{sample}",
        N = IP.N
    threads: 1
    resources:
        mem_gb = 80,
        walltime = 1200
    script:
        "../scripts/prepare_mesc.py"


def get_qtl_file(wildcards):
    folder = checkpoints.preprocess.get(
                                    project=wildcards.project,
                                    sample=wildcards.sample
                                    ).output[0]

    return f"{folder}/{wildcards.sample}.{wildcards.chromosome}.tsv"


rule expression_scores_summary_stats:
    input:
        get_qtl_file
    output:
        multiext("expression_scores_ss/{project}/{sample}/{sample}.{chromosome}", ".gannot.gz", ".G", ".ave_h2cis")
    params:
        software = CFG.software["MESC"],
        prefix = "expression_scores_ss/{project}/{sample}/{sample}",
    conda: "mesc"
    threads: 1
    resources:
        mem_gb = 30,
        walltime = 60
    shell:
        "python {params.software} --eqtl-sumstat {input} --compute-expscore-sumstat --out {params.prefix}"
