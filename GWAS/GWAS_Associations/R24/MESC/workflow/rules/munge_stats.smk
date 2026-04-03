rule munge_stats:
    input:
        GS.file
    output:
        "munge_stats/{disease}/{dataset}.sumstats.gz"
    params:
        software = CFG.software["MUNGE_STATS"],
        N = GS.N,
        prefix = lambda wildcards, output: output[0][:-12]
    conda: "ldsc"
    resources:
        mem_gb = 30,
        walltime = 60
    shell:
        "python {params.software} --sumstats {input} --signed-sumstats Z,0 --ignore SNPID --N {params.N} --out {params.prefix}"