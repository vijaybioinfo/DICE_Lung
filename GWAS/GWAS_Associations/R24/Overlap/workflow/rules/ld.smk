######## Estimate LD ##########
rule LD:
    input:
        gwas = GS.file
    output:
        "LD/{disease}/{dataset}.ld"
    params:
        LD = CFG.LD,
        prefix = lambda wildcards, output: output[0][:-3]
    threads: 1
    resources:
        mem_gb = lambda wildcards, attempt: attempt * 50,
        walltime = 360
    script:
        "../scripts/ld.py"