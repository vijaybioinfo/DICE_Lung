###### Overlap with QTLs ######
rule Overlap:
    input:
        ld = rules.LD.output,
        qtl_file = QTL.file
    output:
        raw = CFG.prefix + "/raw/{project}/{disease}/{dataset}/{cell}.tsv",
        filtered = CFG.prefix + "/overlap/{project}/{disease}/{dataset}/{cell}.tsv"
    params:
        qtl = CFG.type,
        raw = lambda wildcards, output: output.raw[:-4],
        filtered = lambda wildcards, output: output.filtered[:-4]
    threads: 1
    resources:
        mem_gb = 30,
        walltime = lambda wildcards, attempt: attempt * 60
    script:
        "../scripts/overlap.py"