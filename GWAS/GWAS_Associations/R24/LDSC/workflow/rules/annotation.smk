rule single_annotation:
    input:
        variants = BASE.annotation if BASE.has_base_annotation else BASE.get_bim_file,
        bed = IP.file
    output:
        "ldscores/{project}/{cell}.{chr}.annot.gz"
    params:
        type = BASE.has_base_annotation,
        prefix = "ldscores/{project}/{cell}.{chr}"
    resources:
        mem_gb = lambda wildcards, attempt: attempt * 20,
        walltime = 30
    script:
        "../scripts/annotation.py"