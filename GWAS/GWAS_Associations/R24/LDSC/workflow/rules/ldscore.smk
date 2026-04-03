def get_snps(wildcards):
    snps = BASE.snps(wildcards)
    return "--print-snps " + snps if not snps is None else ""

rule LDScore:
    input:
        bfile = BASE.get_plink_files,
        annot = rules.single_annotation.output
    output:
        "ldscores/{project}/{cell}.{chr}.l2.ldscore.gz"
    params:
        bfile = lambda wildcards, input: input.bfile[0][:-4],
        prefix = "ldscores/{project}/{cell}.{chr}",
        snps = get_snps
    conda: "ldsc"
    resources:
        mem_gb = 50,
        walltime = 60
    shell:
        "python /home/jrocha/BioAdHoc/tools/ldsc/ldsc.py --l2 --bfile {params.bfile} --ld-wind-cm 1 --annot {input.annot} --out {params.prefix} {params.snps}"
