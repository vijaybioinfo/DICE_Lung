def get_ldscores(wildcards):
    return expand(rules.LDScore.output, project=wildcards.project, cell=wildcards.cell, chr=BASE.Chromosomes)

rule heritability:
    input:
        sumstats = rules.munge_stats.output,
        model = get_ldscores,
        weights = BASE.weights,
        frq = BASE.frequencies
    output:
        OUT + "/heritability/{project}/{disease}/{dataset}/{cell}.results"
    params:
        ref_ld = lambda wildcards, input: re.sub(r"\.(.+)?(\d+\.l2\.ldscore\.gz)", ".", input.model[0]),
        w_ld = lambda wildcards, input: re.sub(r"\.(([A-z]+)?\d+)?(\d+\.l2\.ldscore\.gz)", ".", input.weights[0]),
        frqfile = lambda wildcards, input: re.sub(r"(([A-z]+)?\d+)?(\.frq)$", "", input.frq[0]),
        prefix = lambda wildcards, output: output[0][:-8]
    conda: "ldsc"
    resources:
        mem_gb = 50,
        walltime = 60
    shell:
        "python /home/jrocha/BioAdHoc/tools/ldsc/ldsc.py --h2 {input.sumstats} --ref-ld-chr {params.ref_ld} --w-ld-chr {params.w_ld} --overlap-annot --frqfile-chr {params.frqfile} --out {params.prefix}"
