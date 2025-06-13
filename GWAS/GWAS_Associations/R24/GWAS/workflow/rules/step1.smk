checkpoint step1:
    input:
        files = GS.file,
        databases = CFG.snp_database.values()
    output:
        directory("step1/{disease}/{dataset}/")
    resources:
        mem_gb = lambda wildcards, attempt: attempt * 50,
        walltime = 180
    params:
        sumstats = GS.sumstats,
        versions = GS.data_genome_versions,
        databases = CFG.snp_database
    script:
        "../scripts/step1.py"