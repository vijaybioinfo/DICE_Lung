def get_file_for_liftover(wildcards):
    folder = checkpoints.step1.get(disease = wildcards.disease, 
                                     dataset = wildcards.dataset
                                    ).output[0]
    
    alt_versions = CFG.get_alternative_versions(wildcards)
    
    for alt in alt_versions:
        if file_exists(f"{folder}/{alt}.tsv"):
            return f"{folder}/{alt}.tsv"

    raise Exception(f"Cannot liftover data ({alt_versions}) to version {wildcards.version}, not chain-map found")

def get_source_version(wildcards):
    file = get_file_for_liftover(wildcards)
    return get_filename(file, ext=False)


rule step2:
    input:
        file = get_file_for_liftover
    output:
        "step2/{version}/{disease}/{dataset}.tsv"
    resources:
        mem_gb = lambda wildcards, attempt: attempt * 50,
        walltime = 60
    params:
        chains = CFG.chain_maps,
        source = get_source_version,
        target = "{version}"
    script:
        "../scripts/step2.py"