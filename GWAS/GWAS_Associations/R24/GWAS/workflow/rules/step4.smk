################################ Functions #################################
def get_datasets(wildcards):
    return expand(rules.step3.output, 
                    version = wildcards.version, 
                    disease = wildcards.disease, 
                    dataset = GS.get_datasets(wildcards)
                )


def get_original_genome_versions(wildcards):
    original_version = {}

    ### Get disease datasets ###
    for dataset_name in GS.get_datasets(wildcards):
        ### Temporay new wildcards ###
        class new_wildcards:
            disease = wildcards.disease
            dataset = dataset_name

        ### Get original genome versions ###
        Gv = GS.data_genome_versions(new_wildcards)

        ### Check if target version exists in original versions ###
        if wildcards.version in Gv:
            original_version[dataset_name] = wildcards.version
            continue
    
        ### Check if versions available for liftover match original versions ###
        for alt in CFG.get_alternative_versions(wildcards):
            if alt in Gv:
                original_version[dataset_name] = alt
                break

    return original_version
        

################################ Rules #################################
rule risk_alleles:
    input:
        get_datasets
    output:
        CFG.OutDir + "/risk_alleles/{version}/{disease}.tsv"
    resources:
        mem_gb = 30,
        walltime = 300
    script:
        "../scripts/step4.py"


rule metadata:
    input:
        files = get_datasets,
        sources = GS.get_source_files,
        sample_table = GS.sample_table_file
    output:
        CFG.OutDir + "/sample_tables/{version}/{disease}.csv"
    resources:
        mem_gb = 20,
        walltime = 180
    params:
        disease = "{disease}",
        version = "{version}",
        original_version = get_original_genome_versions
    script:
        "../scripts/metadata.py"