##########################################
############### Functions ################
##########################################
def get_file_for_step3(wildcards):
    folder = checkpoints.step1.get(disease = wildcards.disease, 
                                     dataset = wildcards.dataset
                                    ).output[0]

    ##### Take original file if genome version exists 
    step1_file = f"{folder}/{wildcards.version}.tsv"
    if file_exists(step1_file): return step1_file

    ##### otherwise take liftover file
    return f"step2/{wildcards.version}/{wildcards.disease}/{wildcards.dataset}.tsv"


def get_db(wildcards): return CFG.snp_database[wildcards.version]


def is_sumstats(wildcards):
    stats = GS.sumstats(wildcards)

    if len(set(stats)) > 1:
        raise Exception(f"Sample {wildcards.dataset} from Disease {wildcards.disease} has more than one sumstats value. It should be either 1 or 0")

    return bool(stats[0])



######################################
############### Rules ################
######################################

rule step3:
    input: 
        file = get_file_for_step3,
        db = get_db
    output:
        CFG.OutDir + "/data/{version}/{disease}/{dataset}.tsv"
    resources:
        mem_gb = lambda wildcards, attempt: get_mem_gb(get_file_for_step3(wildcards), attempt),
        walltime = 180
    params:
        summary_stats = is_sumstats
    script:
        "../scripts/step3.py"