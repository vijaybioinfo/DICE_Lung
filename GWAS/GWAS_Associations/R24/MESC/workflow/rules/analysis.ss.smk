##### Summary statistics #####
include: "munge_stats.smk"
include: "expression_scores.ss.smk"
include: "h2med.ss.smk"
include: "summary.ss.smk"

wildcard_constraints:
    disease = "|".join(UDISEASES),
    dataset = "|".join(UDATASETS),
    project = "|".join(IP.Projects(True)),
    sample = "|".join(IP.Datasets(True)),
    category = "|".join(GS.Categories(True))


#localrules: filter_summary_stats

############## Target Files #####################
TARGET_FILES = expand(expand(OUT +"/summary/{{project}}/{disease}/{disease}_{category}.all.h2med.ss", zip, disease=DISEASES, category=CATEGORIES), project=IP.Projects(True)),
        