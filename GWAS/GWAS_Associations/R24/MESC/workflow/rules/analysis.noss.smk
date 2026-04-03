##### Non Summary statistics #####
include: "munge_stats.smk"
include: "expression_scores.noss.smk"
include: "h2med.noss.smk"
include: "summary.noss.smk"

wildcard_constraints:
    disease = "|".join(UDISEASES),
    dataset = "|".join(UDATASETS),
    project = "|".join(IP.Projects(True)),
    sample = "|".join(IP.Datasets(True)),
    category = "|".join(GS.Categories(True)),
    group = "|".join(IP.Groups(True))



localrules: filter_group, filter_non_summary_stats

############## Target Files #####################
TARGET_FILES = expand(expand(OUT + "/summary/{{project}}/{disease}/{disease}_{category}.all.h2med.noss", zip, disease=DISEASES, category=CATEGORIES), project=IP.Projects(True))
TARGET_FILES = TARGET_FILES + expand(expand(OUT + "/summary/{{project}}/{disease}/{disease}_{category}.all.h2med.groups", zip, disease=DISEASES, category=CATEGORIES), project=IP.Projects(True))