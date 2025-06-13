################ Input Functions ###################
def dataset_merge(wildcards):
	datasets = T2.get_datasets(wildcards.project2)
	return [OUT + "/filtered/{{project2}}/{{project1}}/{{dataset1}}/{dataset2}.tsv".format(dataset2=file) for file in datasets]

def category_merge(wildcards):
	datasets = T1.get_by_category(wildcards)
	return [OUT + "/filtered/{{project2}}/{{project1}}/{dataset1}.tsv".format(dataset1=file) for file in datasets]


####################################### Rules ##################################

### Combine filtered data
rule merge:
	input:
		dataset_merge
	output:
		OUT + "/filtered/{project2}/{project1}/{dataset1}.tsv"
	params:
		trait = "dataset2"
	run:
		from R24.Coloc.module.utils import merge_coloc_files
		merged = merge_coloc_files(input, params.trait)
		merged.to_csv(output[0], sep="\t", index=False)


### summarize 
rule merge_by_category:
	input:
		category_merge
	output:
		table = OUT + "/summary/{project2}/{project1}/{project1}_{category}.tsv",
		matrix = OUT + "/summary/{project2}/{project1}/{project1}_{category}.mtx"
	params:
		column = "dataset1",
		metadata = T2.metadata
	script:
		"../scripts/summary.py"