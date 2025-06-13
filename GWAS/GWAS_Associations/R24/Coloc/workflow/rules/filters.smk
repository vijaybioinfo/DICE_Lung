### Filter raw data
#rule filter:
#	input:
#		rules.coloc.output.summary
#	output:
#		OUT + "/filtered/{project2}/{project1}/{dataset1}/{dataset2}.tsv"
#	run:
#		import pandas as pd
#		data = pd.read_table(input[0])
#		data = data[(data["PPH4"] >= 0.5) & (data["PPH4/PPH3"] >= 5)]
#		data.to_csv(output[0], sep="\t", index=False)


rule filter:
	input:
		rules.coloc.output.summary
	output:
		OUT + "/filtered/{project2}/{project1}/{dataset1}/{dataset2}.tsv"
	threads: 1
	resources:
		mem_gb = 1
	shell:
		"awk -F '\t' '$8>=0.5 && $9>=5' {input} > {output}"