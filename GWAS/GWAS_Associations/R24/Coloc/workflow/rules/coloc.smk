def get_coloc_run_args(wildcards):
	args = {**make_run_args(METHOD, wildcards, T1, "1"), **make_run_args(METHOD, wildcards, T2, "2")}
	args["--out"] = rules.coloc.output.summary.format(**wildcards)[:-4]
	list_args = [x for x in args.items() if x[1]]
	return [y for x in list_args for y in x]


### Coloc
rule coloc:
	input:
		trait1 = T1.file,
		trait2 = T2.file
	output:
		summary = OUT + "/raw/{project2}/{project1}/{dataset1}/{dataset2}.tsv",
		snps = OUT + "/raw/{project2}/{project1}/{dataset1}/{dataset2}_snps.tsv"
	params:
		METHOD,
		get_coloc_run_args
	resources:
		mem_gb = lambda wildcards, attempt: attempt * 80,
		walltime = 120
	conda: "coloc"
	script:
		"../scripts/Coloc.R"