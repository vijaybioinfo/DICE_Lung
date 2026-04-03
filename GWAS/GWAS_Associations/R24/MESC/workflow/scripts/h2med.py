import subprocess


def make_empty_files(prefix):
    all_h2med = "Quantity\tEstimate\tSE(Estimate)\tEstimate_over_h2\tSE(Estimate_over_h2)"
    categories_h2med = "Gene_category\tNum_genes\tProp_genes\tAve_h2cis\th2med\tSE(h2med)\th2med_over_h2med(tot)\tSE(h2med_over_h2med(tot))\th2med_enrichment\tSE(h2med_enrichment)\th2med_enrichment_pvalue\tCoefficient\tSE(Coefficient)"

    with open(f"{prefix}.all.h2med", "w") as file:
        print(all_h2med, file=file)

    with open(f"{prefix}.categories.h2med", "w") as file:
        print(categories_h2med, file=file)
        

def main(snakemake):
    input = snakemake.input
    params = snakemake.params

    exe = subprocess.run(f"{params.software} --h2med {input.gwas} --exp-chr {params.exp} --out {params.out}", shell=True, stderr=subprocess.PIPE)

    if exe.returncode:
        messages = [
                    "ValueError: After merging GWAS summary statistics with LD scores, 0 SNPs remain.", 
                    "ValueError: More dimensions than datapoints.",
                    "ValueError: More blocks than data points."
                ]
        stderr = exe.stderr.decode("utf-8")

        for message in messages:
            if message in stderr:
                make_empty_files(params.out)
                break
        else:
            raise Exception(stderr)


if __name__ == "__main__":
    main(snakemake)


