from R24.Overlap.module import Overlap
from R24.Overlap.module.utils import most_significant, read_file
from R24.Overlap.module.formatting import QTL_FORMATTER, QTL_TYPES

def main(args):
    qtl_fmt, overlap =  QTL_FORMATTER().get(args.qtl), Overlap()

    qtl = qtl_fmt.process(read_file(args.qtl_file))
    gwas = read_file(args.gwas)

    results = overlap(qtl, gwas)
    top = most_significant(results)

    results.to_csv(args.out[0] + ".tsv", sep="\t", index=False)
    top.to_csv(args.out[1] + ".tsv", sep="\t", index=False)


if __name__ == "__main__":
    if not "snakemake" in globals():
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--gwas", required=True, help="GWAS summary statistics file")
        parser.add_argument("--qtl_file", required=True, help="QTL summary statistics filtered file")
        parser.add_argument("--qtl", default="eQTL", choices=list(QTL_TYPES), help="QTL type")
        parser.add_argument("--out", required=True, help="Output prefix")
        #parser.add_argument("--threshold", default=5*10**-8, type=float)
        #parser.add_argument("--genome", choices=["hg38", "hg19"])

        args = parser.parse_args()
        args.out = [args.out, f"{args.out}_most_significant"]
    
    else:
        from collections import namedtuple

        args = namedtuple("args", ["ld", "qtl_file", "out", "qtl"])
        args.gwas = snakemake.input.ld                                                        #type: ignore
        args.qtl_file = snakemake.input.qtl_file                                              #type: ignore  
        args.qtl = snakemake.params.qtl                                                       #type: ignore
        args.out = [snakemake.params.raw, snakemake.params.filtered]                          #type: ignore

    main(args)
