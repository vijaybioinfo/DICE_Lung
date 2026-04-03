import numpy as np

from R24.LD import LDFactory
from R24.Overlap.module.formatting import GWAS
from R24.Overlap.module.utils import Columns, read_file

def main(args):
    if not args.LD is None:
        ld_engine = LDFactory.build("PLINK")
        ld_engine.initialize(args.LD)
        
        fmt = GWAS(ld_engine, 5*10**-8)
    else:
        fmt = GWAS(threshold = 5*10**-8)
        
    gwas = read_file(args.gwas)

    if not Columns.BETA in gwas.columns:
        gwas[Columns.BETA] = np.nan
    
    ld = fmt.process(gwas[[Columns.SNPID, Columns.PVAL, Columns.BETA]])
    ld.to_csv(f"{args.out}.ld", index=False, sep="\t")


if __name__ == "__main__":
    from collections import namedtuple

    args = namedtuple("args", ["gwas", "out", "LD"])
    args.gwas = snakemake.input.gwas                                                      #type: ignore
    args.LD = snakemake.params.LD                                                         #type: ignore
    args.out = snakemake.params.prefix                                                    #type: ignore

    main(args)
