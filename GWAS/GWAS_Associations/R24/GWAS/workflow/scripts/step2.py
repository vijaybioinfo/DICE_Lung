import numpy as np
import pandas as pd

from R24.LiftOver import liftover
from R24.LiftOver.utils import read_file
from R24.GWAS.module.utils.columns import Columns


#########################
##### Main LiftOver #####
#########################
def main(args):
    if args.source == args.target:
        raise Exception(f"Trying to liftover from version {args.source} to {args.target}")
    
    results = liftover(read_file(args.file), args.source, args.target, args.chains)                        #type: ignore

    isna = results[[Columns.CHR, Columns.POS]].isna().any(axis=1)
    results, nans = results[~isna], results[isna]

    results[Columns.SNPID] = results[Columns.CHR].astype(str).str.replace("\.0", "").str[:] + ":" + \
                        results[Columns.POS].astype(str).str.replace("\.0", "").str[:] + ":" + \
                        results[Columns.NEA].fillna("").str[:] + ":" + \
                        results[Columns.EA].fillna("").str[:]

    nans[Columns.SNPID] = np.nan
    results = pd.concat([results, nans])
    
    results.to_csv(args.out, index=False, sep="\t")


#########################
#### Parse Arguments ####
#########################
def parse_snake_args(snake):
    class args:
        file = snake.input.file                                                                          #type: ignore
        out = snake.output[0]                                                                            #type: ignore
        chains = snake.params.chains                                                                     #type: ignore
        source = snake.params.source                                                                     #type: ignore
        target = snake.params.target                                                                     #type: ignore

    return args


if __name__ == "__main__":
    main(parse_snake_args(snakemake))