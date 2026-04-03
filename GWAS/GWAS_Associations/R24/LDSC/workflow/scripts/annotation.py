from collections import namedtuple
from R24.LDSC.ldsc.annotation import main

###################### Parse Arguments ########################
args = namedtuple("args", "bim annot bed labels out chr start end")

if snakemake.params.type:
    args.bim = None
    args.annot = snakemake.input.variants
else:
    args.bim = snakemake.input.variants
    args.annot = None

args.bed = snakemake.input.bed

args.out = snakemake.params.prefix
args.labels = None

args.chr = None
args.start = None
args.end = None


##################### Annotate ###########################
main(args)