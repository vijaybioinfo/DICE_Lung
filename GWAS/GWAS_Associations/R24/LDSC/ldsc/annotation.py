import os
import gzip
import time
import numpy as np


GENETIC_VARIANTS = list[tuple[str, int, str, str]]
ANNOTATION = tuple[int, ...]
COORDINATES = list[tuple[str, int, int]]

#################### Helper Functions #####################
def _get_bed_header(file: str, chr: str = None, start: str = None, end: str = None):
    exists = ((not chr is None), (not start is None), (not end is None))

    if all(exists):
        with open(file) as stream:
            row = stream.readline().strip().split()
            
            if missing := [c for c in (chr, start, end) if not c in row]:
                print(f"Following specified columns are missing: {', '.join(missing)}")

            return (row.index(chr), row.index(start), row.index(end)), True

    if any(exists):
        raise Exception("Columns 'Chromosome', 'Start' and 'End' are partially set, please provide all of them. Alternatively you can unset this options to use three first columns instead.")
    
    return (0, 1, 2), False

def split_by_chr(regions: list) -> dict[str, COORDINATES]:
    chunks = {}

    for region in regions:
        if not region[0] in chunks:
            chunks[region[0]] = []
        
        chunks[region[0]].append(region)
    
    return chunks

#################### Read Functions #################
def read_variants_from_annot(file: str) -> tuple[GENETIC_VARIANTS, dict[str, ANNOTATION]]:
    variants, annot = [], []

    with gzip.open(file, "rt") as stream:
        annot_labels = stream.readline().strip().split()[4:]

        while line := stream.readline().strip().split():
            variants.append( (line[0], int(line[1]), line[2], line[3]) )
            annot.append( line[4:] )

    annot = {l: a for l, a in zip(annot_labels, zip(*annot))}
    return variants, annot

def read_bim(file: str) -> GENETIC_VARIANTS:
    variants = []
    with open(file, "r") as stream:
        while line := stream.readline().strip().split():
            variants.append((line[0], int(line[3]), line[1], line[2]))

    return variants

def read_bed(file: str, chr: str = None, start: str = None, end: str = None) -> COORDINATES:
    regions = []
    H, has_header = _get_bed_header(file, chr, start, end)

    with open(file, "r") as stream:
        if has_header:
            stream.readline()

        while line := stream.readline().strip().split():
            regions.append((line[H[0]], int(line[H[1]]), int(line[H[2]])))

    return regions

################ Intersections #######################
def _base_annotation(variants: GENETIC_VARIANTS):
    return {"base": [1] * len(variants)}

def intersect(variants: GENETIC_VARIANTS, regions: COORDINATES):
    ### Split by Chromosomes ###
    R = {l: (np.array([k[1] for k in i]), np.array([k[2] for k in i])) for l, i in split_by_chr(regions).items()}
    return [ int(any((R[C][0] <= P) & (P <= R[C][1]))) for C, P, _, _ in variants ]


def make_annotation(variants: GENETIC_VARIANTS, annot: dict[str, ANNOTATION], prefix: str):
    headers = ["CHR", "BP", "SNP", "CM"] + list(annot.keys())

    lines = ["\t".join(headers)] + ["\t".join(map(str, V + A)) for V, A in zip(variants, zip(*annot.values()))]
    lines = "\n".join(lines) + "\n"

    with gzip.open(f"{prefix}.annot.gz", "wb") as stream:
        stream.write( lines.encode("utf-8") )


def main(args):
    if args.labels is None:
        args.labels = [os.path.splitext(os.path.basename(bedfile))[0] for bedfile in args.bed]

    if len(args.labels) != len(args.bed):
        raise Exception("Number of annotation labels doesn't match number of bed files...")

    VARIANTS, ANNOT = read_variants_from_annot(args.annot) if args.annot else (read_bim(args.bim), {})

    if conflicting := set(args.labels) & set(ANNOT.keys()):
        raise Exception(f"Following labels are shared between old annotation and new annotations: {', '.join(conflicting)}")

    REGIONS = {label: read_bed(file, args.chr, args.start, args.end) for label, file in zip(args.labels, args.bed)}
    INTERSECTIONS = {label: intersect(VARIANTS, regions) for label, regions in REGIONS.items()}
    
    ANNOT.update(_base_annotation(VARIANTS))
    ANNOT.update(INTERSECTIONS)

    make_annotation(VARIANTS, ANNOT, args.out)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--bim", default=None, help="Genotype file")
    parser.add_argument("--annot", default=None, help="Gzipped file with annotation (use this when you already have an annotation matrix and want to add more).")

    parser.add_argument("--bed", required=True, nargs="+", help="Bed file(s) with genomic coordinates")
    parser.add_argument("--labels", default=None, nargs="+", help="Space separated list of labels to use for each annotation column")
    parser.add_argument("--out", required=True, help="Output prefix")

    parser.add_argument("--chr", default=None, help="Chromosome column name in bed file")
    parser.add_argument("--start", default=None, help="Start column name in bed file")
    parser.add_argument("--end", default=None, help="End column name in bed file")

    args = parser.parse_args()

    
    ########################################################################
    if args.bim is None and args.annot is None:
        raise Exception("You must specify either --bim or --annot")
    
    if not args.bim is None and not args.annot is None:
        raise Exception("You must specify either --bim or --annot not both")

    main(args)