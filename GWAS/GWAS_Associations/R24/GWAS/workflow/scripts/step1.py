from R24.GWAS.module.utils.misc import mkdir
from R24.GWAS.module.Parsers import PARSER_CONSTRUCTOR
from R24.GWAS.module.Sanity.programs import COORDINATES_PROGRAM
from R24.GWAS.module.Sanity.Coordinates import make_coordinates_sanitizer


def main(args):
    
    for file, stats, version in zip(args.files, args.sumstats, args.versions):

        if int(stats):
            PARSER = PARSER_CONSTRUCTOR.construct_sumstats_parser()
            SANITIZER = make_coordinates_sanitizer()
        else:
            PARSER = PARSER_CONSTRUCTOR.construct_catalog_parser()
            SANITIZER = make_coordinates_sanitizer(args.dbs[version])

        SANITIZER.PROGRAM.add(COORDINATES_PROGRAM)

        data = PARSER.parse(file)
        data = SANITIZER.clean(data, dropna=True)

        mkdir(args.out)
        data.to_csv(f"{args.out}/{version}.tsv", sep="\t", index=False)


def parse_snake_args(snake):
    class args:
        files = snake.input.files
        dbs = snake.params.databases
        sumstats = snake.params.sumstats
        versions = snake.params.versions
        out = snake.output[0]

    return args

if __name__ == "__main__":
    main(parse_snake_args(snakemake))



