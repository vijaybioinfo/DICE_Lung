import pandas as pd
from R24.GWAS.module.Sanity import PROTOCOL_CONSTRUCTOR


def main(args):
    print(args.summary_stats)
    if args.summary_stats:
        PROTOCOL = PROTOCOL_CONSTRUCTOR.construct_sumstats_protocol(args.db)
    else:
        PROTOCOL = PROTOCOL_CONSTRUCTOR.construct_catalog_protocol(args.db)

    data = pd.read_table(args.file)
    data = PROTOCOL.run_protocol(data)

    data.to_csv(args.out, sep="\t", index=False)


def parse_snake_args(snake):
    class args:
        file = snake.input.file
        db = snake.input.db
        summary_stats = snake.params.summary_stats
        out = snake.output[0]

    return args


if __name__ == "__main__":
    main(parse_snake_args(snakemake))