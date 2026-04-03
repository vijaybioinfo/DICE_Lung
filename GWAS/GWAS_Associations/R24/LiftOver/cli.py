from R24.LiftOver.utils import read_file
from R24.LiftOver import liftover

def parse_args():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--file", required=True, help="Input file")
    parser.add_argument("--source", required=True, help="Original genome version")
    parser.add_argument("--target", required=True, help="Output genome version")
    parser.add_argument("--chains", required=True, type=str, nargs="+", help="Chain file(s)")
    parser.add_argument("--out", required=True, help="Output prefix")

    args = parser.parse_args()

    return args

def liftover_from_args(args):
    #executor = liftover(args.chains)
    data = read_file(args.file)
    
    return liftover(data, args.source, args.target, args.chains)


def liftover_cli(args = None):
    if not args:
        args = parse_args()
        
    results = liftover_from_args(args)
    results.to_csv(f"{args.out}.tsv", index=False, sep="\t")

if __name__ == "__main__":
    liftover_cli()