from pandas import DataFrame
from .chain_map import ChainFiles
from .LiftOver import LiftOver as liftover_executor

def liftover(data: DataFrame, source: str, target: str, chains: list[str]) -> DataFrame:
    if source == target:
        print("Nothing to be done...")
        return data

    run_liftover = liftover_executor(ChainFiles.from_names(chains))
    return run_liftover(data, source, target)