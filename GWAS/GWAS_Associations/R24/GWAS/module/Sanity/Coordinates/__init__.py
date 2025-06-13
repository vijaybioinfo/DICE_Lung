import re
import pandas as pd

from numpy import nan
from typing import Union
from ...utils.columns import Columns
from ...utils.snpdb import ReferenceSNP
from .sanitizers import make_coordinates_sanitizer

################################################################
##################### Check Genome Version #####################
################################################################

###################### Helper Functions ##########################
def _has(data: pd.DataFrame, cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    '''return 2 dataframes one containing rows with missing data for columns specified, and the otherone with complete data for columns'''                            
    where = data[cols].isna().any(axis=1)
    return data[~where], data[where] 

def has_rsid(data: pd.DataFrame): return _has(data, [Columns.RSID])

def has_missing_coordinates(data: pd.DataFrame): return _has(data, [Columns.CHR, Columns.POS])[::-1] 

###Match
def _match(data: pd.DataFrame, db: str):
    def _look(x):
        reference = ReferenceSNP(db)
        query = reference.search_chromosome_position(x[0], x[1][Columns.POS].astype(int).tolist(), [Columns.CHR, Columns.POS])
        return query.drop_duplicates()

    data[Columns.POS] = data[Columns.POS].astype(int)
    query = pd.concat(map(_look, data.groupby(Columns.CHR)))
    match_ = pd.merge(data, query, how="left", on=[Columns.CHR, Columns.POS], indicator=True)
    match_["_merge"] = match_["_merge"].replace({"both": True, "left_only": False}).astype(bool)
    
    return match_

def match_versions(data: pd.DataFrame, dbs: dict[str, str]):
    versions = data.copy(deep=True)
    for name, db in dbs.items():
        tmp =  _match(data, db)[[Columns.CHR, Columns.POS, "_merge"]].rename(columns={"_merge": name})
        versions = pd.merge(versions, tmp)
        
    return versions

def select_best_version_match(data: pd.DataFrame, versions: list[str]):
    name, best, more = "", -1, []
    
    for v in versions:
        fraction = data[v].astype(int).sum() / data.shape[0]
        
        if fraction > best:
            name, best, more = v, fraction, [v]
        elif fraction == best:
            more.append(v)

    if best < 1:
        print("INFO - Not perfect match found for snps...")
        print("INFO - Assigning best match possible...")

        if best < 0.9:
            print(f"WARNING - Best proportion of snps matching is {best}")

    core_cols = data.columns.difference(versions)
    
    if more:
        mx, name = -1, ""
        for label in more:
            if (num := re.search("\d+", label)):
                num = int(num.group())
            else:
                raise Exception(f"Cannot parse numbers from genome version {label}")
            
            if num > mx:
                mx, name = num, label
            
    return data[list(core_cols) + [name]]

###check version
def find_genome_version(data: pd.DataFrame, dbs: dict[str, str], expected: Union[str, None] = None):
    cols = data.columns.tolist()

    ############ Match and select best #################
    matches = match_versions(data, dbs)
    best = select_best_version_match(matches, list(dbs.keys()))

    ############ Check only one best #############
    version = list(best.columns.difference(cols))
    if len(version) != 1:
        raise ValueError(f"Number of best versions should be 1, but {len(version)} found...")
    version = version[0]
    
    ############ Check if best equals expected ############
    if not expected is None and expected != version:
        #raise Exception(f"ERROR - You are expecting version {expected} but best match is {version}")
        print(f"WARNING - You are expecting version {expected} but best match is {version}")
    
    ############ Check for low match #############
    perc = best[version].sum() / best.shape[0]
    if perc < 0.9 and best.shape[0] >= 20:
        raise Exception(f"ERROR - Proportion of snps matching to selected genome version is too low ({perc}). Cannot identify genome version")
    
    return best

def _core(data: pd.DataFrame, dbs: dict[str, str], expected: Union[str, None] = None):
    if data.empty:
        return pd.DataFrame(columns=list(data.columns) + ["Genome"])

    miss, nomiss = has_missing_coordinates(data)

    if nomiss.empty:
        data["Genome"] = nan
        return data

    if not miss.empty and (nomiss.shape[0] / miss.shape[0]) < 0.5 and data.shape[0] > 20:
        print("WARNING - SNPs with coordinates are less than snps with no coordinates. Cannot confidently find genome version.")
        data["Genome"] = nan
        return data

    nomiss[Columns.CHR] = nomiss[Columns.CHR].astype(str)
    nomiss[Columns.POS] = nomiss[Columns.POS].astype(int)
    nomiss = find_genome_version(nomiss, dbs, expected)
    version = nomiss.columns.difference(data.columns.tolist())[0]
    
    data["Genome"] = version
    return data

def get_gwas_catalog_genome_version(data: pd.DataFrame, dbs: dict[str, str], expected: str):
    withRSID, noRSID = has_rsid(data)

    ###SNPs with RSID
    withRSID = _core(withRSID, dbs, expected)

    ###SNPs with no RSID
    noRSID = _core(noRSID, dbs)

    return pd.concat([withRSID, noRSID])

###############################################################################################################################