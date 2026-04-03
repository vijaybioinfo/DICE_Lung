import numpy as np
import pandas as pd

from .base import SANITIZER

from ..utils.columns import Columns
from ..utils.snpdb import ReferenceSNP

from .Coordinates import make_coordinates_sanitizer
from .Statistics.sanitizers import STATS_BUILDER, NULL_STATS_BUILDER
from .programs import COORDINATES_PROGRAM, ALLELES_PROGRAM, STATS_PROGRAM
from .Alleles.sanitizers import BIALLELE_SANITIZER_BUILDER, MULTIALLELE_SANITIZER_BUILDER

##########################################################################
############################## COLUMNS ################################### 
##########################################################################
INFO = [Columns.CHR, Columns.POS, Columns.NEA, Columns.EA, Columns.SNPID, Columns.RSID, "Original_RSID"]
STATS = [Columns.PVAL, Columns.BETA, Columns.SE, Columns.Z, Columns.AF]
EXTRA = [Columns.RA, Columns.NRA, "Notes", "found"]

##########################################################################
########################## HELPER FUNCTIONS ##############################
##########################################################################

################### STATS ###################
def adjust_sign(where: "pd.Series[bool]", values: "pd.Series[float]") -> "pd.Series[float]":
    return np.where(where, values * -1, values)                                                               #type: ignore

def adjust_AF(where: "pd.Series[bool]", values: "pd.Series[float]") -> "pd.Series[float]":
    return np.where(where, 1 - values, values)                                                                #type: ignore

################## ALLELES ###################
def _redefine_alleles(frame: pd.DataFrame) -> pd.DataFrame:
    where = frame["Notes"] == "Not Found"

    found, nofound = frame[~where], frame[where]
    found = found.drop(columns=[Columns.NEA, Columns.EA])
    nofound = nofound.drop(columns=[Columns.REF, Columns.ALT])

    found = found.rename(columns={Columns.REF: Columns.NEA, Columns.ALT: Columns.EA})       #type: ignore
    return pd.concat([found, nofound])

def get_risk_allele(data: pd.DataFrame) -> pd.DataFrame:
    if (not Columns.BETA in data.columns) and (not Columns.Z in data.columns):
        print("Warning: Missing BETA and Z so risk allele cant be identified")
        data[Columns.RA] = np.nan
        data[Columns.NRA] = np.nan
        return data

    signed_column = Columns.BETA if Columns.BETA in data.columns else Columns.Z
    where = np.array(data[signed_column]) >= 0

    data[Columns.RA] = np.where(where, data[Columns.EA], data[Columns.NEA])
    data[Columns.NRA] = np.where(where, data[Columns.NEA], data[Columns.EA])
    return data

################## IDs #######################
def add_SNPID(data: pd.DataFrame) -> pd.DataFrame:
    #In this case NEA = REF and EA = ALT
    data[Columns.SNPID] = data[Columns.CHR].astype(str).str[:] + ":" + data[Columns.POS].astype(str).str[:] + ":" + \
                          data[Columns.NEA].str[:] + ":" + data[Columns.EA].str[:]   
    
    return data

def add_RSID(data: pd.DataFrame, db: str):
    def _get(x):
        frame = x[1]
        database = ReferenceSNP(db)
        
        chromosome, position = x[0], frame[Columns.POS].astype(int).tolist()
        res = database.search_chromosome_position(chromosome, position, [Columns.SNPID, Columns.RSID])
        res.drop_duplicates(inplace=True)

        return pd.merge(frame, res, how="left", on=Columns.SNPID)

    if Columns.RSID in data.columns:
        data.rename(columns={Columns.RSID: "Original_RSID"}, inplace=True)
    else:
        data["Original_RSID"] = np.nan

    results = pd.concat(map(_get, data.groupby(Columns.CHR)))
    results = pd.concat([results, data[data[Columns.CHR].isna()]])

    return results.drop_duplicates(results.columns.difference(["Original_RSID", Columns.RSID]))

######################### MISC ###############################
def fill_if_empty(cols: list[str], data: pd.DataFrame) -> pd.DataFrame:
    for col in cols:
        if not col in data.columns:
            data[col] = np.nan

    return data

def check_mandatory_columns(mandatory: list[str], data: pd.DataFrame):
    missing = set(mandatory) - set(data.columns)
    if missing:
        raise Exception(f"The following columns are missing: {list(missing)}")
        
def order_columns(data: pd.DataFrame) -> pd.DataFrame:
    order = [col for col in INFO + STATS + EXTRA if col in data.columns]
    return data[order]

##################################################################
######################### PROTOCOL ###############################
##################################################################

class Protocol:
    COORDINATES: SANITIZER
    ALLELES: SANITIZER
    STATS: SANITIZER

    def is_ready(self):
        C = hasattr(self, "COORDINATES") and not self.COORDINATES is None
        P = hasattr(self, "ALLELES") and not self.ALLELES is None
        S = hasattr(self, "STATS") and not self.STATS is None

        return all([C, P, S])

    def _polishing(self, data: pd.DataFrame) -> pd.DataFrame:
        data = data[~data[Columns.PVAL].isna()]

        if data.empty:
            return pd.DataFrame(columns = INFO + STATS + EXTRA)

        data = add_SNPID(data)
        data = add_RSID(data, self.ALLELES.RECOVERY_METHOD.MATCHER.PANEL)                         #type: ignore

        where = data["Notes"].isin(["Flip", "Flip-Complement", "Flip-Fill", "Flip-Complement-Fill"])

        if Columns.BETA in data.columns:
            data[Columns.BETA] = adjust_sign(where, data[Columns.BETA])

        if Columns.Z in data.columns:
            data[Columns.Z] = adjust_sign(where, data[Columns.Z])
        
        if Columns.AF in data.columns:
            data[Columns.AF] = adjust_AF(where, data[Columns.AF])

        data = get_risk_allele(data)
        
        data.sort_values(Columns.PVAL, inplace=True)
        data.drop_duplicates(data.columns.difference([Columns.PVAL]), inplace=True)
        
        data.sort_values([Columns.CHR, Columns.POS], inplace=True)

        return order_columns(data).reset_index(drop=True)

    def run_protocol(self, data: pd.DataFrame) -> pd.DataFrame:
        if not self.is_ready():
            raise Exception("Protocol is not ready, use PROTOCOL_CONSTRUCTOR to properly initialize this class...")

        data = fill_if_empty([Columns.NEA, Columns.EA], data)

        ########## Check Required Columns ###########
        MANDATORY = [Columns.CHR, Columns.POS, Columns.NEA, Columns.EA, Columns.PVAL]
        check_mandatory_columns(MANDATORY, data)

        ########## Programs ###########
        CP = COORDINATES_PROGRAM
        AP = ALLELES_PROGRAM
        SP = [item for item in data.columns if item in STATS_PROGRAM]

        ########## Add Programs ############
        self.COORDINATES.PROGRAM.add(CP)
        self.ALLELES.PROGRAM.add(AP)
        self.STATS.PROGRAM.add(SP)

        ########## Run Programs ############
        data = self.COORDINATES.clean(data, dropna=True)
        data = self.ALLELES.clean(data, dropna=True)
        data = self.STATS.clean(data, dropna=False)

        ######### Final Polishing ########
        data = self._polishing(data)

        return data

##########################################################################
########################## PROTOCOL CONSTRUCTOR ##########################
##########################################################################
def make_alleles_sanitizer(builder, panel: str):
    BUILDER = builder()
    BUILDER.set_wiper_factory()
    BUILDER.set_recovery_method(panel)

    return BUILDER.get_sanitizer()

def make_stats_sanitizer(builder):
    BUILDER = builder()
    BUILDER.set_wiper_factory()
    BUILDER.set_recovery_method()

    return BUILDER.get_sanitizer()

def make_protocol(C, A, S) -> Protocol:
    protocol = Protocol()

    protocol.COORDINATES = C
    protocol.ALLELES = A
    protocol.STATS = S

    return protocol

class PROTOCOL_CONSTRUCTOR:
    @staticmethod
    def construct_catalog_protocol(panel: str):
        ############# BUILD COORDINATES SANITIZER ######################
        coordinates = make_coordinates_sanitizer(panel)

        ############# BUILD ALLELE SANITIZER #################
        alleles = make_alleles_sanitizer(MULTIALLELE_SANITIZER_BUILDER, panel)

        ############# BUILD STATS SANITIZER ##################
        stats = make_stats_sanitizer(NULL_STATS_BUILDER)

        ############# CONSTRUCT PROTOCOL ##################
        return make_protocol(coordinates, alleles, stats)

    @staticmethod
    def construct_sumstats_protocol(panel: str):
        ############# BUILD COORDINATES SANITIZER ######################
        coordinates = make_coordinates_sanitizer()

        ############# BUILD ALLELE SANITIZER #################
        alleles = make_alleles_sanitizer(BIALLELE_SANITIZER_BUILDER, panel)

        ############# BUILD STATS SANITIZER ##################
        stats = make_stats_sanitizer(STATS_BUILDER)

        ############# CONSTRUCT PROTOCOL ##################
        return make_protocol(coordinates, alleles, stats)
