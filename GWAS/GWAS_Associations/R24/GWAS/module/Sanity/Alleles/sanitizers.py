import pandas as pd

from ...utils.columns import Columns
from .wipers import AlleleWiperFactory
from .matchers import PANEL_MATCHER, PANEL_DATABASE_MATCH
from ..base import RECOVERY_METHOD, SANITIZER_BUILDER, RECOVER_METHOD_NOT_READY
from .harmonizers import ALLELE_HARMONIZER, BIALLELIC_HARMONIZER, MULTI_ALLELE_HARMONIZER


##################################################################################
############################## RECOVERY METHOD ###################################
##################################################################################

class ALLELE_RECOVERY_METHOD(RECOVERY_METHOD):
    MATCHER: PANEL_MATCHER
    HARMONIZER: ALLELE_HARMONIZER

    def check(self):
        if not(hasattr(self, "MATCHER") and hasattr(self, "HARMONIZER")):
            raise RECOVER_METHOD_NOT_READY("METHOD not ready, use a builder class to properly initialize...")
        
        #if not isinstance(self.MATCHER, type(PANEL_MATCHER)):
        #    raise ValueError(f"Attribute MATCHER must be of type {type(PANEL_MATCHER)}, but {type(self.MATCHER)} found instead...")
        
        #if not isinstance(self.HARMONIZER, type(ALLELE_HARMONIZER)):
        #    raise ValueError(f"Attribute HARMONIZER must be of type {type(ALLELE_HARMONIZER)}, but {type(self.HARMONIZER)} found instead...")

    def recover(self, data: pd.DataFrame) -> pd.DataFrame:
        ################ Separate NANs from data ##############
        isna = data[Columns.CHR].isna()
        na, ok = data[isna], data[~isna]

        na["Notes"] = "Not Found"
        na["found"] = False

        recovered = []
        ok.rename(columns={Columns.EA: Columns.A1, Columns.NEA: Columns.A2}, inplace=True)

        ############### Recovery Method ###################
        for chr in ok[Columns.CHR].unique():
            print(f"[INFO] - Recovering Chromosome {chr}...")
            tmp = ok[ok[Columns.CHR] == chr]

            print("[INFO] - Matching alleles...")
            found, notfound = self.MATCHER.match(tmp)

            print("[INFO] - Harmonizing alleles...")
            found = self.HARMONIZER.harmonize(found)

            print("[INFO] - Harmonization finished ...")
            notfound["Notes"] = "Not Found"
            notfound["found"] = False

            tmp = pd.concat([found, notfound])
            
            print(f"[INFO] - Finished recovering Chromosome {chr} ...")
            recovered.append(tmp)

        ok = pd.concat(recovered)
        ok.rename(columns={Columns.A2: Columns.NEA, Columns.A1: Columns.EA}, inplace=True)

        return pd.concat([ok, na])

##################################################################################################
########################################### BUILDERS #############################################
##################################################################################################

class BIALLELE_SANITIZER_BUILDER(SANITIZER_BUILDER):
    def set_wiper_factory(self):
        self.product.FACTORY = AlleleWiperFactory()

    def set_recovery_method(self, panel: str):
        recovery = ALLELE_RECOVERY_METHOD()
        
        recovery.MATCHER = PANEL_DATABASE_MATCH()
        recovery.MATCHER.add_panel(panel)

        recovery.HARMONIZER = BIALLELIC_HARMONIZER()

        recovery.check()
        self.product.RECOVERY_METHOD = recovery


class MULTIALLELE_SANITIZER_BUILDER(SANITIZER_BUILDER):
    def set_wiper_factory(self):
        self.product.FACTORY = AlleleWiperFactory()

    def set_recovery_method(self, panel: str):
        recovery = ALLELE_RECOVERY_METHOD()
        
        recovery.MATCHER = PANEL_DATABASE_MATCH()
        recovery.MATCHER.add_panel(panel)

        recovery.HARMONIZER = MULTI_ALLELE_HARMONIZER()

        recovery.check()
        self.product.RECOVERY_METHOD = recovery

