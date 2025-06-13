from pandas import DataFrame
from .wipers import WiperFactory
from .methods import MISSING_REQUIREMENTS, MethodFactory
from ..base import RECOVER_METHOD_NOT_READY, SANITIZER_BUILDER, RECOVERY_METHOD

class STATS_RECOVERY_METHOD(RECOVERY_METHOD):
    METHOD_FACTORY = MethodFactory()

    def check(self):
        if not hasattr(self, "METHOD_FACTORY"):
            raise RECOVER_METHOD_NOT_READY(f"Recovery Method is not ready, please use class builder to properly initialize...")

    def recover(self, data: DataFrame) -> DataFrame:
        options = list(self.METHOD_FACTORY.METHODS.keys())
        possible = [item for item in options if not item in data.columns]
        print(possible)

        for item in possible + possible:
            method = self.METHOD_FACTORY.new(item)
            
            try:
                if not (tmp := method(data)) is None:
                    data = tmp
            except MISSING_REQUIREMENTS:
                continue

        return data
    
class NULL_RECOVERY_METHOD(RECOVERY_METHOD):
    def check(self): pass

    def recover(self, data: DataFrame) -> DataFrame: return data
    
####################################################################################
################################## BUILDER #########################################
####################################################################################

class STATS_BUILDER(SANITIZER_BUILDER):
    def set_wiper_factory(self):
        self.product.FACTORY = WiperFactory()

    def set_recovery_method(self):
        self.product.RECOVERY_METHOD = STATS_RECOVERY_METHOD()


class NULL_STATS_BUILDER(SANITIZER_BUILDER):
    def set_wiper_factory(self):
        self.product.FACTORY = WiperFactory()

    def set_recovery_method(self):
        self.product.RECOVERY_METHOD = NULL_RECOVERY_METHOD()

        
