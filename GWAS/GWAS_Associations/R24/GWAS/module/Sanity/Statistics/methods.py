import numpy as np

from pandas import DataFrame
from scipy.stats import norm
from abc import ABC, abstractmethod
from ...utils.columns import Columns as Col

################################################################################
############################### BASE STAT METHOD ###############################
################################################################################
class MISSING_REQUIREMENTS(Exception):
    pass

class BaseMethod(ABC):
    NAME: str
    REQS: list[tuple[str, ...]]

    def _meets_req(self, data: DataFrame) -> bool:
        meets = []
        for option in self.REQS:
            found = all([col in data.columns for col in option])
            meets.append(found)
        return any(meets)

    def _preflight(self, data: DataFrame):
        if not self._meets_req(data):
            raise MISSING_REQUIREMENTS(f"To estimate {self.NAME} following columns are required: {self.REQS[0]}")

    @abstractmethod
    def _method(self, data: DataFrame, **kwargs) -> DataFrame:
        ...
    
    def __call__(self, data: DataFrame, force: bool = False, **kwargs):
        if not self.NAME in data.columns or force:
            self._preflight(data)
            return self._method(data, **kwargs)
        
        print(f"WARNING - Column {self.NAME} already exists. Set force = True to overwrite.")


#########################################################################################
################################## IMPLEMENTATIONS ######################################
#########################################################################################

class EFFECT_SIZE(BaseMethod):
    NAME = Col.BETA
    REQS = [tuple([Col.OR])]

    def _method(self, data: DataFrame) -> DataFrame:
        '''Check for BETA or OR columns and estimate if missing'''
        # Formulas obtained from https://bogdan.dgsom.ucla.edu/news/item?item_id=54624
        # Estimate BETA if missing and OR column is present
        reqs = self.REQS[0][0]
        data[self.NAME] = np.log(data[reqs]) 
        return data


class ZSCORE(BaseMethod):
    NAME = Col.Z
    REQS = [(Col.BETA, Col.SE), (Col.PVAL, Col.BETA)]

    def _method(self, data: DataFrame) -> DataFrame:
        '''Estimate Zscore if missing'''
        # Formulas obtained from https://bogdan.dgsom.ucla.edu/news/item?item_id=54624 

        if Col.SE in data.columns:
            #If SE in columns then estimate Z score form beta/se
            data[self.NAME] = data[Col.BETA] / data[Col.SE]
        else:
            #Estimate Z score form sign(beta)*inverse_comulative_normal_distirbution(pval/2)
            data[self.NAME] = np.sign(data[Col.BETA]) * norm.ppf(data[Col.PVAL]/2)                       

        return data


class StandardError(BaseMethod):
    NAME = Col.SE
    REQS = [(Col.Z, Col.BETA)]

    def _method(self, data: DataFrame) -> DataFrame:
        '''Estimate SE from BETA and Zscore if missing'''
        #Estimate SE as beta / Z score
        data[self.NAME] = abs(data[Col.BETA] / data[Col.Z])

        return data

#################################################################################
############################# Method Factory ####################################
#################################################################################
class MethodFactory:
    METHODS = {Col.BETA: EFFECT_SIZE(),
               Col.Z: ZSCORE(),
               Col.SE: StandardError()}

    @classmethod
    def new(cls, type: str):
        return cls.METHODS[type]