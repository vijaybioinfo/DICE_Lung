import os
import pandas as pd

from peppy import Project
from subprocess import run
from enum import Enum, auto
from tempfile import gettempdir
from collections import namedtuple

def _internal_decorator(function):
    def wrapper(self, wildcards):
        new_wildcards = namedtuple("wildcards", ["project", "dataset"])

        if hasattr(wildcards, "disease") and hasattr(wildcards, "dataset"):
            new_wildcards.project = wildcards.disease
            new_wildcards.dataset = wildcards.dataset

        if hasattr(wildcards, "project") and hasattr(wildcards, "sample"):
            new_wildcards.project = wildcards.project
            new_wildcards.dataset = wildcards.sample

        return function(self, new_wildcards)

    return wrapper

#######################################################################
############################ Interfaces ###############################
#######################################################################
class AvailableAnalyses(Enum):
    SummaryStats = auto()
    NonSummaryStats = auto()

########################
######### BASE #########
########################
class DATA:
    PROJECT = "project"
    DATASET = "dataset"

    def __init__(self, pep: Project) -> None:
        self.pep = pep

    def _make_sample_table(self):
        table = pd.concat([pd.DataFrame(i.to_dict(), index=range(len(i.file))) for i in self.pep.samples])
        table = table.drop_duplicates().reset_index(drop=False)
        self._sample_table = table

    @property
    def sample_table(self):
        if not hasattr(self, "_sample_table"):
            self._make_sample_table()
        return self._sample_table

    def subset(self, column: str, value: str):
        return self.sample_table[self.sample_table[column] == value]

    def get_column(self, column: str, unique=False):
        if unique:
            return self.sample_table[column].drop_duplicates()
        return self.sample_table[column]

    def Projects(self, unique=False):
        return self.get_column(self.PROJECT, unique)

    def Datasets(self, unique: bool = False):
        return self.get_column("sample_name", unique)
    
    def get_datasets(self, project: str):
        return self.subset(self.PROJECT, project)["sample_name"].tolist()

    def file(self, project: str, dataset: str):
        selected = self.subset(self.PROJECT, project)
        return selected.loc[selected["sample_name"] == dataset, "file"].values[0]

    @property
    def get_metadata(self) -> str:
        return self.pep._config.sample_table

    @_internal_decorator
    def N(self, wildcards):
        selected = self.subset(self.PROJECT, wildcards.project)
        return selected.loc[selected["sample_name"] == wildcards.dataset, "sample_size"].values[0]

########################
######### GWAS #########
########################
class GWAS(DATA):
    PROJECT = "disease"

    def __init__(self, pep):
        super().__init__(pep)

        if "Sum stats" in self.sample_table.columns:
            self.sample_table["MESC"].replace({"True": True, "False": False}, inplace=True)
            self._sample_table = self.sample_table[self.sample_table["MESC"]]

            self.sample_table["sample_size"] = self.sample_table["sample_size"].astype(float)
            self._sample_table = self.sample_table[~self.sample_table["sample_size"].isna()]
            self._sample_table = self.sample_table[self.sample_table["sample_size"] != None]
            self._sample_table = self.sample_table[self.sample_table["sample_size"] != 0]

    def file(self, wildcards):
        return super().file(wildcards.disease, wildcards.dataset)

    def Diseases(self, unique: bool=False):
        return self.Projects(unique)

    def get_datasets(self, wildcards) -> list[str]:
        return super().get_datasets(wildcards.disease)

    def Categories(self, unique=False): return self.get_column("Category", unique)
    
    def get_by_category(self, wildcards):
        selected = self.subset(self.PROJECT, wildcards.disease)
        return selected[selected["Category"] == wildcards.category]["sample_name"]


###########################
######### Project #########
###########################
class InternalProject(DATA):
    def get_project_datasets(self, project: str) -> list[str]:
        return self.get_datasets(project)

    def file(self, wildcards):
        return super().file(wildcards.project, wildcards.sample)

class SummaryStats(InternalProject):
    pass

class NonSummaryStats(InternalProject):
    GROUP = "group"

    def Groups(self, unique=False):
        if unique:
            return self.sample_table["group"].unique().tolist() + ["All"]
        return self.sample_table["group"].tolist() + ["All"]

    def get_group_samples(self, wildcards):
        if wildcards.group == "All":
            return list(set(self.get_project_datasets(wildcards.project)))
        selected = self.subset(self.GROUP, wildcards.group)
        return selected["sample_name"].unique().tolist()


################################################
############### CONFIGURATION ##################
################################################
class CONFIGURATION:
    def __init__(self, config: dict):
        self.config = config

    @property
    def prefix(self):
        return self.config["OutDir"]

    @property
    def software(self): return self.config["SOFTWARE"]

    @property
    def analysis_type(self): return self.config["ANALYSIS_TYPE"]


class SummaryStatsConfiguration(CONFIGURATION):
    @property
    def annotation(self): return self.config["ANNOTATION"]

    @property
    def DB(self): return self.config["SNP_DATABASE"]


class NonSummaryStatsConfiguration(CONFIGURATION):
    @property
    def genotype(self): return self.config["GENOTYPE"]


#################################################################################################################################################

################ Functions #################
def clean_logs(path: str):
    if os.path.exists(path) and os.listdir(path):
        run(f"rm -r {path}/*", shell=True)

def temp_dir(path: str):
    if path is None:
        user = run("echo `whoami`", capture_output=True, shell=True).stdout
        user = user.decode("UTF-8").strip()
        path = os.path.join(gettempdir(), user)
    
    return os.path.join(path, "MESC")

def make_pep_project(config: dict, key: str):
    return Project(config[key]["pepfile"], amendments=config[key].get("amend"))



################ Interface Builders #################
def build_gwas_interface(config):
    return GWAS(make_pep_project(config, "GWAS_PEP"))

def build_configuration_interface(config):
    if not config["ANALYSIS_TYPE"] in AvailableAnalyses.__members__:
        raise ValueError(f"Analysis type: {config['ANALYSIS_TYPE']}, choose one from {AVAILABLE_ANALYSES}")

    if config["ANALYSIS_TYPE"] == AvailableAnalyses.SummaryStats.name:
        return SummaryStatsConfiguration(config)
    
    return NonSummaryStatsConfiguration(config)

def build_project_interface(config):
    if not config["ANALYSIS_TYPE"] in AvailableAnalyses.__members__:
        raise ValueError(f"Analysis type: {config['ANALYSIS_TYPE']}, choose one from {AVAILABLE_ANALYSES}")

    if config["ANALYSIS_TYPE"] == AvailableAnalyses.SummaryStats.name:
        return SummaryStats(make_pep_project(config, "PROJECT_PEP"))
    
    return NonSummaryStats(make_pep_project(config, "PROJECT_PEP"))




########## Analysis type helpers ##########
def test_summary_stats(CFG):
    ### Test annotation ###
    if not CFG.annotation: raise Exception("Annotation is none")

    ### Test annotation ###
    if not CFG.DB: raise Exception("SNP database is none")

def test_non_summary_stats(CFG, PROJECT):
    ### Test Groups ###
    if not PROJECT.GROUP in PROJECT.sample_table.columns:
        raise Exception("group column not present in sample table")

    ### Test Genotype ###
    if not CFG.genotype["COHORT"]: raise Exception("Genotype info is missing")
    if not CFG.genotype["PANEL"]: raise Exception("Reference Panel Genotype info is missing")

    ### Test software ###
    if not CFG.software["PLINK"]: raise Exception("Path to plink software is missing")
    if not CFG.software["META_ANALYSIS"]: raise Exception("Path to MESC meta_analyze_weights software is missing")




######### Initialize ###########
def initialize(config: dict):
    CFG = build_configuration_interface(config)
    GS = build_gwas_interface(config)
    PROJECT = build_project_interface(config)

     ### Check Interface Type ###
    if CFG.analysis_type == AvailableAnalyses.SummaryStats.name:
        test_summary_stats(CFG)
    else:
        test_non_summary_stats(CFG, PROJECT)

    return CFG, GS, PROJECT