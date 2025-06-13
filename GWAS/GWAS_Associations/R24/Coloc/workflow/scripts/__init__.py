import os
import pandas as pd
from copy import deepcopy
from peppy import Project
from subprocess import run
from tempfile import gettempdir

######## Decorators ############
def _internal_decorator(function):
    def wrapper(self, wildcards):
        new_wildcards = deepcopy(wildcards)

        if hasattr(wildcards, "project1") or hasattr(wildcards, "project2"):
            new_wildcards.project = wildcards.project1 if self.name == "Trait1" else wildcards.project2                   #type: ignore

        if hasattr(wildcards, "dataset1") or hasattr(wildcards, "dataset2"):
            new_wildcards.dataset = wildcards.dataset1 if self.name == "Trait1" else wildcards.dataset2                   #type: ignore

        return function(self, new_wildcards)

    return wrapper

##### Data Interfaces #####
### Standard Interface ###
class DATA:
    PROJECT = "project"
    DATASET = "dataset"
    CATEGORY = "Category"

    def __init__(self, pep, attrs: dict, name: str):
        self.name = name
        self.pep = pep
        self.config = attrs

        self._check_type()

    def _make_sample_table(self):
        table = pd.concat([pd.DataFrame(i.to_dict(), index=range(len(i.file))) for i in self.pep.samples])
        table = table.drop_duplicates().reset_index(drop=False)
        self._sample_table = table

    @property
    def sample_table(self):
        if not hasattr(self, "_sample_table"):
            self._make_sample_table()
        return self._sample_table

    @property
    def feature(self):
        return self.config["Feature"]

    def _check_type(self):
        if not self.type in ["quant", "cc"]:
            raise Exception("Type must be either 'quant' or 'cc'")
    
    def _has_categories(self):
        return self.CATEGORY in self.sample_table.columns

    @property
    def type(self):
        return self.config.get("Type", "quant")

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
    
    def Categories(self, unique=False):
        return self.get_column(self.CATEGORY, unique) if self._has_categories() else ["1"]
    
    def get_datasets(self, project: str):
        return self.subset(self.PROJECT, project)["sample_name"].tolist()
    
    @_internal_decorator
    def sample_size(self, wildcards):
        selected = self.subset(self.PROJECT, wildcards.project)

        if "sample_size" in selected.columns:
            n = selected.loc[selected["sample_name"] == wildcards.dataset, "sample_size"].values[0]
            return n if n == n else None
        return None

    @_internal_decorator
    def get_by_category(self, wildcards):
        selected = self.subset(self.PROJECT, wildcards.project)
        if self._has_categories():
            selected = selected[selected[self.CATEGORY] == wildcards.category]
        return selected["sample_name"]
    
    @_internal_decorator
    def get_cc_ratio(self, wildcards):
        if self.type == "cc":
            if not "cc_ratio" in self.sample_table.columns:
                raise Exception("When using 'cc' you need to provide a cc_ratio column in your sample table.")
            
            selected = self.subset(self.PROJECT, wildcards.project)
            value = selected.loc[selected["sample_name"] == wildcards.dataset, "cc_ratio"].values[0]
            return value if value == value else None
        return None

    @_internal_decorator
    def file(self, wildcards):
        selected = self.subset(self.PROJECT, wildcards.project)
        return selected.loc[selected["sample_name"] == wildcards.dataset, "file"].values[0]
    
    @property
    def show_heatmap(self):
        if self.name != "Trait1":
            raise Exception("Heatmap can only be displayed for Trait2...")
        return self.config.get("Heatmap", False)

    @property
    def metadata(self): return self.pep._config.sample_table

### Interface for GWAS data ###        
class GWAS(DATA):
    PROJECT = "disease"

    def __init__(self, pep, config, name):
        super().__init__(pep, config, name)

        if "COLOC" in self.sample_table.columns:
            self.sample_table["COLOC"].replace({"True": True, "False": False}, inplace=True)
            self._sample_table = self.sample_table[self.sample_table["COLOC"]]


################ Type arg makers #################
def _make_abf_args(wds, obj, p: str):
    args = {"d": obj.file(wds),
            "feature": obj.feature,
            "type": obj.type,
            "N": obj.sample_size(wds),
            "cc-ratio": obj.get_cc_ratio(wds)
            }
    return {"--"+key+p: value for key, value in args.items()}


################ Functions #################
def make_run_args(type: str, wildcards, object, prefix: str):
    factories = {"abf": _make_abf_args}

    if not type in factories:
        raise Exception(f"{type} is not a coloc option...\nPlease use one of {list(factories.keys())}")
    
    return factories[type](wildcards, object, prefix)

def clean_logs(path: str):
    if os.path.exists(path) and os.listdir(path):
        run(f"rm -r {path}/*", shell=True)
        
def cp(path: str):
    if not os.path.exists(path):
        raise Exception(f"{path} doesn't exists...")
    cmd = f"cp {path} ./" if os.path.isfile(path) else f"cp -r {path} ./"
    run(cmd, shell=True)

def temp_dir(path: str):
    if path is None:
        user = run("echo `whoami`", capture_output=True, shell=True).stdout
        user = user.decode("UTF-8").strip()
        path = os.path.join(gettempdir(), user)
    
    return os.path.join(path, "Coloc")

def create_interface(name: str, config: dict[str, str]):
    pep = Project(config["pepfile"], amendments=config.get("amend"))                          #type: ignore

    if config.get("Trait-Type") == "GWAS":
        interface = GWAS
    else:
        interface = DATA

    return interface(pep, config, name)

def initialize(config):
    T1 = create_interface("Trait1", config["Trait1"])
    T2 = create_interface("Trait2", config["Trait2"])

    OUT = os.path.join(config['OutDir'], "Coloc")
    return OUT, T1, T2