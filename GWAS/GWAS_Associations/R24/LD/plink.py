import os
import math
import random
import string
import pandas as pd

from subprocess import run
from functools import partial
from typing import Any, Protocol
from collections import defaultdict

shell = partial(run, capture_output=True, shell=True)

################ Helper Functions #################
def split_snps_by_chr(snps: list[str], keep_non_canonical: bool = False) -> dict[str, list[str]]:
    chromosomes = defaultdict(list)
    for snp in snps:
        if not snp: continue
        chrom = snp.split(':')[0]

        canonical = not(chrom.startswith("NW") or chrom.endswith("alt") or chrom.endswith("random"))
        if canonical or keep_non_canonical:                                                                        
            chromosomes[chrom].append(snp)

    return chromosomes

def get_temp_file(prefix: str) -> str:
    ID = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    return f"{prefix}.{ID}"

def check_no_valid_variants_error(error: str):
    message = "Error: No valid variants specified"
    return message in error

def write_mock_file(name: str):
    with open(name, "w") as stream:
        print("CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR2", file=stream)

def write_input_snps(snps: list[str], file: str):
    with open(file, "w") as stream:
        for snp in set(snps):
            print(snp, file=stream)

def get_cmd(file, chromosome, cfg):
    commands = {}
    for pop, popfile in cfg['populations'].items():
        cmd = f"plink --bfile {cfg['bfile']} --r2 "#dprime "  ### modified but needs more work
        cmd = cmd + f"--ld-snp-list {file} "
        cmd = cmd + f"--keep {popfile} " 
        cmd = cmd + f"--ld-window-kb {math.ceil(cfg['window'] / 1000)} "
        cmd = cmd + f"--ld-window {cfg['ld_window']} "
        cmd = cmd + f"--ld-window-r2 {cfg['ld_window_r2']} "
        cmd = cmd + f"--chr {chromosome} --out {file}.{pop}"

        commands[pop] = cmd
    return commands

def create_missing(missing: list[str], pop_cols: list[str]):
    data = {"Lead": list(missing), "LD": list(missing)}
    pops = {pop: [1] * len(missing) for pop in pop_cols}
    data = pd.DataFrame({**data, **pops})
    return data

###################### Merge LD files ############################
def _make_core(files: list[str], threshold: float = 0.8):
    snps = set()

    for file in files:
        with open(file, "r") as stream:
            stream.readline()
            while line := stream.readline().strip().split():
                if float(line[-1]) >= threshold:
                    snps.add((line[2], line[5]))

    return snps

def _merge(files: dict[str, str], threshold: float = 0.8):
    snps = _make_core(files.values(), threshold)
    core = pd.DataFrame(snps, columns=["SNP_A", "SNP_B"])

    for pop, file in files.items():
        tmp = []

        with open(file, "r") as stream:
            stream.readline()
            while line := stream.readline().strip().split():
                if (line[2], line[5]) in snps:
                    tmp.append((line[2], line[5], float(line[-1])))

        #shell(f"rm {file}")
        tmp = pd.DataFrame(tmp, columns=["SNP_A", "SNP_B", pop])
        core = pd.merge(core, tmp, how="left")
    
    return core
    
    
###################### LD Analysis #######################
class PLINK_NOT_FOUND(Exception):
    pass

def check_plink_installation():
    exe = shell("which plink")

    if exe.returncode:
        raise PLINK_NOT_FOUND("plink command line tool not found, please make sure to install plink ...")


def _run_plink(file: str, params: dict[str, Any], chromosome: str):
    ### check plink ###
    check_plink_installation()

    ### estimate ld ###
    ld_files = {}

    for pop, cmd in get_cmd(file, chromosome, params).items():
        exe = shell(cmd)
        ldfile = f"{file}.{pop}.ld"

        if exe.returncode != 0:

            if check_no_valid_variants_error(exe.stderr.decode("UTF-8")):                                                                      #type: ignore
                write_mock_file(ldfile)
            else:
                raise Exception(f"Command failed: {cmd} with returncode {exe.returncode}\n\n" + exe.stderr.decode("UTF-8"))                        #type: ignore
        
        if not os.path.exists(ldfile):
            raise Exception(f"Command finished successfully, but output file does not exist. Command: {cmd} Output: {ldfile}")

        ld_files[pop] = ldfile
        shell(f"rm {file}.{pop}.nosex; \n rm {file}.{pop}.log")

    return ld_files








#####################################################################################################
########################################## LD INTERFACE #############################################
#####################################################################################################
class PLINK:
    TEMPDIR: str = os.path.join(os.getcwd(), "tmp_LD")
    
    #### Configuration ####
    PARAMS = ("bfile", "populations", "window", "ld_window", "ld_window_r2")

    REQUIRED = tuple(["bfile"])
    TYPES = {"bfile": str, "populations": dict, "window": int, "ld_window": int, "ld_window_r2": float}
    DEFAULT = {"populations": None, "window": 250_000, "ld_window": 999_999, "ld_window_r2": 0}

    def __init__(self):
        self.THRESHOLD = 0.8
        self._params = self.DEFAULT

    def initialize(self, params: dict[str, Any] = None):
        self.set_params(params)
        self.ready()

    def _check_requirements(self):
        return [p for p in self.REQUIRED if not p in self._params ]
    
    @property
    def is_ready(self): return not(bool(self._check_requirements()))

    def ready(self):
        if not self.is_ready:
            raise Exception(f"Following arguments must be set before running analysis {self._check_requirements()}. Use set_params method")
        self.set_tempdir()
        
    def set_params(self, params: dict[str, Any] = None):
        new = {p: v for p, v in params.items() if p in self.PARAMS} if not params is None else {}

        #Check Types
        for p, value in new.items():
            if type(value) != self.TYPES[p]:
                raise TypeError(f"Parameter {p} expects type {self.TYPES[p]} but got type {type(value)}")

        self._params.update(new)

    def set_tempdir(self, tempdir: str = None):
        if not tempdir is None:
            self.TEMPDIR = tempdir
        os.makedirs(self.TEMPDIR, exist_ok=True)


    def _method(self, snps: list[str], chromosome: str):
        ##### Write input snps ######
        snp_file = get_temp_file(os.path.join(self.TEMPDIR, chromosome))
        write_input_snps(snps, snp_file)

        ##### Run LD analysis #####
        ld_files = _run_plink(snp_file, self._params, chromosome)
        ld = _merge(ld_files, self.THRESHOLD)
        
        ld = ld.rename(columns={"SNP_A": "Lead", "SNP_B": "LD"})
        ld = ld[["Lead", "LD"] + list(self._params["populations"].keys())]

        shell(f"rm {snp_file}")

        return ld

    def __call__(self, snps: list[str], threshold: float = 0.8, keep_non_canonical: bool = False):
        self.THRESHOLD = threshold
        self.ready()

        if variants := split_snps_by_chr(snps, keep_non_canonical):
            ld = pd.concat([self._method(snps, chromosome) for chromosome, snps in variants.items()])
        else:
            ld = pd.DataFrame(columns=["Lead", "LD"] + list(self._params["populations"].keys()))

        #### Add missing SNPs ########
        missing = list(set(snps) - set(ld["Lead"]))
        
        if missing:
            data = create_missing(missing, self._params["populations"])
            ld = pd.concat([ld, data])

        return ld