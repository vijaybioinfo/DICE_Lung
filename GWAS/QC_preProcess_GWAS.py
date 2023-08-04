############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   Quality control for pre-processed GWAS    -----    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

############    -----------------------------------------    ############
### ----------------------- About the script ------------------------ ###
############    -----------------------------------------    ############
# Author: Elizabeth Marquez-Gomez
# Date: 2023-June
# Version: 1
# Subversion: 0
# Using the output GWAS pre-processed files for each disease, this script will: review datatypes and content to be qualified as useful GWAS dataset

#-----> Usage
# python3 QC_preProcess_GWAS.py --disease [name_of_disease(comma separated)]


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############
import pandas as pd
import datatable as dt
import os
import re
import math
import argparse
from datetime import datetime
from R24.Utils import Utilities as Utils
from R24.GWAS.utils.columns import Columns as Cols
import pandas.api.types as types


############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############
## Get diseases to run as argument
parser = argparse.ArgumentParser()
parser.add_argument('--disease', type=str, required=True)
args = parser.parse_args()
diseases = [item for item in args.disease.split(',')]

sourcePathGWAS = '../GWAS'

## Values for CHRs
rangeChr = [str(i) for i in list(range(1,25))]
altChr = set(rangeChr + [True,'True','X','x','Y','y'])


############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############
## Download GWAS datasets by disease
for disease in diseases:
    reviewStudies = 0

    info = dt.fread(f'{sourcePathGWAS}/Metadata_tables/{disease}.csv').to_pandas()
    info.replace({True: 1}, inplace=True)
    info.replace({False: 0}, inplace=True)
    info.fillna(0, inplace=True)
    info['Sum stats'] = info['Sum stats']

    ## Set log file to report disease status
    fileLog = f'{sourcePathGWAS}/Validation/Pre_processed/log_preProcess_GWAS_{disease}.txt'
    open(fileLog, 'a').close()
    now = datetime.now()
    with open(fileLog, 'a') as logFile:
        logFile.write(f'\n\n\t+++++++++++++ ==== STARTING NEW ROUND ===== {now.strftime("%d/%m/%Y %H:%M:%S")} +++++++++++++\n')
        logFile.write(f'\n\n\t-------------- Working on {disease} --------------\n\n')
    print(f'-------------- Working on {disease} --------------')

    ## Set studies column
    info['Studies'] = info['Sample name'].str[:] + "_" + info['Genome version'].astype(str).str[:]

    ## Clean messy datasets
    info.drop(info.index[info['Messy dataset']==1], inplace=True)
    ## Clean non-summary stats
    info.drop(info.index[info['Sum stats']==0], inplace=True)

    for study in info['Studies']:
        checkFlag = False

        subsetStudy = info[info['Studies'] == study]
        ## Get genome version
        gnmVersion = list(subsetStudy['Genome version'])[0]

        with open(fileLog, 'a') as logFile:
            logFile.write(f'------------------------- ==== Working on: {sourcePathGWAS}/{gnmVersion}/pre_process/{disease}/{study}.ssf ===== -------------------------\n')

        print(f'------------------------- ==== Working on: {sourcePathGWAS}/{gnmVersion}/pre_process/{disease}/{study}.ssf ===== -------------------------')
        data = dt.fread(f'{sourcePathGWAS}/{gnmVersion}/pre_process/{disease}/{study}.ssf').to_pandas()

        ## Check CHR column
        if not types.is_integer_dtype(data[Cols.CHR]):
            if data[Cols.CHR].isnull().values.any():
                data.drop(data.index[data[Cols.CHR].isnull().values], inplace=True)
                try:
                    data[Cols.CHR] = data[Cols.CHR].astype(int)
                    checkFlag = True
                    with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING CHR column is int but NAN values present!!! +++++++++++++\n{dataTrashy}\n')
                except: #ValueError
                    checkFlag = True
                    dataTrashy = data[data[Cols.CHR].astype(str).str.match("[^\d]")].head()
                    with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING CHR column not int and NAN values present!!! +++++++++++++\n{dataTrashy}\n')
            else:
                chrs = set(data[Cols.CHR])
                if len(chrs-altChr):
                    checkFlag = True
                    with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING CHR column not int, check values!!! +++++++++++++\n')
        else:
            if (data[Cols.CHR].min() < 1 or data[Cols.CHR].max() > 24):
                checkFlag = True
                with open(fileLog, 'a') as logFile:
                    logFile.write(f'=====\nWARNING CHR column out of range(1,24)!!! +++++++++++++\n{data[Cols.CHR].unique()}\n')

        ## Check POS column
        if not types.is_integer_dtype(data[Cols.POS]):
            if data[Cols.POS].isnull().values.any():
                data.drop(data.index[data[Cols.POS].isnull().values], inplace=True)
                try:
                    data[Cols.POS] = data[Cols.POS].astype(int)
                    checkFlag = True
                    with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING POS column is int but NAN values present!!! +++++++++++++\n{dataTrashy}\n')
                except: #ValueError
                    checkFlag = True
                    dataTrashy = data[data[Cols.POS].astype(str).str.match("[^\d]")].head()
                    with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING POS column not int and NAN values present!!! +++++++++++++\n{dataTrashy}\n')
            else:
                try:
                    data[Cols.POS] = data[Cols.POS].astype(int)
                except: #ValueError
                    checkFlag = True
                    dataTrashy = data[data[Cols.POS].astype(str).str.match("[^\d]")].head()
                    with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING POS column not int!!! +++++++++++++\n{dataTrashy}\n')

        ## Check PVAL column
        if not types.is_float_dtype(data[Cols.PVAL]):
            try:
                data[Cols.PVAL] = data[Cols.PVAL].astype(float)
            except: #ValueError
                checkFlag = True
                with open(fileLog, 'a') as logFile:
                    logFile.write(f'=====\nWARNING PVAL column not float!!! +++++++++++++\n')
        else:
            if (data[Cols.PVAL].min() < 0 or data[Cols.PVAL].max() > 1):
                checkFlag = True
                with open(fileLog, 'a') as logFile:
                    logFile.write(f'=====\nWARNING PVAL column out of range(0,1)!!! +++++++++++++\n')

        ## Check BETA/OR column
        if Cols.BETA in data.columns:
            if not types.is_float_dtype(data[Cols.BETA]):
                checkFlag = True
                with open(fileLog, 'a') as logFile:
                    logFile.write(f'=====\nWARNING BETA column not float!!! +++++++++++++\n')
        elif Cols.OR in data.columns:
            if not types.is_float_dtype(data[Cols.OR]):
                checkFlag = True
                with open(fileLog, 'a') as logFile:
                    logFile.write(f'=====\nWARNING OR column not float!!! +++++++++++++\n')
        else:
            checkFlag = True
            with open(fileLog, 'a') as logFile:
                logFile.write(f'=====\nWARNING missing BETA/OR columns!!! +++++++++++++\n')

        ## Check SE column
        if Cols.SE in data.columns:
            if not types.is_float_dtype(data[Cols.SE]):
                checkFlag = True
                with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING SE column not float!!! +++++++++++++\n')
            else:
                if data[Cols.SE].min() < 0:
                    checkFlag = True
                    with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING SE column <= 0!!! +++++++++++++\n')

        ## Check Z column
        if Cols.Z in data.columns:
            if not types.is_float_dtype(data[Cols.Z]):
                checkFlag = True
                with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING Z column not float!!! +++++++++++++\n')

        ## Check AF column
        if Cols.AF in data.columns:
            if not types.is_float_dtype(data[Cols.AF]):
                checkFlag = True
                with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING AF column not float!!! +++++++++++++\n')
            else:
                if (data[Cols.AF].min() < 0 or data[Cols.AF].max() > 1):
                    checkFlag = True
                    with open(fileLog, 'a') as logFile:
                        logFile.write(f'=====\nWARNING AF column out of range(0,1)!!! +++++++++++++\n')

        if checkFlag:
            reviewStudies += 1

    with open(fileLog, 'a') as logFile:
        logFile.write(f'All done in {disease}! %%%%%%%%%%%%%\n')
        logFile.write(f'Review: {reviewStudies}! :( \n')
    print(f'All done in {disease}! %%%%%%%%%%%%%')
