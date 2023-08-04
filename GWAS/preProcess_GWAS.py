############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   GWAS dataset downloading and parsing for pre-process    -----    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

############    -----------------------------------------    ############
### ----------------------- About the script ------------------------ ###
############    -----------------------------------------    ############
# Author: Elizabeth Marquez-Gomez
# Date: 2022-September
# Version: 1
# Subversion: 0
# Using the metadata tables for each disease, this script will: call "download.sh" script to retrieve the dataset using the url and prepare it for pre-processing. Determining missing columns, status of data and correction of format.

#-----> Usage
# python3 preProcess_GWAS.py --disease [name_of_disease(comma separated)]


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############

import pandas as pd
import pandas.api.types as types
import datatable as dt
import os
import re
import math
import argparse
from datetime import datetime
from R24.Utils import Utilities as Utils
from R24.GWAS.utils.columns import Columns as Cols


############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############

### --------------------------- Download --------------------------- ###
## Dictionary of possible file extensions for GWAS downloads
extensions = {'gz':1, 'txt':2, 'tsv':2, 'linear':2, 'zip':4, 'bgz':3, 'xlsx':4}
def Download(arguments, out):
    '''
        Function to call "download.sh" script to retrieve data
        :param arguments: dictionary; url address and local download path
        :param out: string; destination location for pre-proccessed GWAS
        :return: process execution
    '''
    ## Get file name from download
    if arguments['url'].startswith('/'):
        print('Already downloaded -----------')
        name = arguments['url']
        extFlag = 5 ## already downloaded
    else:
        ## Special downloading
        if '.tar' in arguments['url']:
            extFlag = 4
        ## Determining class of file to download
        else:
            name = arguments['url'].split('/')[-1]
            ext = name.split('.')[-1]
            extFlag = extensions[ext] if ext in extensions.keys() else 5

    ## Declare arguments
    script = f"./download.sh"
    url = f"-url {arguments['url']}"
    location = f"-location {arguments['location']}"
    name = f"-name {name}"
    ext = f"-ext {extFlag}"
    out = f"-out {out}"

    cmd = f"sh {script} {url} {location} {name} {ext} {out}"
    print('Downloading =======================================')
    print(cmd)

    if os.system(cmd) != 0:
        raise Exception(f"Command {cmd} failed!")


### --------------------------- GetGWAS --------------------------- ###
## Dictionary of variations for column names
rawColumns = {Cols.CHR:['Chr','CHROM','CHR','Chromosome','chr','chromosome','#CHR','#CHROM'],
    Cols.POS:['pos','POS','BP','Position','Position_b37','base_pair_location','Position (hg19)','position','base_pair_location_grch37','BP_hg19'],
    Cols.RSID:['rsID','dbSNP rsID','variant_id','SNP','#SNP','rsid','db SNP RS ID/Marker','Variant','RSID','ID','rs_id_all','rs_id'],
    Cols.PVAL:['P','p-value','P-value','Pvalue','P.value_association','P_BOLT_LMM','P_value','p_value','pval_meta','P_GC','p.khor','all_inv_var_meta_p','GC-adjusted P','Logistic-adjusted P-value','frequentist_add_pvalue','p.value','PVALUE','pval','P_Value','p'],
    Cols.BETA:['beta','Beta','BETA','beta_meta','beta_0','all_inv_var_meta_beta','EFFECT','EFFECT1','Effect'],
    Cols.OR:['OR','or.khor','odds_ratio','or_meta','all_OR'],
    Cols.SE:['se','SE','standard_error','se_meta','se.khor','se_0','StdErr','all_inv_var_meta_sebeta','STDERR','sebeta','SE_GC'],
    # Cols.REF:['other_allele','NEA','REF','NonEffect_Allele','A2','ALLELE0','allele_A','ref','Non_coded','REF_allele','Major allele'],
    # Cols.ALT:['effect_allele','ALT','EA','Effect_Allele','A1','ALLELE1','allele_B','alt','Coded','ALT_allele','Minor allele','Allele1'],
    # Cols.REF:['REF','ref','REF_allele'],
    # Cols.ALT:['ALT','alt','ALT_allele'],
    Cols.NEA:['other_allele','NEA','NonEffectAllele','NonEffect_Allele','Non_Effect_allele','non_effect_allele','A2','a2','Allele2','ALLELE0','Non_coded','allele_A','Major allele','Major_allele','major allele','major_allele','ALT','alt','ALT_allele','OTHER_ALLELE'],
    Cols.EA:['effect_allele','EA','EffectAllele','Effect_Allele','Effect_allele','A1','a1','ALLELE1','Coded','Allele1','allele_B','Minor allele','Minor_allele','minor allele','minor_allele','REF','ref','REF_allele','effect_allel','EFFECT_ALLELE'],
    Cols.AF:['EAF','FREQ_Effect_Allele','MAF','A1FREQ','maf1','effect_allele_frequency','Coded_freq','minor_allele_frequency','Freq1', 'AF','FREQ1','af'],
    Cols.Z:['z','Z','Zscore']}

def GetGWAS(file, arguments, out):
    '''
        Function to pre-process GWAS summary statistics datasets
        :param file: string; raw download local file address
        :param arguments: dictionary; destination location pre-process phase, log file to report status by disease, source GWAS database, population ancestry ID
        :param out: string, destination location
        :return: True if succesful process, else False
    '''
    data = dt.fread(file).to_pandas()
    columns = data.columns.tolist()
    ## List to encapsulate columns to extract
    extract = []
    ## Columns to check data type
    checkDtypes = set([Cols.PVAL, Cols.BETA, Cols.OR, Cols.SE, Cols.Z, Cols.AF])

    if data.empty:
        with open(arguments['log'], 'a') as logFile:
            logFile.write(f'==== {file} =====\nERROR:Empty dataset, JUMPING +++++++++++++\n')
        print('Empty dataset, JUMPING +++++++++++++')
        return(False)


    ## Special parsing for UKBB files
    if arguments['source'] == 'Pan-UK Biobank':
        columns = set(data.columns)
        ## The p-values in the summary statistic flat files are natural logged p-values (i.e. expect negative values; see File Information below for paths).
        natLogMeta = set(['beta_meta','pval_meta','se_meta'])
        natLogEur = set(['pval_EUR','beta_EUR','se_EUR'])
        ## The p-values in the summary statistic flat files are negative log10 p-values. P-value columns have been renamed with a "neglog10" suffix to indicate this.
        negLog10Meta = set(['neglog10_pval_meta','beta_meta','se_meta'])
        negLog10Eur = set(['neglog10_pval_EUR','beta_EUR','se_EUR'])

        if natLogMeta&columns == natLogMeta:
            if arguments['pops'] == 'EUR':
                with open(arguments['log'], 'a') as logFile:
                    logFile.write(f'==== {file} =====\nERROR: UKBB meta is only EUR +++++++++++++\n')
                print('UKBB meta is only EUR')
            data = data[['chr','pos','pval_meta','beta_meta','se_meta','ref','alt']]
            data.rename(columns={'chr':Cols.CHR,'pos':Cols.POS,'pval_meta':Cols.PVAL,'beta_meta':Cols.BETA,'se_meta':Cols.SE,'ref':Cols.NEA,'alt':Cols.EA}, inplace=True)
            data[Cols.PVAL] = [math.exp(x) for x in data[Cols.PVAL]]
            extract = list(data.columns)

        elif natLogEur&columns == natLogEur:
            if arguments['pops'] != 'EUR':
                with open(arguments['log'], 'a') as logFile:
                    logFile.write(f'==== {file} =====\nERROR: UKBB EUR is not only EUR +++++++++++++\n')
                print('UKBB EUR is not only EUR')
            data = data[['chr','pos','pval_EUR','beta_EUR','se_EUR','ref','alt']]
            data.rename(columns={'chr':Cols.CHR,'pos':Cols.POS,'pval_EUR':Cols.PVAL,'beta_EUR':Cols.BETA,'se_EUR':Cols.SE,'ref':Cols.NEA,'alt':Cols.EA}, inplace=True)
            data[Cols.PVAL] = [math.exp(x) for x in data[Cols.PVAL]]
            extract = list(data.columns)

        elif negLog10Meta&columns == negLog10Meta:
            if arguments['pops'] == 'EUR':
                with open(arguments['log'], 'a') as logFile:
                    logFile.write(f'==== {file} =====\nERROR: UKBB meta is only EUR +++++++++++++\n')
                print('UKBB meta is only EUR')
            data = data[['chr','pos','neglog10_pval_meta','beta_meta','se_meta','ref','alt']]
            data.rename(columns={'chr':Cols.CHR,'pos':Cols.POS,'neglog10_pval_meta':Cols.PVAL,'beta_meta':Cols.BETA,'se_meta':Cols.SE,'ref':Cols.NEA,'alt':Cols.EA}, inplace=True)
            data[Cols.PVAL] = [10**(-(x)) for x in data[Cols.PVAL]]
            extract = list(data.columns)

        elif negLog10Eur&columns == negLog10Eur:
            if arguments['pops'] != 'EUR':
                with open(arguments['log'], 'a') as logFile:
                    logFile.write(f'==== {file} =====\nERROR: UKBB EUR is not only EUR +++++++++++++\n')
                print('UKBB EUR is not only EUR')
            data = data[['chr','pos','neglog10_pval_EUR','beta_EUR','se_EUR','ref','alt']]
            data.rename(columns={'chr':Cols.CHR,'pos':Cols.POS,'neglog10_pval_EUR':Cols.PVAL,'beta_EUR':Cols.BETA,'se_EUR':Cols.SE,'ref':Cols.NEA,'alt':Cols.EA}, inplace=True)
            data[Cols.PVAL] = [10**(-(x)) for x in data[Cols.PVAL]]
            extract = list(data.columns)


    ## Special parsing for COVID initiative files
    if arguments['source'] == 'COVID-19 HGI':
        data = data[['#CHR','POS','REF','ALT','SNP','all_inv_var_meta_beta','all_inv_var_meta_sebeta','all_inv_var_meta_p','all_meta_AF','rsid']]
        data.rename(columns={'#CHR':Cols.CHR,'POS':Cols.POS,'REF':Cols.NEA,'ALT':Cols.EA,'SNP':Cols.SNPID,'all_inv_var_meta_beta':Cols.BETA,
            'all_inv_var_meta_sebeta':Cols.SE,'all_inv_var_meta_p':Cols.PVAL,'all_meta_AF':Cols.AF,'rsid':Cols.RSID}, inplace=True)
        extract = list(data.columns)


    ## Parsing for NHGRI-EBI Catalog files
    if arguments['source'] == 'NHGRI-EBI Catalog':
        ## Iterate over all columns and retrieve the useful ones, determine whether there are NULL columns
        for key,items in rawColumns.items():
            for item in items:
                if item in columns:
                    if not data[item].isnull().all():
                        data.rename(columns={item:key}, inplace=True)
                        extract.append(key)
                        print(f'Adding column: {extract}')

    ## Special datasets contain chr{chr}:{pos} in column named MarkerName
    if 'MarkerName' in columns:
        if Cols.CHR not in extract and Cols.POS not in extract:
            data[Cols.CHR] = [i.split(':')[0] for i in data.MarkerName]
            data[Cols.CHR] = [i.replace('chr','') for i in data[Cols.CHR]]
            data[Cols.POS] = [i.split(':')[1] for i in data.MarkerName]
            data[Cols.POS] = data[Cols.POS].astype(int)
            extract += [Cols.CHR, Cols.POS]
            with open(arguments['log'], 'a') as logFile:
                logFile.write(f'==== {file} =====\nWARNING: MarkerName used for CHR/POS +++++++++++++\n')
            print('MarkerName used for CHR/POS')

    ## Check for exceptions and incorrect values
    if Cols.BETA not in extract and Cols.OR not in extract:
        with open(arguments['log'], 'a') as logFile:
            logFile.write(f'==== {file} =====\nERROR: No BETA/OR values +++++++++++++\n')
        print('No BETA/OR')

    ## BETA has priority over OR
    if Cols.BETA in extract and Cols.OR in extract:
        extract.remove(Cols.OR)

    if Cols.NEA in extract:
        if types.is_numeric_dtype(data[Cols.NEA]):
            extract.remove(Cols.NEA)
            with open(arguments['log'], 'a') as logFile:
                logFile.write(f'==== {file} =====\nERROR: NEA values are numeric +++++++++++++\n')
            print('NEA values are numeric')
        else:
            data[Cols.NEA] = data[Cols.NEA].str.upper()

    if Cols.EA in extract:
        if types.is_numeric_dtype(data[Cols.EA]):
            extract.remove(Cols.EA)
            with open(arguments['log'], 'a') as logFile:
                logFile.write(f'==== {file} =====\nERROR: EA values are numeric +++++++++++++\n')
            print('EA values are numeric')
        else:
            data[Cols.EA] = data[Cols.EA].str.upper()

    ## Turn to NaN if they don't provide valid alleles
    if Cols.NEA in extract and Cols.EA in extract:
        if ((data[Cols.NEA]=='REF') | (data[Cols.EA]=='REF')).any():
            with open(arguments['log'], 'a') as logFile:
                logFile.write(f'==== {file} =====\nERROR: Wrong allele values +++++++++++++\n')
            print('Wrong alleles values')

    if Cols.NEA in extract and Cols.EA not in extract:
        data[Cols.EA] = float('nan')
        extract.append(Cols.EA)
    if Cols.EA in extract and Cols.NEA not in extract:
        data[Cols.NEA] = float('nan')
        extract.append(Cols.NEA)

    ## Subset useful data
    data = data[extract]

    ## Check datatype from columns [Cols.PVAL, (Cols.BETA | Cols.OR), Cols.SE, Cols.Z, Cols.AF]
    checkDtypes = list(checkDtypes&set(data.columns))
    for element in checkDtypes:
        if not types.is_float_dtype(data[element]):
            if data[element].isnull().values.any():
                dataTrashy = data[data[element].isnull().values]
                data.drop(data.index[data[element].isnull().values], inplace=True)
                data = pd.concat([data, dataTrashy])
                with open(arguments['log'], 'a') as logFile:
                    logFile.write(f'==== {file} =====\nWARNING dType column {element}: NaN values +++++++++++++\n')
                print(f'dType column {element}: NaN values')
            else:
                with open(arguments['log'], 'a') as logFile:
                    logFile.write(f'==== {file} =====\nWARNING dType column {element}: Weird values, CHECK!!! +++++++++++++\n')
                print(f'dType column {element}: Weird values, CHECK!!!')

    ## Check datatypes for POS & CHR, and standarize data
    if Cols.POS in extract:
        if not types.is_integer_dtype(data[Cols.POS]):
            if data[Cols.POS].isnull().values.any():
                dataTrashy = data[data[Cols.POS].isnull().values]
                data.drop(data.index[data[Cols.POS].isnull().values], inplace=True)
                data = pd.concat([data, dataTrashy])
                with open(arguments['log'], 'a') as logFile:
                    logFile.write(f'==== {file} =====\nWARNING dType column {Cols.POS}: NaN values +++++++++++++\n')
                print(f'dType column {Cols.POS}: NaN values')
            else:
                dataTrashy = data[data[Cols.POS].astype(str).str.match('[^\d]')]
                data.drop(data.index[data[Cols.POS].astype(str).str.match('[^\d]')], inplace=True)
                data = pd.concat([data, dataTrashy])
                with open(arguments['log'], 'a') as logFile:
                    logFile.write(f'==== {file} =====\nWARNING dType column {Cols.POS}: Weird values, CHECK!!! +++++++++++++\n')
                print(f'dType column {Cols.POS}: Weird values, CHECK!!!')

    if Cols.CHR in extract:
        chrs = set(data[Cols.CHR])
        if 'True' in chrs:
            data[Cols.CHR] = data[Cols.CHR].str.replace('True', '1')
        if 'X' in chrs or 'x' in chrs:
            data[Cols.CHR] = data[Cols.CHR].str.upper()
            data[Cols.CHR] = data[Cols.CHR].str.replace('X', '23')
        if 'Y' in chrs or 'y' in chrs:
            data[Cols.CHR] = data[Cols.CHR].str.upper()
            data[Cols.CHR] = data[Cols.CHR].str.replace('Y', '24')

        ## Order values
        data[Cols.CHR] = data[Cols.CHR].astype('category')
        order = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        data[Cols.CHR].cat.set_categories(order)
        data.sort_values(by=[Cols.CHR,Cols.POS], inplace=True)

    ## Report missing columns at the end of cleaning
    required = set([Cols.CHR, Cols.POS, Cols.RSID, Cols.PVAL, Cols.BETA, Cols.OR, Cols.SE, Cols.NEA, Cols.EA, Cols.AF, Cols.Z])
    missing = required-set(extract)
    with open(arguments['log'], 'a') as logFile:
        logFile.write(f'==== {file} =====\nMISSING columns: {missing} +++++++++++++\n')
    print(f'Missing columns: {missing}')
    if missing:
        with open(arguments['log'], 'a') as logFile:
            logFile.write(f'COLUMNS: {columns} +++++++++++++\n')
        print(columns)

    ## Write dataset in new location "pre_process" phase
    data.to_csv(f'{arguments["location"]}/{out}', sep='\t', index=False)
    return(True)


### --------------------------- GetNonSummStats --------------------------- ###
def GetNonSummStats(file, arguments, out):
    '''
        Function to pre-process GWAS non-summary statistics datasets
        :param file: string; raw download local file address
        :param arguments: dictionary; destination location pre-process phase, log file to report status by disease
        :param out: string, destination location
        :return: True if succesful process, else False
    '''
    with open(arguments['log'], 'a') as logFile:
        logFile.write(f'Non-summary: {file} --> {arguments["location"]}/{out} +++++++++++++\n')
    print(f'Non-summary: {file} --> {arguments["location"]}/{out}')

    os.system(f'cp {file} {arguments["location"]}/{out}')
    return(True)


############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############
parser = argparse.ArgumentParser()
parser.add_argument('--disease', type=str, required=True)
args = parser.parse_args()
diseases = [item for item in args.disease.split(',')]

rootPathGWAS = '../GWAS'

############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############
## Download GWAS datasets by disease
for disease in diseases:
    print(f'-------------- Working on {disease} --------------\n')
    info = dt.fread(f'{rootPathGWAS}/Metadata_tables/{disease}.csv').to_pandas()
    info.replace({True: 1}, inplace=True)
    info.replace({False: 0}, inplace=True)
    info.fillna(0, inplace=True)
    info['Sum stats'] = info['Sum stats'].astype(str)
    info[['n_cases','n_controls','sample_size','AFR_Pop','AMR_Pop','EUR_Pop','EAS_Pop','SAS_Pop','Category']] = info[['n_cases','n_controls','sample_size','AFR_Pop','AMR_Pop','EUR_Pop','EAS_Pop','SAS_Pop','Category']].astype(int)

    ## Set genome versions to work on
    versions = [version.replace('hg','') for version in list(set(info['Genome version']))]
    ## Drop empty rows
    versions = list(filter(None, versions))

    ## Create structure by genome version
    for version in versions:
        print(f'Creating structure hg{version}-------------------------')
        Utils.mkdir(f'{rootPathGWAS}/hg{version}/raw/{disease}')
        Utils.mkdir(f'{rootPathGWAS}/hg{version}/pre_process/{disease}')

    ## Set log file to report disease status
    logFile = f'{rootPathGWAS}/Logs/log_GWAS_{disease}.txt'
    open(logFile, 'a').close()
    now = datetime.now()
    with open(logFile, 'a') as fileLog:
        fileLog.write(f'\n\n\t+++++++++++++ ==== STARTING NEW ROUND ===== {now.strftime("%d/%m/%Y %H:%M:%S")} +++++++++++++\n')

    ## Set file names column
    info['Studies'] = info['Sample name'].str[:] + "_" + info['Genome version'].astype(str).str[:]

    for study in info['Studies']:
        print('Subset info -------------------------')
        subsetStudy = info[info['Studies'] == study]
        print(subsetStudy)

        ## Check if data available
        subsetStudy['Messy dataset'] = subsetStudy['Messy dataset'].astype(str)
        if (list(subsetStudy['Messy dataset'])[0] == '1') or (list(subsetStudy['Messy dataset'])[0] == 'nc'):
            print('Messy dataset -- JUMPING dataset ---------------')
        else:
            ## Download original file and parse, determine genome version and type of summary statistics dataset
            gnmVersion = list(subsetStudy['Genome version'])[0].replace('hg','')
            extension = 'ssf' if list(subsetStudy['Sum stats'])[0] == '1' else 'nss'
            nameOut = f'{study}.tsv'

            location = f'{rootPathGWAS}/hg{gnmVersion}/raw/{disease}'
            print(f'{location}/{nameOut}')
            if not os.path.exists(f'{location}/{nameOut}'):
                url = list(subsetStudy['Link'])[0]
                arguments = {'url':url, 'location':location}
                print(f'Downloading: {arguments} --> {nameOut}')
                Download(arguments, nameOut)

            ## File variable will save path to downloaded file
            file = f'{location}/{nameOut}'
            nameOut = nameOut.replace('tsv', extension)

            location = f'{rootPathGWAS}/hg{gnmVersion}/pre_process/{disease}'

            if not os.path.exists(f'{location}/{nameOut}'):
                if os.path.exists(file):
                    ## Parse GWAS file and determine protocol for parsing
                    if extension == 'ssf':
                        source = list(subsetStudy['Source'])[0]
                        pops = list(subsetStudy['populations'])[0]
                        arguments = {'location':location, 'log':logFile, 'source':source, 'pops':pops}
                        print(f'Parsing GWAS: {file} == {arguments} --> {nameOut}')
                        if not GetGWAS(file, arguments, nameOut):
                            print(f'JUMPING {file}')
                    else:
                        arguments = {'location':location, 'log':logFile}
                        print(f'Transferring GWAS: {file} == {arguments} --> {nameOut}')
                        if not GetNonSummStats(file, arguments, nameOut):
                            print(f'JUMPING {file}')
                else:
                    with open(logFile, 'a') as fileLog:
                        fileLog.write(f'==== {file} =====\nERROR: Download does not exist\n{url} +++++++++++++\n')

    with open(logFile, 'a') as fileLog:
        fileLog.write(f'All done in {disease}! %%%%%%%%%%%%%')


print('All done! %%%%%%%%%%%%%')
