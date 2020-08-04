##########################################################
## ccAF:  test_ccAF.py                                  ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor          ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

## 0. Load up packages
# General
from os.path import exists
import numpy as np
import pandas as pd
#import pickle

# Single Cell Packages
import scanpy as sc

# Import ccAF for classification of cell cycle phase
import ccAF

## 1. Functions
# Make a confusion matrix
def confusionMatrix(a1, b1):
    v1 = set(a1)
    v2 = set(b1)
    dfTmp = pd.DataFrame([pd.Series(list(a1)),pd.Series(list(b1))]).T
    df1 = pd.DataFrame(index=v1, columns = v2)
    for i in v1:
        for j in v2:
            df1.loc[i,j] = dfTmp.loc[(dfTmp[0]==i) & (dfTmp[1]==j)].shape[0]
    return df1

## 2. Load up conversion dictionaries for gene symbol to Ensembl gene ID
# HGNC -> downlaoded from HGNC website (https://www.genenames.org/download/custom/)
hgncEnsmbl = pd.read_csv('data/hgnc_geneSymbols_ensmbl.txt', index_col=1, header=0, sep='\t')
hgncEnsmbl = hgncEnsmbl.loc[~hgncEnsmbl['Ensembl ID(supplied by Ensembl)'].isnull()]

ensmblHgnc = pd.Series(hgncEnsmbl.index)
ensmblHgnc.index = list(hgncEnsmbl['Ensembl ID(supplied by Ensembl)'])

hgncPrevEnsmbl = {}
for i in hgncEnsmbl.loc[~hgncEnsmbl['Previous symbols'].isnull()].index:
    splitUp = hgncEnsmbl.loc[i,'Previous symbols'].split(', ')
    ensmbl = hgncEnsmbl.loc[i,'Ensembl ID(supplied by Ensembl)']
    for j in splitUp:
        hgncPrevEnsmbl[j] = ensmbl

hgncAliasEnsmbl = {}
for i in hgncEnsmbl.loc[~hgncEnsmbl['Alias symbols'].isnull()].index:
    splitUp = hgncEnsmbl.loc[i,'Alias symbols'].split(', ')
    ensmbl = hgncEnsmbl.loc[i,'Ensembl ID(supplied by Ensembl)']
    for j in splitUp:
        hgncAliasEnsmbl[j] = ensmbl


## 3. Load up data
# Load up loom file
print('\nLoading WT data...')
adata1 = sc.read_loom('data/WT.loom')


## 4. Convert to Ensembl gene IDs
# Convert from gene symbols to Ensmbl gene IDs
adata1 = adata1[:,[True if i in hgncEnsmbl.index or i in hgncPrevEnsmbl or i in hgncAliasEnsmbl else False for i in adata1.var_names]]
convertMe_hgnc = pd.DataFrame(np.nan, index=adata1.var_names, columns=['ensembl'])
numConverted_hgnc = 0
missed = []
for j in adata1.var_names:
    if j in hgncEnsmbl.index:
        convertMe_hgnc.loc[j,'ensembl'] = hgncEnsmbl.loc[j,'Ensembl ID(supplied by Ensembl)']
    elif j in hgncPrevEnsmbl:
        convertMe_hgnc.loc[j,'ensembl'] = hgncPrevEnsmbl[j]
    elif j in hgncAliasEnsmbl:
        convertMe_hgnc.loc[j,'ensembl'] = hgncAliasEnsmbl[j]
    else:
        missed.append(j)

adata1.var_names = pd.Index(convertMe_hgnc['ensembl'])


## 5. ACTINN Classification
# Run ccAF ACTINN classifier against data
print('\nACTINN...')
predictedLabels = ccAF.ccAF.predict_labels(adata1)
adata1.obs['ccAF'] = predictedLabels
adata1.write(filename='WT_ccAF.h5ad')


## 6. Downstream analyses
# Dataframe of true labels, predictions, probabilities for all iterations
DF = pd.concat([adata1.obs, pd.DataFrame({'Predictions':predictedLabels},index=adata1.obs_names) ], axis=1)
DF.to_csv('ccAF_results.csv')

# Make a confusion matrix table
t1 = pd.DataFrame(pd.Series(predictedLabels).value_counts()).T
t1.index = ['Counts']
t2 = confusionMatrix(adata1.obs['clusts_named'], predictedLabels)
table1 = pd.concat([t1, t2[list(t1.columns)]])
table1.to_csv('confusionMatrix_table.csv')

