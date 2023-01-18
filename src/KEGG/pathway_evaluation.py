'''
=====
Distributed by: Computational Science Initiative, Brookhaven National Laboratory (MIT Liscense)
- Associated publication:
url: 
doi: 
github: 
=====
'''
import copy
import numpy as np
import pandas as pd

from scipy import stats
from pprint import pprint
from scipy.stats import norm
from scipy.stats import ttest_1samp


#####################################
### Step 1: load data
#####################################
data = pd.read_csv('GSE6978_gsnz.csv')
colNameDic = {}
for idx, colName in enumerate(list(data.columns)):
    colNameDic[colName] = idx
dose = data['Dose']

#####################################
### Step 2: load KEGG pathways
#####################################
pathways = pd.read_csv('kegg.pathway.gene.tsv')
pathwayDic = {}
for pathway in pathways.values:
    temp = copy.deepcopy(pathway)
    temp = temp[0]
    temp = temp.split('\t')
    key = temp[0]
    pathwayDic[key] = []
    for val in temp[1:]:
        pathwayDic[key].append(int(val))

#####################################
### Step 3: filter genes
#####################################
removePath = []
for key in pathwayDic:
    geneLst = []
    for gene in pathwayDic[key]:
        if str(gene) in colNameDic:
            geneLst.append(colNameDic[str(gene)])
    if geneLst:
        pathwayDic[key] = geneLst
    else:
        removePath.append(key)

for key in removePath:
    del pathwayDic[key]

sepLineS = '-'*37
sepLineE = '-'*37+'\n'

print(sepLineS)
print('Pathway identification summary')
for item in pathwayDic:
    print('Pathway {} has {} genes'.format(item, len(pathwayDic[item])))
print(sepLineE)

#####################################
### Step 3: pathway activation score
#####################################
idxExp, idxNoneExp = [], []
for idx, item in enumerate(dose):
    if item == '0Gy':
        idxExp.append(idx)
    else:
        idxNoneExp.append(idx)

tscores, pvalues = [], []
data = data.to_numpy()
for idxPathway in pathwayDic:
    pathway = pathwayDic[idxPathway]

    # initialize the log-likelihood ratio for each gene
    activeScore = [0. for _ in range(data.shape[0])]

    for gene in pathway:
        # for the selected gene, take the column out (all samples)
        values = data[:, gene]
        
        # construct conditional distributions dependent on the class label
        dataExp = values[idxExp]
        dataNoneExp = values[idxNoneExp]
        muExp = np.mean(dataExp)
        muNoneExp = np.mean(dataNoneExp)
        stdExp = np.std(dataExp)
        stdNoneExp = np.std(dataNoneExp)

        # compute the log-likelihood ratio
        for idx, value in enumerate(values):
            activeScore[idx] += (np.log(norm(muExp, stdExp).pdf(value)/norm(muNoneExp, stdNoneExp).pdf(value)))
    
    tscore = np.mean(activeScore)/(np.sqrt(np.var(activeScore, ddof=1))*np.sqrt(1/data.shape[0]))
    pvalue = 1 - stats.t.cdf(tscore, df=data.shape[0]-1)

    print('Analysis of pathway {} is completed'.format(idxPathway))

    tscores.append((abs(tscore), idxPathway))
    pvalues.append((pvalue, idxPathway))

tscores = sorted(tscores, key=lambda x:x[0], reverse=True)
pvalues = sorted(pvalues, key=lambda x:x[0])

print(sepLineS)
print('Pathway activation summary')
print('t Statistic')
pprint(tscores)
print('p Value')
pprint(pvalues)
print(sepLineE)

np.savetxt('tscores.txt', tscores)
np.savetxt('pvalues.txt', pvalues)