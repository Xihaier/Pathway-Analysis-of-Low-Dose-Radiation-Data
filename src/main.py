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
import logging
import numpy as np
import pandas as pd

from scipy import stats
from pprint import pprint
from scipy.stats import norm


#####################################
### Step 0: set up the logging
#####################################
tasks = '43151_all'
logging.basicConfig(filename='loggings.log', level=logging.INFO)
sepLineS = '-'*37
sepLineE = '-'*37+'\n'


for task in tasks:
    #####################################
    ### Step 1: load data
    #####################################
    if '6978' in task:
        data = pd.read_csv('GSE/GSE6978_gsnz.csv')
    elif '43151' in task:
        data = pd.read_csv('GSE/GSE43151_gs.csv')
    else:
        raise TypeError
    colNameDic = {}
    for idx, colName in enumerate(list(data.columns)):
        colNameDic[colName] = idx
    dose = data['Dose']

    logging.info(sepLineS)
    logging.info('Task {} begins'.format(task))
    logging.info(sepLineE)


    #####################################
    ### Step 2: load KEGG pathways
    #####################################
    pathways = pd.read_csv('GSE/kegg.pathway.gene.tsv')
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

    logging.info(sepLineS)
    logging.info('Pathway identification summary')
    for item in pathwayDic:
        logging.info('Pathway {} has {} genes'.format(item, len(pathwayDic[item])))
    logging.info(sepLineE)


    #####################################
    ### Step 4: pathway activation score
    #####################################
    idxExp, idxLowExp, idxHighExp = [], [], []

    if task == '43151_all':
        for idx, item in enumerate(dose):
            if item == '0.Gy':
                idxExp.append(idx)
            elif item == '0.005Gy':
                idxLowExp.append(idx)
            elif item == '0.01Gy':
                idxLowExp.append(idx)
            elif item == '0.025Gy':
                idxLowExp.append(idx)
            elif item == '0.05Gy':
                idxLowExp.append(idx)
            elif item == '0.1Gy':
                idxHighExp.append(idx)
            elif item == '0.5Gy':
                idxHighExp.append(idx)
            else:
                raise TypeError
    else:
        raise TypeError

    tscores_low, pvalues_low = [], []
    tscores_high, pvalues_high = [], []

    data = data.to_numpy()
    for idxPathway in pathwayDic:
        pathway = pathwayDic[idxPathway]

        # initialize the log-likelihood ratio for each gene
        activeScoreLow = [0. for _ in range(data.shape[0])]
        activeScoreHigh = [0. for _ in range(data.shape[0])]

        for gene in pathway:
            # for the selected gene, take the column out (all samples)
            values = data[:, gene]
            
            # construct conditional distributions dependent on the class label
            dataExp = values[idxExp]
            dataLowExp = values[idxLowExp]
            dataHighExp = values[idxHighExp]

            muExp = np.mean(dataExp)
            muLowExp = np.mean(dataLowExp)
            muHighExp = np.mean(dataHighExp)

            stdExp = np.std(dataExp)
            stdLowExp = np.std(dataLowExp)
            stdHighExp = np.std(dataHighExp)

            # compute the log-likelihood ratio
            for idx, value in enumerate(values):
                activeScoreLow[idx] += (np.log(norm(muExp, stdExp).pdf(value)/norm(muLowExp, stdLowExp).pdf(value)))
                activeScoreHigh[idx] += (np.log(norm(muExp, stdExp).pdf(value)/norm(muHighExp, stdHighExp).pdf(value)))
        
        tscore = np.mean(activeScoreLow)/(np.sqrt(np.var(activeScoreLow, ddof=1))*np.sqrt(1/data.shape[0]))
        pvalue = 1 - stats.t.cdf(tscore, df=data.shape[0]-1)
        tscores_low.append((abs(tscore), idxPathway))
        pvalues_low.append((pvalue, idxPathway))

        tscore = np.mean(activeScoreHigh)/(np.sqrt(np.var(activeScoreHigh, ddof=1))*np.sqrt(1/data.shape[0]))
        pvalue = 1 - stats.t.cdf(tscore, df=data.shape[0]-1)
        tscores_high.append((abs(tscore), idxPathway))
        pvalues_high.append((pvalue, idxPathway))

        logging.info('Analysis of pathway {} is completed'.format(idxPathway))
        print('Analysis of pathway {} is completed'.format(idxPathway))

    tscores = sorted(tscores_low, key=lambda x:x[0], reverse=True)
    pvalues = sorted(pvalues_low, key=lambda x:x[0])

    logging.info(sepLineS)
    logging.info('Low-dose pathway activation results summary')
    logging.info('t Statistic')
    logging.info(tscores)
    logging.info('p Value')
    logging.info(pvalues)
    logging.info(sepLineE)

    print(sepLineS)
    print('Low-dose pathway activation results summary')
    print('t Statistic')
    pprint(tscores)
    print('p Value')
    pprint(pvalues)
    print(sepLineE)

    tscores = sorted(tscores_high, key=lambda x:x[0], reverse=True)
    pvalues = sorted(pvalues_high, key=lambda x:x[0])

    logging.info(sepLineS)
    logging.info('High-dose pathway activation results summary')
    logging.info('t Statistic')
    logging.info(tscores)
    logging.info('p Value')
    logging.info(pvalues)
    logging.info(sepLineE)

    print(sepLineS)
    print('High-dose pathway activation results summary')
    print('t Statistic')
    pprint(tscores)
    print('p Value')
    pprint(pvalues)
    print(sepLineE)