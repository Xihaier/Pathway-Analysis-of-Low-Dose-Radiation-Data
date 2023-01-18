'''
=====
Distributed by: Computational Science Initiative, Brookhaven National Laboratory (MIT Liscense)
- Associated publication:
url: 
doi: 
github: 
=====
'''
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.family': 'serif', 'font.size': 37})    


# -------------------------------------
# helper function
def get_num(score, idces):
    res = []
    for idx in idces:
        res.append(score[idx])
    mu = np.mean(res)
    sigma = np.std(res)
    return [num for num in res if num>mu-3*sigma and num<mu+3*sigma]


GSE6978 = False
if GSE6978:
    # -------------------------------------
    # GSE6978
    dic_LLR = np.load('results/GSE6978/dic_LLR.npy', allow_pickle=True)

    dic_idx_0 = np.load('results/GSE6978/dic_idx_0.npy', allow_pickle=True)
    dic_idx_low = np.load('results/GSE6978/dic_idx_low.npy', allow_pickle=True)
    dic_idx_high = np.load('results/GSE6978/dic_idx_high.npy', allow_pickle=True)

    idx_0 = dic_idx_0.item().get('6978')
    idx_low = dic_idx_low.item().get('6978')
    idx_high = dic_idx_high.item().get('6978')

    tar_low_pathways = ['hsa04650', 'hsa04115', 'hsa03030', 'hsa03420', 'hsa01524', 'hsa00471', 'hsa04210', 'hsa03430', 'hsa03410', 'hsa04215', 'hsa05130', 'hsa05132', 'hsa05142', 'hsa04979', 'hsa05417', 'hsa05200', 'hsa05135', 'hsa04962', 'hsa04928']

    for pathway in tar_low_pathways:
        dat = dic_LLR.item().get('6978_{}'.format(pathway))
        fig, ax = plt.subplots(figsize=(7, 12))
        plt.scatter(np.ones_like(get_num(dat, idx_0))*0, get_num(dat, idx_0), s=270, alpha=0.3, label='0')
        plt.scatter(np.ones_like(get_num(dat, idx_low))*1, get_num(dat, idx_low), s=270, alpha=0.3, label='0.05')
        plt.scatter(np.ones_like(get_num(dat, idx_high))*2, get_num(dat, idx_high), s=270, alpha=0.3, label='0.5')
        plt.xticks(np.arange(3), ['0', '0.05', '0.5'])
        plt.xlabel('Dose level')
        plt.savefig('{}.png'.format(pathway), bbox_inches='tight')
else:
    # -------------------------------------
    # GSE43151
    dic_LLR = np.load('results/GSE43151/dic_LLR.npy', allow_pickle=True)

    dic_idx_0 = np.load('results/GSE43151/dic_idx_0.npy', allow_pickle=True)
    dic_idx_low = np.load('results/GSE43151/dic_idx_low.npy', allow_pickle=True)
    dic_idx_high = np.load('results/GSE43151/dic_idx_high.npy', allow_pickle=True)

    idx_0 = dic_idx_0.item().get('43151_005')
    idx_low_005 = dic_idx_low.item().get('43151_005')
    idx_low_01 = dic_idx_low.item().get('43151_01')
    idx_low_025 = dic_idx_low.item().get('43151_025')
    idx_low_05 = dic_idx_low.item().get('43151_05')
    idx_low_1 = dic_idx_low.item().get('43151_1')
    idx_high = dic_idx_high.item().get('43151_005')

    tar_low_pathways = ['hsa04650', 'hsa04520', 'hsa00600', 'hsa05165', 'hsa00601', 'hsa05130', 'hsa04612', 'hsa05332', 'hsa04630', 'hsa04215', 'hsa05130', 'hsa05132', 'hsa05142', 'hsa04979', 'hsa05202', 'hsa05203', 'hsa04145', 'hsa05320', 'hsa04218', 'hsa04660', 'hsa04940', 'hsa05332', 'hsa05330', 'hsa05416']

    for pathway in tar_low_pathways:
        dat = dic_LLR.item().get('43151_all_{}'.format(pathway))
        fig, ax = plt.subplots(figsize=(21, 11.6))
        plt.scatter(np.ones_like(get_num(dat, idx_0))*0, get_num(dat, idx_0), s=270, alpha=0.3, label='0.005')
        plt.scatter(np.ones_like(get_num(dat, idx_low_005))*1, get_num(dat, idx_low_005), s=270, alpha=0.3, label='0.01')
        plt.scatter(np.ones_like(get_num(dat, idx_low_01))*2, get_num(dat, idx_low_01), s=270, alpha=0.3, label='0.025')
        plt.scatter(np.ones_like(get_num(dat, idx_low_025))*3, get_num(dat, idx_low_025), s=270, alpha=0.3, label='0.05')
        plt.scatter(np.ones_like(get_num(dat, idx_low_05))*4, get_num(dat, idx_low_05), s=270, alpha=0.3, label='0.1')
        plt.scatter(np.ones_like(get_num(dat, idx_low_1))*5, get_num(dat, idx_low_1), s=270, alpha=0.3, label='0.1')
        plt.scatter(np.ones_like(get_num(dat, idx_high))*6, get_num(dat, idx_high), s=270, alpha=0.3, label='0.1')
        plt.xticks(np.arange(7), ['0', '0.005', '0.01', '0.025', '0.05', '0.1', '0.5'])
        plt.xlabel('Dose level')
        plt.savefig('{}.png'.format(pathway), bbox_inches='tight')
        plt.close()