'''
The ntf2 protein 1ogz has many strcutures. This script trying to summarize the results. 
'''
import os
import pandas as pd
import re

workdir = '/mnt/e/DesignData/ligands/LigandBB/ntf2/family_3vsy_align/'

all_pdb_win_clu = []

for folder in os.listdir(workdir):
    if '.' in folder:
        continue
    info = folder.split('_')
    title = info[2]
    ps = info[3] + info[4]
    file = workdir + folder + '/' + '_summary.tsv'

    win_clu = {}
    with open(file, 'r') as f:
        for line in f.readlines():
            if 'Wins' in line:
                continue
            kk = line.split('\t')
            key = kk[0]
            _clu = re.split('_|-', kk[1])
            clu = _clu[0] + '-' + _clu[2] + '-' + _clu[4]
            if key not in win_clu.keys():
                win_clu[key]= []
            win_clu[key].append(clu)
    all_pdb_win_clu.append((title, ps, win_clu))


# print all keys
keys = set()

for win_clu in all_pdb_win_clu:
    for k in win_clu[2].keys():
        keys.add(k)
print(keys)
#{'13_17_25', '18_26_55', '17_25_54', '15_19_27', '19_27_56', '14_18_26', '13_17_54'}
A = {'13_17_25', '15_19_27', '14_18_26'}
B = {'13_17_54'}
C = {'17_25_54', '18_26_55', '19_27_56'}
with open(workdir + '_summary.tsv', 'w') as f:
    f.write('Name\tPS\twin_13-17-25\twin_13-17-54\twin_17-25-54\n')
    for win_clu in all_pdb_win_clu:
        f.write(win_clu[0] + '\t' + win_clu[1] + '\t')
        for k in win_clu[2].keys():
            if k in A:
                a = '||'.join(win_clu[2][k]) + '\t'
            else:
                a = 'NULL\t'
            if k in B:
                b = '||'.join(win_clu[2][k]) + '\t'
            else:
                b = 'NULL\t'
            if k in C:
                c = '||'.join(win_clu[2][k]) + '\n'
            else:
                c = 'NULL\n'
        f.write( a + b + c)
        