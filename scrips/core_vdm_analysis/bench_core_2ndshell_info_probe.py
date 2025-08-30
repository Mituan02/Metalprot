'''
Extract info of 3 contact core 2ndshell from Probe calculation.
The info is in '_probe2ndshell.csv' files.
'''

import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
import pickle
import pandas as pd

database_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_2ndshell_date_reps_probe2ndshell/pickle_noCYS/'

with open(database_dir + 'all_2ndshell_vdms.pkl', 'rb') as f:
    all_2ndshell_vdms = pickle.load(f)

with open(database_dir + 'all_2ndshell_vdm_id_dict.pkl', 'rb') as f:
    all_2ndshell_vdm_id_dict = pickle.load(f)

#Transform the key so that the bench core could be indexed.
list(all_2ndshell_vdm_id_dict.keys())[0]
vdm_id_dict = {}
for key in all_2ndshell_vdm_id_dict.keys():
    newkeysStr = key.split('_centroid_')[1]
    vdm_id_dict[newkeysStr] = all_2ndshell_vdm_id_dict[key]
list(vdm_id_dict.keys())[0]

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_3contact_noCYS/"
files = []
for file in os.listdir(workdir):
    if '.pdb' not in file:
        continue
    files.append(file)


infos = []

database_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_2ndshell_date_reps_probe2ndshell/probe_infos/'
for file in files:
    name = file.split('.')[0]
    df_hbs = pd.read_csv(database_dir + name + '_probe2ndshell.csv', header=0,index_col=0)

    total_bb = 0
    score = 0
    seen = set()
    for i in range(len(df_hbs)):
        chain1 = df_hbs.iloc[i]['chain1']
        resnum1 = int(df_hbs.iloc[i]['resnum1'])
        chain2 = df_hbs.iloc[i]['chain2']
        resnum2 = int(df_hbs.iloc[i]['resnum2'])
        is_bb_or_sc = '_sc_'
        if df_hbs.iloc[i]['name2'] == 'N' or df_hbs.iloc[i]['name2'] == 'O' or df_hbs.iloc[i]['name2'] == 'H':
            is_bb_or_sc = '_bb_'
            total_bb += 1
        tag = is_bb_or_sc +  '_'.join([str(x) for x in df_hbs.iloc[i]])
        key = name + tag

        seen.add((chain1, resnum1))
        try:
            id = vdm_id_dict[key]
            vdm = all_2ndshell_vdms[id]
            score += vdm.score
        except:
            print(key + ' is not in vdm_id_dict.')

        empty = 3 - len(seen)
        total = len(df_hbs)

    infos.append((name, total, total_bb, empty, score))


with open(workdir + '_probeHB_summary.tsv', 'w') as f:
    f.write('Title\tTotal2nd\tTotalbb2nd\tEmpty\tTotal2ndScore\n')

    for info in infos:
        f.write('\t'.join([str(x) for x in info]) + '\n')