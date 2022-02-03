'''
Here we want to know how the score of 1st shell and 2nd shell are related.

There is two direction to do this.
For each 1st shell vdm, if there is no 2ndshell, the score is -10. If there are 2nshells, select the one with highest score.
For each 2ndshell, find the 1st shell score.
'''

import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
import pickle
import pandas as pd

### get 1st shell vdm info.
query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/20220128_1stshell/20220128_selfcenter/pickle_noCYS/'

with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)

_1stvdm_id_dict = {}
for q in all_querys:
    ids = q.query.getTitle().split('_centroid_')[1].split('_')
    id = '_'.join([ids[x] for x in [0, 1, 2, 3, 5, 6]])
    _1stvdm_id_dict[id] = q
testkey = list(_1stvdm_id_dict.keys())[0]

### get 2nd shell vdm info.
database_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_2ndshell_date_reps_probe2ndshell/pickle_noCYS/'

with open(database_dir + 'all_2ndshell_vdms.pkl', 'rb') as f:
    all_2ndshell_vdms = pickle.load(f)

with open(database_dir + 'all_2ndshell_vdm_id_dict.pkl', 'rb') as f:
    all_2ndshell_vdm_id_dict = pickle.load(f)

#Transform the key so that the bench core could be indexed.
testkey = list(all_2ndshell_vdm_id_dict.keys())[0]

_2ndvdm_id_dict = {}
for key in all_2ndshell_vdm_id_dict.keys():
    ids = key.split('_centroid_')[1].split('_')
    id = '_'.join([ids[x] for x in [0, 1, 2, 3, 6, 7]])

    if id not in _2ndvdm_id_dict.keys():
        _2ndvdm_id_dict[id] = [all_2ndshell_vdms[all_2ndshell_vdm_id_dict[key]]]
    else:
        _2ndvdm_id_dict[id].append(all_2ndshell_vdms[all_2ndshell_vdm_id_dict[key]])
testkey = list(_2ndvdm_id_dict.keys())[0]


### 1st shell with 2nd shell.
outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/reason/'

with open(outdir + '_2ndshell_1stshell_score_relation.tsv', 'w') as f:
    f.write('_1stvdm\t_1stvdm_aa\t_1stvdm_rank\t_1stvdm_score\t_2ndvdm\t_2ndvdm_aa\t_2ndvdm_rank\t_2ndvdm_score\n')
    for key in _1stvdm_id_dict.keys():
        _1stvdm = _1stvdm_id_dict[key]
        if key in _2ndvdm_id_dict.keys():
            for _2ndvdm in _2ndvdm_id_dict[key]:
                f.write(_1stvdm.query.getTitle() + '\t')
                f.write(_1stvdm.aa_type + '\t')
                f.write(str(_1stvdm.clu_rank) + '\t')
                f.write(str(_1stvdm.score) + '\t')
                f.write(_2ndvdm.query.getTitle() + '\t')
                f.write(_2ndvdm.query.select('resindex 1').getResnames()[0] + '\t')
                f.write(str(_2ndvdm.clu_rank) + '\t')
                f.write(str(_2ndvdm.score) + '\t')
                f.write('\n')
        else:
            f.write(_1stvdm.query.getTitle() + '\t')
            f.write(_1stvdm.aa_type + '\t')
            f.write(str(_1stvdm.clu_rank) + '\t')
            f.write(str(_1stvdm.score) + '\t')
            f.write('\n')


