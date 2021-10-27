'''

'''


import os
import sys
from numpy.lib.shape_base import split
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.search import extract_vdm
from metalprot.basic import hull
from metalprot.basic import SecondShellVDM
import pickle

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/search_2ndshell/database_2ndshell_pickle.py
'''

vdm_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211018/20211026_2ndshell_selfcenter/'


centroid_vdms = extract_vdm.extract_all_centroid(vdm_dir, summary_name = '_summary.txt', file_name_includes = ['AAMetalPhiPsi', 'cluster'], file_name_not_includes=['CYS'], score_cut = 0, clu_num_cut = 0)


all_2ndshell_vdms = []

max_clu_num = 0
for vdm_id in range(len(centroid_vdms)):
    print(vdm_id)
    _vdm = centroid_vdms[centroid_vdms]
    vdm = SecondShellVDM(_vdm.query.copy(), score = _vdm.score, clu_num = _vdm.clu_num, clu_total_num = _vdm.total_num)
    
    if 'cluster_0' in vdm.query.getTitle():
        max_clu_num = vdm.clu_num

    vdm.id = vdm_id
    vdm.clu_rank = int(vdm.query.getTitle().split('_')[1]) 
    vdm.max_clu_num = max_clu_num
    all_2ndshell_vdms.append(vdm)

outdir = vdm_dir + 'pickle_noCYS/'

with open(outdir + 'all_2ndshell_vdms.pkl', 'wb') as f:
    pickle.dump(all_2ndshell_vdms, f)