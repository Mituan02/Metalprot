'''
Transform the 2ndshell cluster database into pickle file.
'''


import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.search import extract_vdm
from metalprot.basic import vdmer_2ndshell
import pickle

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/search_2ndshell/database_2ndshell_pickle.py
'''

vdm_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/20220116_selfcenter_bb2ndshell_notconnect/'
outdir = vdm_dir + 'pickle_noCYS/'
os.makedirs(outdir, exist_ok=True)

centroid_vdms = extract_vdm.extract_all_centroid(vdm_dir, summary_name = '_summary.txt', file_name_includes = ['AA2ndS', 'cluster'], file_name_not_includes=['CYS'], score_cut = 0, clu_num_cut = 0)

### Generate all_2ndshell_vdms.
root_vdm = centroid_vdms[0]

all_2ndshell_vdms = []
all_2ndshell_vdm_id_dict = {}

max_clu_num = 0
for vdm_id in range(len(centroid_vdms)):
    _vdm = centroid_vdms[vdm_id]
    #print(vdm_id)
    try:
        transform = pr.calcTransformation(_vdm.query.select('resindex 1 and name N C CA'), root_vdm.query.select('resindex 1 and name N C CA'))
        transform.apply(_vdm.query)
    except:
        print(_vdm.query.getTitle())
        continue

    vdm = vdmer_2ndshell.SecondShellVDM(_vdm.query, score = _vdm.score, clu_num = _vdm.clu_num, clu_total_num = _vdm.clu_total_num)
    
    if 'cluster_0' in vdm.query.getTitle():
        max_clu_num = vdm.clu_num

    vdm.id = vdm_id
    vdm.clu_rank = int(vdm.query.getTitle().split('_')[2]) 
    vdm.max_clu_num = max_clu_num
    all_2ndshell_vdms.append(vdm)
    all_2ndshell_vdm_id_dict[vdm.query.getTitle()] = vdm_id


### Creat allInOne2ndShellVdm
allInOne2ndShellVdm = pr.AtomGroup('AllInOne2ndShellVdm')
coords = []
chids = []
names = []
resnames = []
resnums = []
last = 0
for _vdm in all_2ndshell_vdms:
    qs = _vdm.query.select('name N C CA O MG MN FE CO NI CU ZN')
    coords.extend(qs.getCoords())
    names.extend(qs.getNames())
    resnames.extend(qs.getResnames())
    resnums.extend([last, last, last, last, last + 1, last + 1, last+ 1, last + 1, last + 2])
    last += 3
    if len(qs) != 9:
        print('Coords len is wrong.' + _vdm.query.getTitle())
allInOne2ndShellVdm.setCoords(coords)
chids = ['A' for i in range(len(coords))]
allInOne2ndShellVdm.setNames(names)
allInOne2ndShellVdm.setResnums(resnums)
allInOne2ndShellVdm.setResnames(resnames)
allInOne2ndShellVdm.setChids(chids)
pr.writePDB(outdir + allInOne2ndShellVdm.getTitle() + '.pdb', allInOne2ndShellVdm)

### Output 
with open(outdir + 'all_2ndshell_vdms.pkl', 'wb') as f:
    pickle.dump(all_2ndshell_vdms, f)

with open(outdir + 'all_2ndshell_vdm_id_dict.pkl', 'wb') as f:
    pickle.dump(all_2ndshell_vdm_id_dict, f)

with open(outdir + 'allInOne2ndShellVdm.pkl', 'wb') as f:
    pickle.dump(allInOne2ndShellVdm, f)