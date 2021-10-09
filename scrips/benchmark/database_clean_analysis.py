'''
The old database didn't consider bfactor, some weird vdMs may come from cores that contain large bfactors.
The reason why we focus on analyze 3HIScontact core is because we assume there is some effect of the two oxygen in ASP/GLU.
'''

import os
import sys
import shutil
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
from metalprot.database import database_evaluate
from metalprot.basic import get_metal_contact_atoms
import itertools
import numpy as np
from numpy.core.fromnumeric import argmin


root_dir = '/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20210624/'

workdir = root_dir + '_Seq_core_date_3HIScontact/'

pdbs = database_extract.get_all_pbd_prody(workdir)


filtered_pdbs = []

cut_off = 45

for pdb in pdbs:
    cts = get_metal_contact_atoms(pdb)

    bs = []
    for c in cts:
        bs.append(c.getBeta())

    if all([b <= cut_off for b in bs]):
        filtered_pdbs.append(pdb)


### copy the filter pdb into a new folder.
outdir = root_dir + '_Seq_core_date_3HIScontact_B45/'
os.makedirs(outdir, exist_ok=True)
filtered_pdbs_dict = set([pdb.getTitle() for pdb in filtered_pdbs])
for p in os.listdir(workdir):
    if '.pdb' not in p:
        continue
    if p.split('.')[0] in filtered_pdbs_dict:
        shutil.copyfile(workdir + p, outdir + p)

dist_angle = []

for pdb in filtered_pdbs:
    pair_infos = database_evaluate.calc_pair_geometry(pdb)
    for info in pair_infos:
        dist_angle.append((info[0][0], info[0][1], info[1], pdb.getTitle(), info[2], info[3]))


with open(workdir + '_pairwise_angle_dist_info_filteredB45.tsv', 'w') as f:
    f.write('aa\tinds\tTitle\tDist\tAngle\n')
    for da in dist_angle:
        x = da[0] + '-' + da[1]
        if da[0] < da[1]:
            x = da[1] + '-' + da[0]
        f.write(x + '\t' + '\t'.join([str(x) for x in da[2:]]) + '\n')