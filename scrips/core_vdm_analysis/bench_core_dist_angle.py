'''
Calculate the pairwise infomation of the 3-contact binding (HIS, ASP, GLU) core pdbs in _Seq_core_date_3contact
'''


import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract, database_evaluate
import itertools
import numpy as np


'''
Extract binary metal binding dists and angles. 
'''
metal_sel = 'name NI MN ZN CO CU MG FE' 

root_dir = '/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20210624/'

workdir = root_dir + '_Seq_core_date_reps/'

outdir = root_dir + '_Seq_core_date_3contact/'

workdir = outdir

pdbs = database_extract.get_all_pbd_prody(workdir)

dist_angle = []

for pdb in pdbs:
    pair_infos = database_evaluate.calc_pair_geometry(pdb)
    for info in pair_infos:
        dist_angle.append((info[0][0], info[0][1], info[1], pdb.getTitle(), info[2], info[3]))


with open(workdir + '_pairwise_angle_dist_info.tsv', 'w') as f:
    f.write('aa\tinds\tTitle\tDist\tAngle\n')
    for da in dist_angle:
        x = da[0] + '-' + da[1]
        if da[0] < da[1]:
            x = da[1] + '-' + da[0]
        f.write(x + '\t' + '\t'.join([str(x) for x in da[2:]]) + '\n')

