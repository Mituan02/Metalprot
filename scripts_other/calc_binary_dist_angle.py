import os
import sys
from metalprot.apps.search_struct import Query, Comb
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import ligand_database

'''
Extract binary metal binding dists and angles. 
'''

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb/20210527/M5_2aa_sep_cores_reps/"

pdbs = ligand_database.get_all_pbd_prody(workdir)

dist_angle = []
count = 0
for pdb in pdbs:
    count += 1
    if count%2 == 0:
        continue
    query = Query(pdb, is_bivalent=True)
    comb = Comb([query])
    comb.calc_pair_geometry()
    if not comb.pair_dists:
        dist_angle.append((comb.querys[0].query.getTitle(), None, None))
    else:
        dist_angle.append((comb.querys[0].query.getTitle(), comb.pair_dists[0], comb.pair_angles[0]))


with open(workdir + '_summary.txt', 'w') as f:
    f.write('Title\tDist\tAngle\n')
    for da in dist_angle:
        f.write('\t'.join([str(x) for x in da]) + '\n')


