import os
import sys

from metalprot.basic.quco import Query, Comb
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract

'''
Extract binary metal binding dists and angles from bivalence vdMs. 
'''

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/AAdAAMetal_aa_aa_reps/"

pdbs = database_extract.get_all_pbd_prody(workdir)

dist_angle = []

for pdb in pdbs:
    try:
        query = Query(pdb, is_bivalent=True)
        comb = Comb([query])
        comb.calc_pair_geometry()
        if not comb.pair_dists:
            dist_angle.append((comb.querys[0].query.select('name CA').getResnames()[0], comb.querys[0].query.select('name CA').getResnames()[1], comb.querys[0].query.getTitle(), None, None))
        else:
            dist_angle.append((comb.querys[0].query.select('name CA').getResnames()[0], comb.querys[0].query.select('name CA').getResnames()[1], comb.querys[0].query.getTitle(), comb.pair_dists[0], comb.pair_angles[0]))
    except:
        print(pdb.getTitle())

with open(workdir + '_pairwise_angle_dist_info.tsv', 'w') as f:
    f.write('aa\tTitle\tDist\tAngle\n')
    for da in dist_angle:
        x = da[0] + '-' + da[1]
        if da[0] < da[1]:
            x = da[1] + '-' + da[0]
        f.write(x + '\t' + '\t'.join([str(x) for x in da[2:]]) + '\n')


