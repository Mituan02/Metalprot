'''
Design metal + ligand bindin protein.
Need to run search_selfcenter function first. 
'''

from metalprot.search import search_selfcenter
from metalprot.basic import filter
import metalprot.basic.constant as constant
from metalprot.combs import position_ligand
import pickle
import time
import prody as pr
import os


#ligand positioning for the first represents.
# lig_path = workdir + 'Zn_Phenol.pdb'
# lig = pr.parsePDB(lig_path)
# lig_connects = ['OH', 'ZN']
# ro1 = ['ZN', 'OH']
# ro2 = ['OH', 'CZ']


# lig_path = workdir + 'ff3_zn_120.pdb'
# lig = pr.parsePDB(lig_path)
# lig_connects = ['N3', 'ZN']
# ro1 = ['ZN', 'N3']
# ro2 = ['N3', 'S2']


# lig_path = workdir + 'AG4_ZN.pdb'
# lig = pr.parsePDB(lig_path)
# lig_connects = ['N9', 'ZN']
# ro1 = ['ZN', 'N9']
# ro2 = ['N9', 'S8']


lig_path = workdir + 'anh_zn.pdb'
lig = pr.parsePDB(lig_path)
lig_connects = ['NE2', 'ZN']
ro1 = ['ZN', 'NE2']
ro2 = ['NE2', 'CE1']


key = list(ss.best_aa_comb_dict.keys())[0]


ideal_geo = ss.best_aa_comb_dict[key].ideal_geo

ideal_geo_o = constant.tetrahydra_geo_o
min_geo_struct, min_rmsd = filter.Search_filter.get_min_geo(ideal_geo, ideal_geo_o) 
pr.writePDB(ss.workdir + min_geo_struct.getTitle(), min_geo_struct)

all_ligs = position_ligand.generate_rotated_ligs(lig, ro1, ro2, rotation_degree = 10)

position_ligand.lig_2_ideageo(all_ligs, lig_connects, min_geo_struct)


os.makedirs(ss.workdir + 'all_ligs/', exist_ok=True)
for lg in all_ligs:
    pr.writePDB(ss.workdir + 'all_ligs/' +  lg.getTitle(), lg)

'''
nature_lig = pr.parsePDB(lig_path)
min_RMSD = 100
min_lg = None
for lg in all_ligs:
    rmsd = pr.calcRMSD(lg, nature_lig)
    if rmsd < min_RMSD:
        min_RMSD = rmsd
        min_lg = lg
print(min_RMSD)
print(min_lg.getTitle())


pr.writePDB(workdir + '_min_' + min_lg.getTitle(), min_lg)
'''

filtered_ligs = position_ligand.ligand_clashing_filter(all_ligs, ss.target, dist = 3)

len(filtered_ligs)

os.makedirs(ss.workdir + 'filtered_ligs/', exist_ok=True)
for lg in filtered_ligs:
    pr.writePDB(ss.workdir + 'filtered_ligs/' +  lg.getTitle(), lg)