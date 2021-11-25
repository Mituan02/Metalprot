from typing import Tuple
from metalprot.search import search_selfcenter
from metalprot.basic import filter
import metalprot.basic.constant as constant
from metalprot.search import position_ligand
import pickle
import time
import prody as pr
import os

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/search_selfcenter/run_selfcenter_search.py

'''
start_time = time.time()

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/pickle_all/'

with open(query_dir + 'all_metal_vdm.pkl', 'rb') as f:
    query_all_metal = pickle.load(f)

with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)

with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
    cluster_centroid_dict = pickle.load(f)

print(len(all_querys))


### run Search_struct

workdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/'

outdir = workdir + 'output_selfcenter/'

target_path = workdir + '5od1_zn.pdb'

#win_filter = []
win_filter = [35,  61,  65]



metal_metal_dist = 0.45

num_contact_vdms = [3]

allowed_aa_combinations = [['H', 'H', 'H']]

_filter = filter.Search_filter(filter_abple = False, filter_phipsi = True, max_phipsi_val = 25, 
    filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
    after_search_filter = True, pair_angle_range = [92, 116], pair_aa_aa_dist_range = [3.0, 3.5], pair_metal_aa_dist_range = None,
    filter_qt_clash = True, write_filtered_result = False, selfcenter_filter_member_phipsi=True)

ss =  search_selfcenter.Search_selfcenter(target_path, outdir, all_querys, cluster_centroid_dict, query_all_metal, 
    num_contact_vdms, metal_metal_dist, win_filter, validateOriginStruct = True, search_filter= _filter, density_radius = 0.6,
    allowed_aa_combinations = allowed_aa_combinations)

search_selfcenter.run_search_selfcenter(ss)

end_time = time.time()
print(end_time - start_time, "seconds")


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


lig_path = workdir + 'AG4_ZN.pdb'
lig = pr.parsePDB(lig_path)
lig_connects = ['N9', 'ZN']
ro1 = ['ZN', 'N9']
ro2 = ['N9', 'S8']


key = list(ss.best_aa_comb_dict.keys())[0]


ideal_geo = ss.best_aa_comb_dict[key].ideal_geo

ideal_geo_o = constant.tetrahydra_geo_o
min_geo_struct, min_rmsd = filter.Search_filter.get_min_geo(ideal_geo, ideal_geo_o) 
pr.writePDB(ss.workdir + min_geo_struct.getTitle(), min_geo_struct)

all_ligs = position_ligand.generate_rotated_ligs(lig, ro1, ro2, rotation_degree = 5)

position_ligand.lig_2_ideageo(all_ligs, lig_connects, min_geo_struct)

for lg in all_ligs:
    os.makedirs(ss.workdir + 'all_ligs/', exist_ok=True)
    pr.writePDB(ss.workdir + 'all_ligs/' +  lg.getTitle(), lg)

'''
nature_lig = pr.parsePDB(workdir + 'AG4_ZN.pdb')
min_RMSD = 100
min_lg = None
for lg in all_ligs:
    rmsd = pr.calcRMSD(lg, nature_lig)
    if rmsd < min_RMSD:
        min_RMSD = rmsd
        min_lg = lg
print(min_RMSD)
print(min_lg.getTitle())
#3.065111315604363
#ff3_zn_120_ZN-N3_50_N3-S2_340

pr.writePDB(workdir + '_min_' + min_lg.getTitle(), min_lg)
'''

filtered_ligs = position_ligand.ligand_clashing_filter(all_ligs, ss.target, dist = 2.5)

len(filtered_ligs)

for lg in filtered_ligs:
    os.makedirs(ss.workdir + 'filtered_ligs/', exist_ok=True)
    pr.writePDB(ss.workdir + 'filtered_ligs/' +  lg.getTitle(), lg)