'''
Special for porphyrin.
Design metal + ligand bindin protein. 
Need to run search_selfcenter function first. 
'''

from metalprot.search import search_selfcenter
from metalprot.basic import filter
import metalprot.basic.constant as constant
from metalprot.search import position_ligand, position_porphyrin
import pickle
import time
import prody as pr
import os


'''
python /mnt/e/GitHub_Design/Metalprot/scrips/porphyrin/run_position_porphyrin.py

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

workdir = '/mnt/e/DesignData/ligands/LigandBB/1a6k_hem/'

outdir = workdir + 'output_selfcenter/'

target_path = workdir + '1a6k.pdb'

win_filter = [93]


geometry_path = None
#geometry_path = workdir + 'tetrahydral_geo.pdb'

metal_metal_dist = 0.45

num_contact_vdms = [1]

allowed_aa_combinations = [['H']]
#allowed_aa_combinations = []

_filter = filter.Search_filter(filter_abple = False, filter_phipsi = True, max_phipsi_val = 25, 
    filter_vdm_score = True, min_vdm_score = 0, filter_vdm_count = True, min_vdm_clu_num = 200,
    after_search_filter_geometry = False, filter_based_geometry_structure = False, angle_tol = 15, aa_aa_tol = 0.3, aa_metal_tol = 0.2,
    pair_angle_range = [85, 130], pair_aa_aa_dist_range = [2.8, 4], pair_metal_aa_dist_range = None,
    after_search_filter_qt_clash = True, vdm_vdm_clash_dist = 2.7, vdm_bb_clash_dist = 2.2, 
    after_search_open_site_clash = False, open_site_dist = 3.0, 
    write_filtered_result = False, selfcenter_filter_member_phipsi=True)

ss =  search_selfcenter.Search_selfcenter(target_path,  outdir, all_querys, cluster_centroid_dict, query_all_metal, 
    num_contact_vdms, metal_metal_dist, win_filter, validateOriginStruct = True, search_filter= _filter, geometry_path = None,
    density_radius = 0.6, allowed_aa_combinations = allowed_aa_combinations)

#ss.run_selfcenter_search()
search_selfcenter.run_search_selfcenter(ss)

end_time = time.time()
print(end_time - start_time, "seconds")

target = pr.parsePDB(workdir + '1a6k.pdb')

vdm_best = None
min_dist = 10
for key in ss.neighbor_comb_dict.keys():
    if 'H' not in str(key):
        continue
    dist = pr.calcDistance(target.select('name FE').getCoords(), list(ss.neighbor_comb_dict[key].centroid_dict.values())[0].get_metal_coord())
    if dist[0] < min_dist:
        min_dist = dist[0]
        vdm_best = list(ss.neighbor_comb_dict[key].centroid_dict.values())[0]


lig_path = workdir + 'hem_n.pdb'
lig = pr.parsePDB(lig_path)
lig_connects = ['NE2', 'FE']
ro1 = ['NA', 'FE']
ro2 = ['NE2', 'FE']


all_ligs = position_porphyrin.generate_rotated_porphyrins(lig, ro1, ro2, rotation_degree1 = 1, rotation_degree2 = 10)
position_porphyrin.porphyrin2vdm(all_ligs, lig_connects, vdm_best)


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

filtered_ligs = position_ligand.ligand_clashing_filter(all_ligs, ss.target, dist = 2.7)

len(filtered_ligs)

os.makedirs(ss.workdir + 'filtered_ligs/', exist_ok=True)
for lg in filtered_ligs:
    pr.writePDB(ss.workdir + 'filtered_ligs/' +  lg.getTitle(), lg)