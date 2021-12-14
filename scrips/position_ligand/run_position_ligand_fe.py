'''
Design metal + ligand bindin protein.
Need to run search_selfcenter function first. 
'''

from metalprot.search import search_selfcenter
from metalprot.basic import filter
import metalprot.basic.constant as constant
from metalprot.search import position_ligand
import pickle
import time
import prody as pr
import os


start_time = time.time()

query_dir = '/mnt/e/DesignData/ligands/FE_rcsb/20211206/20211206_selfcenter/pickle_all/'

with open(query_dir + 'all_metal_vdm.pkl', 'rb') as f:
    query_all_metal = pickle.load(f)

with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)

with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
    cluster_centroid_dict = pickle.load(f)

print(len(all_querys))


### run Search_struct

workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/ntf2/'

outdir = workdir + 'output_selfcenter/'

target_path = workdir + '1cqs.pdb'

win_filter = [16, 17, 84]

metal_metal_dist = 0.6

num_contact_vdms = [3]

allowed_aa_combinations = [['H', 'H', 'D'], ['H', 'H', 'E']] 
#allowed_aa_combinations = []

geometry_path = None
geometry_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/fe_geo.pdb'

_filter = filter.Search_filter(filter_abple = False, filter_phipsi = True, max_phipsi_val = 30, 
    filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
    after_search_filter_geometry = True, filter_based_geometry_structure = True, angle_tol = 15, aa_aa_tol = 0.35, aa_metal_tol = 0.25,
    pair_angle_range = [92, 116], pair_aa_aa_dist_range = [3.0, 3.5], pair_metal_aa_dist_range = None,
    after_search_filter_qt_clash = True, vdm_vdm_clash_dist = 2.7, vdm_bb_clash_dist = 2.2, 
    write_filtered_result = False, selfcenter_filter_member_phipsi=True)

ss =  search_selfcenter.Search_selfcenter(target_path,  outdir, all_querys, cluster_centroid_dict, query_all_metal, 
    num_contact_vdms, metal_metal_dist, win_filter, validateOriginStruct = False, search_filter= _filter, geometry_path = geometry_path,
    density_radius = 0.6, allowed_aa_combinations = allowed_aa_combinations)


search_selfcenter.run_search_selfcenter(ss)

########################################################################################

dbdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/'
lig_path = dbdir + 'tts_fe.pdb'
lig = pr.parsePDB(lig_path)
lig_connects = ['O1','O3', 'FE1']
ro1 = ['C8', 'C7']
rest1 = ['C9', 'O1', 'O3', 'O4', 'FE1']
ro2 = ['C7', 'C6']
rest2 = ['H6', 'H7', 'C8', 'C9', 'O1', 'O3', 'O4', 'FE1']


all_ligs = position_ligand.generate_rotated_ligs(lig, [ro1, ro2], [rest1, rest2], [20, 20])


key = list(ss.best_aa_comb_dict.keys())[0]
ideal_geo = ss.best_aa_comb_dict[key].ideal_geo
#ideal_geo = pr.parsePDB(dbdir + 'fe_geo.pdb')
ideal_geo_o = pr.parsePDB(dbdir + 'fe_geo_o.pdb')

min_geo_struct, min_rmsd = filter.Search_filter.get_min_geo(ideal_geo, ideal_geo_o) 
pr.writePDB(ss.workdir + min_geo_struct.getTitle(), min_geo_struct)

position_ligand.lig_2_ideageo(all_ligs, lig_connects, min_geo_struct, geo_sel = 'name OE2 OK1 FE')

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

filtered_ligs = position_ligand.ligand_clashing_filter(all_ligs, ss.target, dist = 2.5)

len(filtered_ligs)

os.makedirs(ss.workdir + 'filtered_ligs/', exist_ok=True)
for lg in filtered_ligs:
    pr.writePDB(ss.workdir + 'filtered_ligs/' +  lg.getTitle(), lg)