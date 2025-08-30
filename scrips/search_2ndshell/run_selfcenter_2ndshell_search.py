import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.search import search_selfcenter
from metalprot.basic import filter
import pickle

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/search_2ndshell/run_selfcenter_2ndshell_search.py
'''

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/pickle_noCYS/'

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

win_filter = [35,  61,  65]


# workdir = '/mnt/e/DesignData/ligands/LigandBB/6dwv/'

# outdir = workdir + 'output_selfcenter/'

# target_path = workdir + '6dwv_core.pdb'

# win_filter = [4, 6, 15]


geometry_path = None
#geometry_path = workdir + 'tetrahydral_geo.pdb'

metal_metal_dist = 0.3

num_contact_vdms = [3]

allowed_aa_combinations = [['H', 'H', 'H']]
allowed_aa_combinations = []

_filter = filter.Search_filter(filter_abple = False, filter_phipsi = True, max_phipsi_val = 25, 
    filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
    after_search_filter_geometry = True, filter_based_geometry_structure = False, angle_tol = 15, aa_aa_tol = 0.3, aa_metal_tol = 0.2,
    pair_angle_range = [85, 130], pair_aa_aa_dist_range = [2.8, 4], pair_metal_aa_dist_range = None,
    after_search_filter_qt_clash = True, vdm_vdm_clash_dist = 2.7, vdm_bb_clash_dist = 2.2, 
    after_search_open_site_clash = True, open_site_dist = 3.0, 
    write_filtered_result = False, selfcenter_filter_member_phipsi=True)


ss =  search_selfcenter.Search_selfcenter(target_path,  outdir, all_querys, cluster_centroid_dict, query_all_metal, 
    num_contact_vdms, metal_metal_dist, win_filter, validateOriginStruct = True, search_filter= _filter, geometry_path = None,
    density_radius = 0.6, allowed_aa_combinations = allowed_aa_combinations, output_wincomb_overlap=True)

#ss.run_selfcenter_search()
search_selfcenter.run_search_selfcenter(ss)

from metalprot.search import __search_2ndshell

_2nd_query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/20220116_selfcenter_bb2ndshell_notconnect/pickle_noCYS/'

with open(_2nd_query_dir + 'all_2ndshell_vdms.pkl', 'rb') as f:
    all_2ndshell_vdms = pickle.load(f)
with open(_2nd_query_dir + 'allInOne2ndShellVdm.pkl', 'rb') as f:
    allInOne2ndShellVdm = pickle.load(f)

print(len(all_2ndshell_vdms))

ss.secondshell_vdms = all_2ndshell_vdms
len(list(ss.best_aa_comb_dict))

__search_2ndshell.run_search_2ndshell(ss.best_aa_comb_dict, ss.target, ss.secondshell_vdms, allInOne2ndShellVdm, rmsd_2ndshell = 1.0)

__search_2ndshell.write_2ndshell(ss, ss.workdir, ss.best_aa_comb_dict)

