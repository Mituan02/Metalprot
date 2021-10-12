import os
from re import T
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.search import search_selfcenter
from metalprot.basic import filter
import pickle

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/search_selfcenter/run_selfcenter_search.py
'''

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211008/20211009_selfcenter/pickle_noCYS_alignBB/'


with open(query_dir + 'all_metal_vdm.pkl', 'rb') as f:
    all_metal_vdm = pickle.load(f)

with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_vdms = pickle.load(f)

with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
    cluster_centroid_dict = pickle.load(f)



print(len(all_vdms))


### run Search_struct

workdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/'

outdir = workdir + 'output_selfcenter/'

target_path = workdir + '5od1_zn.pdb'


# workdir = '/mnt/e/DesignData/ligands/LigandBB/2cab/'

# outdir = workdir + 'output_selfcenter_test3/'

# target_path = workdir + '2cab.pdb'

rmsd_cuts = 0.45

num_iters = [3]

win_filter = None


#win_filter = [30, 31, 32, 33, 34, 35, 36, 37, 38, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
win_filter = [34,  60,  64]
#win_filter = [89,  91,  114]

_filter = filter.Search_filter(filter_abple = False, filter_phipsi = True, max_phipsi_val = 15, 
    filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
    after_search_filter = True, pair_angle_range = [92, 116], pair_aa_aa_dist_range = [3.0, 3.5], pair_metal_aa_dist_range = None,
    filter_qt_clash = True, write_filtered_result = True)

ss =  search_selfcenter.Search_selfcenter(target_path, outdir, all_vdms, cluster_centroid_dict, all_metal_vdm, 
    num_iters, rmsd_cuts, win_filter, validateOriginStruct = True, search_filter= _filter, parallel = False)

ss.run_selfcenter_search()
#ss.run_selfcenter_search2()

