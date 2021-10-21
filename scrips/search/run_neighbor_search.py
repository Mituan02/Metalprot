from math import fabs
import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.search import search, search_eval
from metalprot.basic import filter
import pickle

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/search/run_neighbor_search.py
'''

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_category/pickle_noCYS/'

with open(query_dir + 'all_metal_vdm.pkl', 'rb') as f:
    query_all_metal = pickle.load(f)

with open(query_dir + 'all_vdms.pkl', 'rb') as f:
    all_querys = pickle.load(f)

with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
    cluster_centroid_dict = pickle.load(f)

print(len(all_querys))


### run Search_struct

# workdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/'

# outdir = workdir + 'output_neighbor/'

# target_path = workdir + '5od1_zn.pdb'

#win_filter = [34,  60,  64]


# workdir = '/mnt/e/DesignData/ligands/LigandBB/2cab/'

# outdir = workdir + 'output_neighbor/'

# target_path = workdir + '2cab.pdb'

# rmsd_cuts = 0.45

# num_iters = [3]

# win_filter = [89, 91, 114]

workdir = '/mnt/e/DesignData/ligands/LigandBB/6dwv/'

outdir = workdir + 'output_neighbor/'

target_path = workdir + '6dwv_core.pdb'

rmsd_cuts = 0.45

num_iters = [3]

win_filter = [4, 6, 15]


_filter = filter.Search_filter(filter_abple = False, filter_phipsi = True, max_phipsi_val = 25, 
    filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
    after_search_filter = True, pair_angle_range = [92, 116], pair_aa_aa_dist_range = [3.0, 3.5], pair_metal_aa_dist_range = None,
    filter_qt_clash = True, write_filtered_result = True)

ss =  search.Search_vdM(target_path, outdir, all_querys, cluster_centroid_dict, 
    query_all_metal, num_iters, rmsd_cuts, 
    win_filter, validateOriginStruct = True, search_filter= _filter, parallel = False)

ss.run_neighbor_search()




