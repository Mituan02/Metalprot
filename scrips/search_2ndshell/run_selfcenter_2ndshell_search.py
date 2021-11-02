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

win_filter = [34,  60,  64]


# workdir = '/mnt/e/DesignData/ligands/LigandBB/6dwv/'

# outdir = workdir + 'output_selfcenter/'

# target_path = workdir + '6dwv_core.pdb'

# win_filter = [4, 6, 15]


rmsd_cuts = 0.45

num_iters = [3]


_filter = filter.Search_filter(filter_abple = False, filter_phipsi = True, max_phipsi_val = 15, 
    filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
    after_search_filter = True, pair_angle_range = [92, 116], pair_aa_aa_dist_range = [3.0, 3.5], pair_metal_aa_dist_range = None,
    filter_qt_clash = True, write_filtered_result = False, selfcenter_filter_member_phipsi=True)

ss =  search_selfcenter.Search_selfcenter(target_path, outdir, all_querys, cluster_centroid_dict, query_all_metal, 
    num_iters, rmsd_cuts, win_filter, validateOriginStruct = True, search_filter= _filter, parallel = False)

ss.run_selfcenter_search()


from metalprot.basic import vdmer_2ndshell
from metalprot.search import search_2ndshell

_2nd_query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211018/20211026_2ndshell_selfcenter/pickle_noCYS/'

with open(_2nd_query_dir + 'all_2ndshell_vdms.pkl', 'rb') as f:
    all_2ndshell_vdms = pickle.load(f)

print(len(all_2ndshell_vdms))

ss.secondshell_vdms = all_2ndshell_vdms
len(list(ss.best_aa_comb_dict))

search_2ndshell.run_search_2ndshell(ss.best_aa_comb_dict, ss.target, ss.secondshell_vdms, ss.rmsd_2ndshell)

search_2ndshell.write_2ndshell(ss.workdir, ss.best_aa_comb_dict)
