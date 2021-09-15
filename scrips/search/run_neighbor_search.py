import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.search import search, search_eval
import pickle

'''
python /mnt/e/GitHub_Design/Metalprot/metalprot/scrips/run_neighbor_search.py
'''

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/20210826/pickle/'


with open(query_dir + 'AllMetalQuery.pkl', 'rb') as f:
    query_all_metal = pickle.load(f)

with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)

with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
    cluster_centroid_dict = pickle.load(f)

with open(query_dir + 'id_cluster_dict.pkl', 'rb') as f:
    id_cluster_dict = pickle.load(f)

with open(query_dir + 'cluster_centroid_origin_dict.pkl', 'rb') as f:
    cluster_centroid_origin_dict = pickle.load(f)

print(len(all_querys))


### run Search_struct
'''
workdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/'

outdir = workdir + 'output_neighbor_eval2/'

os.makedirs(outdir, exist_ok=True)

target_path = workdir + '5od1_zn.pdb'
'''

workdir = '/mnt/e/DesignData/ligands/LigandBB/2cab/'

outdir = workdir + 'output_neighbor_eval2/'

os.makedirs(outdir, exist_ok=True)

target_path = workdir + '2cab.pdb'

rmsd_cuts = 0.45

num_iters = [3]

win_filter = None


#win_filter = [30, 31, 32, 33, 34, 35, 36, 37, 38, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
#win_filter = [34,  60,  64]

#ss =  search.Search_vdM(target_path, outdir, all_querys, id_cluster_dict, cluster_centroid_dict, query_all_metal, cluster_centroid_origin_dict, num_iters, rmsd_cuts, win_filter, validateOriginStruct = True, filter_abple = False, filter_phipsi = True, filter_phipsi_val = 15, parallel = False)
#ss.run_neighbor_search()

ss =  search_eval.Search_eval(target_path, outdir, all_querys, id_cluster_dict, cluster_centroid_dict, query_all_metal, cluster_centroid_origin_dict, num_iters, rmsd_cuts, win_filter, validateOriginStruct = True, filter_abple = False, filter_phipsi = True, filter_phipsi_val = 15, parallel = False)
ss.eval_search()
#wins, combs = ss.eval_get_comb()
