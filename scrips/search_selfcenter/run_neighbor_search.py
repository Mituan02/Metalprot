import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.search import search_selfcenter
import pickle

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/search_selfcenter/run_neighbor_search.py
'''

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/20210911/pickle/'


with open(query_dir + 'AllMetalQuery.pkl', 'rb') as f:
    query_all_metal = pickle.load(f)

with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)

with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
    cluster_centroid_dict = pickle.load(f)

with open(query_dir + 'id_cluster_dict.pkl', 'rb') as f:
    id_cluster_dict = pickle.load(f)

cluster_centroid_origin_dict = None 

print(len(all_querys))


### run Search_struct

workdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/'

outdir = workdir + 'output_selfcenter_test/'

os.makedirs(outdir, exist_ok=True)

target_path = workdir + '5od1_zn.pdb'

rmsd_cuts = 0.45

num_iters = [3]

win_filter = None


#win_filter = [30, 31, 32, 33, 34, 35, 36, 37, 38, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
win_filter = [34,  60,  64]

ss =  search_selfcenter.Search_selfcenter(target_path, outdir, all_querys, id_cluster_dict, cluster_centroid_dict, query_all_metal, cluster_centroid_origin_dict, num_iters, rmsd_cuts, win_filter, validateOriginStruct = True, filter_abple = False, filter_phipsi = True, filter_phipsi_val =15, parallel = False)

#ss.run_neighbor_search()

ss.neighbor_generate_query_dict()

ss.neighbor_generate_pair_dict()

win_combs = ss.neighbor_get_win_combs()

comb_dict = ss.neighbor_construct_comb(win_combs[0])

ss.neighbor_extract_query(comb_dict)
ss.neighbor_calc_geometry(comb_dict)

ss.comb_overlap(comb_dict)
ss.neighbor_write_win(comb_dict)

#ss.eval_search()

