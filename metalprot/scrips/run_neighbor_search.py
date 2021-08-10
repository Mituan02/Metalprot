import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import extract_vdm
from metalprot import hull
from metalprot import search
import pickle

'''
python /mnt/e/GitHub_Design/Metalprot/metalprot/scrips/run_neighbor_search.py
'''

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/pickle/'


with open(query_dir + 'AllMetalQuery.pkl', 'rb') as f:
    query_all_metal = pickle.load(f)

with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)

with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
    cluster_centroid_dict = pickle.load(f)

with open(query_dir + 'id_cluster_dict.pkl', 'rb') as f:
    id_cluster_dict = pickle.load(f)

print(len(all_querys))


### run Search_struct

workdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/'

outdir = workdir + 'output_neighbor_val_test5/'

target_path = workdir + '5od1_zn.pdb'

rmsd_cuts = 0.25

num_iter = 3

win_filter = None

'''
#For testing natural metal binding proteins.
_pdb = pr.parsePDB(target_path)

cores = ligand_database.get_metal_core_seq(_pdb, metal_sel = 'ion or name NI MN ZN CO CU MG FE' , extend = 4)

print(cores[0][1].select('name CA').getSequence())

#pr.writePDB(target_path + '_core.pdb', cores[0][1])

win_filter = [x for x in np.unique(cores[0][1].getResindices())][0:-1]

print(win_filter)
'''

#win_filter = [30, 31, 32, 33, 34, 35, 36, 37, 38, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
win_filter = [34,  60,  64]

ss =  search.Search_vdM(target_path, outdir, all_querys, id_cluster_dict, cluster_centroid_dict, query_all_metal, num_iter, rmsd_cuts, win_filter, validateOriginStruct = True, filter_abple = True)

ss.run_neighbor_search()

'''
ss.neighbor_generate_query_dict()

ss.neighbor_generate_pair_dict()

ss.neighbor_search_wins()

ss.neighbor_extract_query()

ss.neighbor_calc_comb_score()

ss.neighbor_calc_geometry()ssssss

ss.neighbor_write()

ss.neighbor_write_summary()

'''


#Test Graph
'''
win_comb = [34, 60, 64]

graph = search.Graph(win_comb, len(ss.querys))

graph.calc_pair_connectivity(ss.neighbor_pair_dict)

graph.get_paths()

'''
