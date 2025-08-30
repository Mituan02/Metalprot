import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import search
import pickle
import multiprocessing as mp

'''
python /mnt/e/GitHub_Design/Metalprot/metalprot/scrips/run_neighbor_search_all_structs.py
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


# run Search_struct
def run_search(workdir, target_file, query_all_metal, all_querys, cluster_centroid_dict, id_cluster_dict):

    #workdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/'

    outdir = workdir + 'output_neighbor_val_' + target_file + '/'

    target_path = workdir + target_file

    print(target_path)

    rmsd_cuts = 0.25

    num_iters = [3, 4]

    win_filter = None


    ss =  search.Search_vdM(target_path, outdir, all_querys, id_cluster_dict, cluster_centroid_dict, query_all_metal, num_iters, rmsd_cuts, win_filter, validateOriginStruct = True, filter_abple = True)
      
    win_search = set()

    try:
        ss.run_neighbor_search()
    except:
        return (target_file + ' Error', win_search)

    for k in ss.neighbor_comb_dict.keys():
        win_search.add(k[0])
    
    return (target_file, win_search)


num_cores = int(mp.cpu_count())
pool = mp.Pool(num_cores)

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/_test_full_pdbs_sub/'

target_files = []
for target_file in os.listdir(workdir):
    if target_file.endswith('.pdb'):
        target_files.append(target_file)

results = [pool.apply_async(run_search, args=(workdir, target_file, query_all_metal, all_querys, cluster_centroid_dict, id_cluster_dict)) for target_file in target_files]
results = [p.get() for p in results]

with open(workdir + '_summary.txt', 'w') as f:
    f.write('target_file\twin_extract\n')
    for r in results:
        try:
            f.write(r[0] + '\t')
            f.write(str(r[1]) + '\t')
        except:
            f.write(r[0] + '\t\n')

