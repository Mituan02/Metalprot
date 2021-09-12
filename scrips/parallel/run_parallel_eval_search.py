'''
Usage:
    ./run_parallel_eval_search.py <workdir> <file> [options]

Arguments:
    <workdir>
        Name of the job.

    <file>
        Job script for the analysis.

Options:

'''


import os
import sys
import prody as pr
import numpy as np
from metalprot import search
from metalprot import search_struct, extract_vdm, ligand_database
import pickle
import multiprocessing as mp
import time
import docopt


def run_search_struct(workdir, pdb_file, all_querys, id_cluster_dict, cluster_centroid_dict, query_all_metal, cluster_centroid_origin_dict):

    print('Start search ' + pdb_file)

    start_time = time.time()

    outdir = workdir + 'output_dist45_phipsi30_' + pdb_file.split('.')[0]

    os.makedirs(outdir, exist_ok=True)

    target_path = workdir + pdb_file

    rmsd_cuts = 0.45

    num_iters = [3]

    win_filter = None

    ss =  search.Search_vdM(target_path, outdir, all_querys, id_cluster_dict, cluster_centroid_dict, query_all_metal, cluster_centroid_origin_dict, num_iters, rmsd_cuts, win_filter, validateOriginStruct = False, filter_abple = False, filter_phipsi = True, filter_phipsi_val = 30, parallel = False)

    #ss.run_neighbor_search()
    try:
        ss.eval_search()
        end_time = time.time()
    except:
        end_time = time.time()
        return ('Failed' + pdb_file, int(end_time - start_time))
    return (pdb_file, int(end_time - start_time))

'''
python /wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/data/zn_eval_bench_mark/run_neighbor_parallel_eval_search.py
'''
def main(arguments):

    workdir = arguments['<workdir>']

    pdb_file = arguments['<file>'] 


    print(workdir)
    print(pdb_file)

    query_dir = '/wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/data/pickle/'

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

    run_search_struct(workdir, pdb_file, all_querys, id_cluster_dict, cluster_centroid_dict, query_all_metal, cluster_centroid_origin_dict)

if __name__=='__main__':
    arguments = docopt.docopt(__doc__)
    main(arguments)


'''
### run parallel Search_struct

workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/data/zn_eval_bench_mark/2013_2014/'

target_files = []
for target_file in os.listdir(workdir):
    if target_file.endswith('.pdb'):
        summaryfile = workdir + 'output_dist45_phipsi30_' + target_file.split('.')[0] + '/_summary.tsv'
        if os.path.isfile(summaryfile):
            continue
        target_files.append(target_file)



#num_cores = int(mp.cpu_count()/2)
num_cores = 15
pool = mp.Pool(num_cores)
print('There are total {} files'.format(len(target_files)))
print(target_files[:15])

results = [pool.apply_async(run_search_struct, args=(workdir, target_file)) for target_file in target_files]

output = [p.get() for p in results]
'''
