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
#To run in develop mode 
python run_selfcenter_search.py
'''

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/pickle_noCYS/'

workdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/'

win_filters = [[34,  60,  64] , [4, 6, 15]]

def main():
    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])
    ind = int(sys.argv[1]) -1
    target_path = workdir + pdb_files[ind]
    win_filter = win_filters[ind]
    outdir = workdir + 'output_selfcenter_' + pdb_files[ind].split('.')[0] + '_/'
    
    run(target_path, win_filter, outdir)
    return

def run(target_path, win_filter, outdir):

    with open(query_dir + 'all_metal_vdm.pkl', 'rb') as f:
        query_all_metal = pickle.load(f)

    with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
        all_querys = pickle.load(f)

    with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
        cluster_centroid_dict = pickle.load(f)

    print(len(all_querys))

    ### run Search_struct

    rmsd_cuts = 0.45

    num_iters = [3]

    _filter = filter.Search_filter(filter_abple = False, filter_phipsi = True, max_phipsi_val = 25, 
        filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
        after_search_filter = True, pair_angle_range = [92, 116], pair_aa_aa_dist_range = [3.0, 3.5], pair_metal_aa_dist_range = None,
        filter_qt_clash = True, write_filtered_result = False, selfcenter_filter_member_phipsi=True)

    ss =  search_selfcenter.Search_selfcenter(target_path, outdir, all_querys, cluster_centroid_dict, query_all_metal, 
        num_iters, rmsd_cuts, win_filter, validateOriginStruct = True, search_filter= _filter, parallel = False)

    ss.run_selfcenter_search()

if __name__=='__main__':
    main()