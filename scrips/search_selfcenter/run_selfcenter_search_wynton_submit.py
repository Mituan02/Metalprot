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

win_filters = [[] , []]

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

    metal_metal_dist = 0.45

    num_contact_vdms = [3]

    allowed_aa_combinations = [['H', 'H', 'H']]

    geometry_path = None
    #geometry_path = workdir + 'tetrahydral_geo.pdb'

    _filter = filter.Search_filter(filter_abple = False, filter_phipsi = True, max_phipsi_val = 25, 
        filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
        after_search_filter_geometry = True, filter_based_geometry_structure = False, angle_tol = 15, aa_aa_tol = 0.3, aa_metal_tol = 0.2,
        pair_angle_range = [92, 116], pair_aa_aa_dist_range = [3.0, 3.5], pair_metal_aa_dist_range = None,
        after_search_filter_qt_clash = True, vdm_vdm_clash_dist = 2.7, vdm_bb_clash_dist = 2.2, 
        write_filtered_result = False, selfcenter_filter_member_phipsi=True)

    ss =  search_selfcenter.Search_selfcenter(target_path,  outdir, all_querys, cluster_centroid_dict, query_all_metal, 
        num_contact_vdms, metal_metal_dist, win_filter, validateOriginStruct = True, search_filter= _filter, geometry_path = None,
        density_radius = 0.6, allowed_aa_combinations = allowed_aa_combinations)

    search_selfcenter.run_search_selfcenter(ss)

if __name__=='__main__':
    main()