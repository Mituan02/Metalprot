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
import importlib.machinery

'''
# To run in develop mode
python run_selfcenter_search.py 0 /mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta_86-88-101/_rosetta_r3/eval/run_selfcenter_search_paras.py 1
'''


def run(target_path, query_dir, outdir, win_filter, para):

    with open(query_dir + 'all_metal_vdm.pkl', 'rb') as f:
        query_all_metal = pickle.load(f)

    with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
        all_querys = pickle.load(f)

    with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
        cluster_centroid_dict = pickle.load(f)

    print(len(all_querys))

    ### run Search_struct

    geometry_path = para.geometry_path
    metal_metal_dist = para.metal_metal_dist
    num_contact_vdms = para.num_contact_vdms
    allowed_aa_combinations = para.allowed_aa_combinations


    _filter = filter.Search_filter(filter_abple = False, filter_phipsi = False, max_phipsi_val = 30, 
        filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
        after_search_filter_geometry = True, filter_based_geometry_structure = True, angle_tol = 35, aa_aa_tol = 0.35, aa_metal_tol = 0.25,
        pair_angle_range = [92, 116], pair_aa_aa_dist_range = [3.0, 3.5], pair_metal_aa_dist_range = None,
        after_search_filter_qt_clash = True, vdm_vdm_clash_dist = 2.7, vdm_bb_clash_dist = 2.2, 
        write_filtered_result = False, selfcenter_filter_member_phipsi=True)

    ss =  search_selfcenter.Search_selfcenter(target_path,  outdir, all_querys, cluster_centroid_dict, query_all_metal, 
        num_contact_vdms, metal_metal_dist, win_filter, validateOriginStruct = False, search_filter= _filter, geometry_path = geometry_path,
        density_radius = 0.6, allowed_aa_combinations = allowed_aa_combinations)

    search_selfcenter.run_search_selfcenter(ss)
    return 


def main():
    on_wynton = bool(int(sys.argv[1]))
    path = sys.argv[2]
    para = importlib.machinery.SourceFileLoader('para', path).load_module()
    print(path)
    query_dir, workdir = para.get_file_path(on_wynton)
    print('on_wynton: ' + str(on_wynton))

    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    ind = int(sys.argv[3]) -1
    if ind > len(pdb_files) -1:
        return
    target_path = workdir + pdb_files[ind]
    win_filter = para.win_filters[ind]
    outdir = workdir + 'output_' + pdb_files[ind].split('.')[0] + '_/'
    
    run(target_path, query_dir, outdir, win_filter, para)
    return


if __name__=='__main__':
    main()