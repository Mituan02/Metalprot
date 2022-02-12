'''
Design metal + ligand bindin protein.
For wynton submit. Please put all protein backbone in one folder and the folder cannot contain other .pdb files.
'''

from metalprot.search import search_selfcenter
from metalprot.basic import filter
import metalprot.basic.constant as constant
from metalprot.combs import position_ligand, search_ligand
import pickle
import time
import prody as pr
import os

import sys
# sys.path.append(r'/wynton/home/degradolab/lonelu/GitHub_Design/Combs2')
# import combs2

from metalprot.combs import search_cg_vdms as metalprot_scv
import pandas as pd
import os
pd.set_option("display.max_columns", None)

start_time = time.time()


def run_metal(query_dir, workdir, outdir, target_path, win_filter):
    '''
    Search possible metal binding positions.
    '''
    with open(query_dir + 'all_metal_vdm.pkl', 'rb') as f:
        query_all_metal = pickle.load(f)

    with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
        all_querys = pickle.load(f)

    with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
        cluster_centroid_dict = pickle.load(f)


    ### run Search_struct

    allowed_aa_combinations = [['H', 'H', 'D'], ['H', 'H', 'E']] 
    #allowed_aa_combinations = []

    geometry_path = None
    geometry_path = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe/fe_geo.pdb'
    geometry_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/fe_geo.pdb'

    _filter = filter.Search_filter(filter_abple = False, filter_phipsi = True, max_phipsi_val = 30, 
        filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
        after_search_filter_geometry = True, filter_based_geometry_structure = True, angle_tol = 15, aa_aa_tol = 0.35, aa_metal_tol = 0.25,
        pair_angle_range = [92, 116], pair_aa_aa_dist_range = [3.0, 3.5], pair_metal_aa_dist_range = None,
        after_search_filter_qt_clash = True, vdm_vdm_clash_dist = 2.7, vdm_bb_clash_dist = 2.2, 
        write_filtered_result = False, selfcenter_filter_member_phipsi=True)

    ss =  search_selfcenter.Search_selfcenter(target_path,  outdir, all_querys, cluster_centroid_dict, query_all_metal, 
        num_contact_vdms = [3], metal_metal_dist = 0.6, win_filtered = win_filter, validateOriginStruct = False, search_filter= _filter, geometry_path = geometry_path,
        density_radius = 0.6, allowed_aa_combinations = allowed_aa_combinations)


    search_selfcenter.run_search_selfcenter(ss)
    ss.write_for_combs()

    print('Search metal finished.')
    return ss


def run_ligand(ss, lig, ro1, ro2, rest1, rest2, ideal_geo_o_path, clash_dist):
    '''
    Generate all potential ligands for each binding position. 
    '''
    all_ligs = position_ligand.generate_rotated_ligs(lig, [ro1, ro2], [rest1, rest2], [5, 5])

    search_ligands = search_ligand.prepare_search_ligand(ss.workdir + 'represents_combs/', ideal_geo_o_path)

    for sl in search_ligands:
        sl.generate_ligands(all_ligs, ss.target, lig_connects,  clash_dist = clash_dist)
        sl.write_ligands(write_all_ligands = False)
    return search_ligands


def load_vdms(sl, resnums, path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/vdMs/'):
    '''
    Load vdMs on the current protein target.
    '''
    input_dir = sl.outdir 
    os.makedirs(input_dir, exist_ok = True)

    pdb_gly = sl.all_gly
    pdb_ala = sl.all_ala

    template = combs2.design.template.Template(pdb_gly)
    template.set_alpha_hull(pdb_ala, alpha=9)

    getResnums = pdb_ala.ca.getResnums() #grabs all residue numbers contained in the allALA backbone.
    if len(resnums) <= 0:
        resnums = getResnums.tolist() #converts all residue numbers to a list format that can be read
    chains = ['A'] * len(resnums)
    segs = [''] * len(resnums)
    segs_chains_resnums = zip(segs, chains, resnums)

    # for x in segs_chains_resnums:
    #     print(x)

    cgs = ['coo', 'bb_cco', 'phenol']
    cg_max_dists = dict(ph=0.9)
    outpath = input_dir

    combs2.design.functions.write_resfile(template, CGs=cgs,
                                            outpath=outpath,
                                            filename='resfile', tag='',
                                            resindices=None, segs_chains_resnums=segs_chains_resnums,
                                            pikaa_dict=None, bb_dep=1,
                                            use_enriched_vdMs=True, CA_burial_distance=None, exclude_exposed=False,
                                            exclude_intermed=False,
                                            exclude_buried=False, top_exposed=None, top_intermed=None, top_buried=None,
                                            alpha_hull_radius=10,
                                            use_propensities=True,
                                            propensity_threshold=0.9, use_abple=True, use_dssp=False,
                                            path_to_pdb_for_dssp=None,
                                            allowed_exposed='ACDEFGHIKLMNPQRSTVWY', allowed_intermed='ACDEFGHIKLMNPQRSTVWY',
                                            allowed_buried='ACDEFGHIKLMNPQRSTVWY',
                                            hb_only_residues='', all_contact_residues='')


    path_to_resfile= input_dir + 'resfile.txt'
    
    sc = combs2.design._sample.Sample(**dict(path_to_resfile=outpath + 'resfile.txt',
                                            path_to_database=path_to_database))
    sc.read_resfile()

    sc.load_vdms(template, filter_by_phi_psi=False, run_parallel=True)
    return sc


def run_combs(ss, sl, sc, input_dict):
    '''
    Search vdms which are close to the ligands cgs.
    '''
    ligands = sl.filtered_ligands

    labels_cgs = {}
    df_cgs = {}
    dist_ind_cgs = {}

    cg_ids = input_dict.keys()
    for cg_id in cg_ids:
        metalprot_scv.search_vdm(sc.cg_dict, ligands, cg_id, input_dict, labels_cgs, df_cgs, dist_ind_cgs, rmsd = 1.5)

    outdir = input_dir + 'combs_out/'
    os.makedirs(outdir, exist_ok= True)
    CgCombInfoDict = metalprot_scv.construct_vdm_write(outdir, ligands, labels_cgs, input_dict, df_cgs, dist_ind_cgs, clash_radius = 2.7)

    metalprot_scv.write_summary(outdir, CgCombInfoDict, name = '_summary.tsv')
    return


########################################################################################
#Define paramters here.

#query_dir = '/wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/data/database/pickle_noCYS_mn_fe_co_220119/'
query_dir = '/mnt/e/DesignData/ligands/all/20220116_FE_MN_CO/20220116_selfcenter/pickle_noCYS/'

workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_1st/'
workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/'

predefined_win_filters = [[15, 19, 27]]

lig_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/tts_fe.pdb'
ideal_geo_o_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/fe_geo_o.pdb'

lig_connects = ['O1','O3', 'FE1']
ro1 = ['C8', 'C7']
rest1 = ['C9', 'O1', 'O3', 'O4', 'FE1']
ro2 = ['C7', 'C6']
rest2 = ['H6', 'H7', 'C8', 'C9', 'O1', 'O3', 'O4', 'FE1']
clash_dist = 3.0

input_dict = {
    (0, 0):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O4'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CB', 'CG', 'OD1', 'OD2']
    },
    (0, 1):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O4'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CB', 'CG', 'OD2', 'OD1']
    },    
    (0, 2):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O4'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CB', 'CG', 'OE1', 'OE2']
    },
    (0, 3):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O4'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CB', 'CG', 'OE2', 'OE1']
    },
    (1, 0):{
        'cg' : 'phenol',
        'lgd_sel' : ['C2', 'C3', 'C4', 'O2'],
        'represent_name' : 'OH',
        'correspond_resname' : 'TYR',
        'correspond_names' : ['CE1', 'CZ', 'CE2', 'OH']
    },
    (1, 1):{
        'cg' : 'phenol',
        'lgd_sel' : ['C2', 'C3', 'C4', 'O2'],
        'represent_name' : 'OH',
        'correspond_resname' : 'TYR',
        'correspond_names' : ['CE2', 'CZ', 'CE1', 'OH']
    },
    (2, 0):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C7', 'C8', 'O1'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (2, 1):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C7', 'C8', 'O1'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (2, 2):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C7', 'C8', 'O1'],
        'represent_name' : 'O',
        'correspond_resname' : 'PRO',
        'correspond_names' : ['CA', 'C', 'O']
    }
}

#predefined_resnumss = [[20, 28, 31, 32, 36, 38, 40, 47, 53, 56, 57, 60, 61, 66, 68, 86, 88, 99, 101, 103, 105, 113, 116, 118, 120]]
predefined_resnumss = [[11, 12, 15, 16, 19, 24, 27, 30, 31, 35, 37, 39, 46, 52, 55, 56, 60, 65, 67, 69, 83, 85, 87, 98, 100, 102, 104, 106, 112, 115, 117, 119]]
########################################################################################

def main():
    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    if len(predefined_win_filters) <= 0:
        win_filters = [[] * len(pdb_files)]
    else:
        win_filters = predefined_win_filters

    if len(predefined_resnumss) <= 0:
        resnumss = [[] * len(pdb_files)]
    else:
        resnumss = predefined_resnumss

    ind = int(sys.argv[1]) -1
    if ind > len(pdb_files) -1:
        return
    target_path = workdir + pdb_files[ind]

    if len(win_filters) == 1:
        win_filter = win_filters[0]
    else:
        win_filter = win_filters[ind]
    outdir = workdir + 'output_selfcenter_' + pdb_files[ind].split('.')[0] + '_/'
    
    if len(resnumss) == 1:
        resnums = resnumss[0]
    else:
        resnums = resnumss[ind]
    run(target_path, win_filter, outdir, resnums, input_dict)
    return


def run(target_path, win_filter, outdir, resnums, input_dict):

    ss = run_metal(query_dir, workdir, outdir, target_path, win_filter)
    if len(list(ss.best_aa_comb_dict.keys())) <= 0:
        return

    lig = pr.parsePDB(lig_path)

    #search_ligands = run_ligand(ss, lig, ro1, ro2, rest1, rest2, ideal_geo_o_path, clash_dist)

    ### Only keep search_ligands where the filtered_ligands is not None.
    # sls = [sl for sl in search_ligands if sl.filtered_ligands]
    # print(len(sls))
    # if len(sls) <= 0:
    #    return

    # for sl in sls:
    #    sc = load_vdms(sl, resnums)
    #    run_combs(ss, sl, sc, input_dict)    

    return


if __name__=='__main__':
    main()