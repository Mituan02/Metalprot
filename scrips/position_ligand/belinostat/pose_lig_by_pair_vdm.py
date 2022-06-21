'''
The idea is to pose lig based on a pair of vdms each with different cgs of the ligs.
The method could really solve the lever arm effect.
The method could also in theory remove a lot of single contacts.
Using belinostat as an example, I can use two cgs phenol (C7 C8 C9 C13 C14 C5) and conh2 (C11 C12 O3 N2). 
After put vdms on bb, I can check the distances between the two cgs and find the pairs that can match the lig. 
Please check metalprot.combs.sample_dist.py
'''

import sys
import os
import pandas as pd
import numpy as np
import prody as pr
from sklearn.neighbors import NearestNeighbors
import itertools
from datetime import datetime
from scipy.sparse import csr_matrix, lil_matrix

from metalprot.basic import constant, utils
from metalprot.combs import search_lig_indep, search_lig_indep_inpair

'''
python /mnt/e/Github_Design/Metalprot/scrips/position_ligand/belinostat/pose_lig_by_pair_vdm.py

'''
class Para:

    resnums = [3, 7, 10, 14, 17, 18, 21, 24, 25, 
        51, 54, 58, 61, 65, 68, 69, 72, 77, 81, 84, 88, 91, 92, 95, 99, 
        125, 128, 132, 135, 139, 142, 146]
    #resnums = [61]
    predefined_resnums = [('A', r) for r in resnums]

    use_enriched = True
    use_abple=True

    rmsd = 0.6

    # vdm_cg_aa_atommap_dict_a = {
    #     ('ph_0'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CG', 'CD1', 'CD2', 'CZ'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_1'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CG', 'CD2', 'CD1', 'CZ'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_2'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CG',
    #         'correspond_names' : ['CD1', 'CG', 'CE1', 'CE2'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_3'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CG',
    #         'correspond_names' : ['CD1', 'CE1', 'CG', 'CE2'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_4'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CE1', 'CD1', 'CZ', 'CD2'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_5'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CE1', 'CZ', 'CD1', 'CD2'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_6'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CZ', 'CE1', 'CE2', 'CG'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_7'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CZ', 'CE2', 'CE1', 'CG'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_8'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CE2', 'CZ', 'CD2', 'CD1'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_9'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CE2', 'CD2', 'CZ', 'CD1'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_10'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CG',
    #         'correspond_names' : ['CD2', 'CE2', 'CG', 'CE1'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_11'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CG',
    #         'correspond_names' : ['CD2', 'CG', 'CE2', 'CE1'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     }
    # }

    # vdm_cg_aa_atommap_dict_a = {
    #     ('ph_0'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C9', 'C14'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CG', 'CE1', 'CE2'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    # }

    vdm_cg_aa_atommap_dict_a ={
        ('coo_0'):{
            'cg' : 'coo',
            'lgd_sel' : ['CG', 'OD1', 'OD2'],
            'correspond_resname' : 'ASP',
            'represent_name' : 'CG',
            'correspond_names' : ['S', 'O1', 'O2'],
            'aas' : 'STYH', #'HDE' is also provide good hb but we want to get rid of confusing the metal binding.
            'filter_hb' : True,
            'filter_cc' : False
        },         
    }


    vdm_cg_aa_atommap_dict_b = {
        ('conh2_0'):{
            'cg' : 'conh2',
            'lgd_sel' : ['O3', 'C12', 'N2'],
            'correspond_resname' : 'ASN',
            'represent_name' : 'CG',
            'correspond_names' : ['OD1', 'CG', 'ND2'],
            'aas' : 'STYH', #'HDE' is also provide good hb but we want to get rid of confusing the metal binding.
            'filter_hb' : True,
            'filter_cc' : False
        },  
    }


def run_local():
    workdir = '/mnt/e/DesignData/Metalloenzyme/belinostat/'

    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

    lig = pr.parsePDB(workdir + 'ligs/meo_50g_amber14eht_md_out/50g_md_0.pdb')

    target = pr.parsePDB(workdir + 'targets/01_f63440_nick_ala.pdb')

    outdir = workdir + 'results/' 
    os.makedirs(outdir)

    para = Para()

    select_chidres_keys = search_lig_indep_inpair._select_chidres_keys(target, lig, para, path_to_database)
    
    for key_a, key_b, chidres_a, chidres_b, abple_a, abple_b in select_chidres_keys:
        search_lig_indep_inpair.search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)
    return

def run_wynton():
    workdir = '/wynton/home/degradolab/lonelu/DesignData/Metalloenzyme/'

    path_to_database='/wynton/home/degradolab/lonelu/DesignData/Database/vdMs/'

    lig = pr.parsePDB(workdir + 'meo_50g_amber14eht_md_out/50g_md_0.pdb')

    target = pr.parsePDB(workdir + 'targets/01_f63440_nick_ala.pdb')

    outdir = workdir + 'results/' 
    os.makedirs(outdir)

    para = Para()

    select_chidres_keys = search_lig_indep_inpair.select_chidres_keys(target, lig, para, path_to_database)

    ind = int(sys.argv[1]) -1

    key_a, key_b, chidres_a, chidres_b, abple_a, abple_b = select_chidres_keys[ind]
    search_lig_indep_inpair._search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)
    return

def run_wynton_multifile():
    workdir = '/wynton/home/degradolab/lonelu/DesignData/Metalloenzyme/'

    path_to_database='/wynton/home/degradolab/lonelu/DesignData/Database/vdMs/'

    ligs = [pr.parsePDB(workdir + 'meo_50g_amber14eht_md_out/' + x) for x in os.listdir(workdir + 'meo_50g_amber14eht_md_out/') if '.pdb' in x]

    targets = [pr.parsePDB(workdir + 'targets/' + x) for x in os.listdir(workdir + 'targets/') if '.pdb' in x]
    
    para = Para()

    ind = int(sys.argv[1]) -1
    select_chidres_keys = search_lig_indep_inpair._select_chidres_keys(targets[0], ligs[0], para, path_to_database)
    key_a, key_b, chidres_a, chidres_b, abple_a, abple_b = select_chidres_keys[ind]
    for target in targets:
        for lig in ligs:           
            outdir = workdir + 'results_' + target.getTitle() + '_' + lig.getTitle() + '/'
            os.makedirs(outdir, exist_ok = True)
            search_lig_indep_inpair.search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)

    return

if __name__=='__main__':
    #run_wynton_multifile()
    run_local()



