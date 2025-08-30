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
import datetime
from scipy.sparse import csr_matrix, lil_matrix

from metalprot.basic import constant, utils
from metalprot.combs import search_lig_indep, search_lig_indep_inpair

'''
python /mnt/e/Github_Design/Metalprot/scrips/position_ligand/saha/run_pose_lig_by_pair_vdm.py local

'''
class Para:

    resnums = [3, 7, 10, 14, 17, 18, 21, 24, 25, 51, 54, 58, 61, 65, 68, 69, 72, 77, 81, 84, 88, 91, 92, 95, 99, 125, 128, 132, 135, 139, 142, 146]

    predefined_resnums_a = [('A', r) for r in resnums]

    predefined_resnums_b = [('A', r) for r in resnums]

    use_enriched = True
    use_abple=True

    rmsd = 0.75

    vdm_cg_aa_atommap_dict_a = {
        ('bb_cnh_0'):{
            'cg' : 'bb_cnh',
            'lgd_sel' : ['C9', 'N1', 'H13'],
            'correspond_resname' : 'GLY',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'N', 'H'],
            'aas' : 'STYHQNDEAGP', 
            'filter_hb' : True,
            'filter_cc' : False
        },   
        ('bb_cnh_1'):{
            'cg' : 'bb_cnh',
            'lgd_sel' : ['C9', 'N1', 'H13'],
            'correspond_resname' : 'ALA',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'N', 'H'],
            'aas' : 'STYHQNDEAGP', 
            'filter_hb' : True,
            'filter_cc' : False
        },   
        ('bb_cnh_2'):{
            'cg' : 'bb_cnh',
            'lgd_sel' : ['C9', 'N1', 'H13'],
            'correspond_resname' : 'PRO',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'N', 'H'],
            'aas' : 'STYHQNDEAGP', 
            'filter_hb' : True,
            'filter_cc' : False
        },
    }

    vdm_cg_aa_atommap_dict_b = {
        ('conh2_0'):{
            'cg' : 'conh2',
            'lgd_sel' : ['O1', 'C7', 'N1'],
            'correspond_resname' : 'ASN',
            'represent_name' : 'CG',
            'correspond_names' : ['OD1', 'CG', 'ND2'],
            'aas' : 'STYWHQNDE',
            #'aas' : 'E',
            'filter_hb' : True,
            'filter_cc' : False
        },  
        ('conh2_1'):{
            'cg' : 'conh2',
            'lgd_sel' : ['O1', 'C7', 'N1'],
            'correspond_resname' : 'GLN',
            'represent_name' : 'CG',
            'correspond_names' : ['OE1', 'CD', 'NE2'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },  
        ('bb_cco_0'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C5', 'C7', 'O1'],
            'correspond_resname' : 'GLY',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_1'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C5', 'C7', 'O1'],
            'correspond_resname' : 'ALA',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_2'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C5', 'C7', 'O1'],
            'correspond_resname' : 'PRO',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        }
    }


def run_local():
    #>>>
    workdir = '/mnt/e/DesignData/Metalloenzyme/SAHA_Vorinostat/'

    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

    lig = pr.parsePDB(workdir + 'SAHA_5311.pdb')

    target = pr.parsePDB(workdir + 'targets/01_f63440_nick_ala.pdb')
    
    outdir = workdir + 'out_pair/'

    os.makedirs(outdir, exist_ok = True)

    para = Para()

    select_chidres_keys = search_lig_indep_inpair._select_chidres_keys_inpair(target, para)
    print(len(select_chidres_keys))
    #for _x in select_chidres_keys:
    #    print(_x)

    for key_a, key_b, chidres_a, chidres_b, abple_a, abple_b in select_chidres_keys[0:10]:
        search_lig_indep_inpair.search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)
    return

def run_wynton():
    #>>>

    workdir = '/wynton/home/degradolab/lonelu/DesignData/Metalloenzyme/SAHA_Vorinostat/'

    path_to_database='/wynton/home/degradolab/lonelu/DesignData/Database/vdMs/'

    lig = pr.parsePDB(workdir + 'SAHA_5311.pdb')

    target = pr.parsePDB(workdir + 'targets/01_f63440_nick_ala.pdb')

    para = Para()

    select_chidres_keys = search_lig_indep_inpair._select_chidres_keys_inpair(target, para)


    #outdir = workdir + 'out_pair_' + datetime.datetime.now().strftime('%Y%m%d-%H%M%S') + '/'
    outdir = workdir + 'out_pair/'
    os.makedirs(outdir, exist_ok = True)

    ind = int(sys.argv[2]) -1
    _range = int(sys.argv[3])

    for _ind in range(ind*_range, ind*_range + _range):
        key_a, key_b, chidres_a, chidres_b, abple_a, abple_b = select_chidres_keys[_ind]
        search_lig_indep_inpair.search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)
    return

def run_wynton_multifile():
    workdir = '/wynton/home/degradolab/lonelu/DesignData/Metalloenzyme/belinostat/'

    path_to_database='/wynton/home/degradolab/lonelu/DesignData/Database/vdMs/'

    ligs = [pr.parsePDB(workdir + 'meo_50g_amber14eht_md_out/' + x) for x in os.listdir(workdir + 'meo_50g_amber14eht_md_out/') if '.pdb' in x]

    targets = [pr.parsePDB(workdir + 'targets/' + x) for x in os.listdir(workdir + 'targets/') if '.pdb' in x]
    
    para = Para()

    ind = int(sys.argv[2]) -1
    select_chidres_keys = search_lig_indep_inpair._select_chidres_keys_inpair(targets[0], para)
    key_a, key_b, chidres_a, chidres_b, abple_a, abple_b = select_chidres_keys[ind]
    for target in targets:
        for lig in ligs:           
            outdir = workdir + 'results_' + target.getTitle() + '_' + lig.getTitle() + '/'
            os.makedirs(outdir, exist_ok = True)
            search_lig_indep_inpair.search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)

    return


if __name__=='__main__':
    if sys.argv[1] == 'wynton_multi': 
        run_wynton_multifile()
    elif sys.argv[1] == 'wynton_single': 
        run_wynton()
    elif sys.argv[1] == 'local':
        run_local()



