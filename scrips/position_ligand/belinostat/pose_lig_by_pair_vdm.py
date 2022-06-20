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
from nbformat import write
import pandas as pd
import numpy as np
import prody as pr
from sklearn.neighbors import NearestNeighbors
import itertools
from datetime import datetime
from scipy.sparse import csr_matrix, lil_matrix

from metalprot.combs import sample_dist
from metalprot.basic import constant, utils
from metalprot.combs import search_lig_indep


class Para:

    resnums = [3, 7, 10, 14, 17, 18, 21, 24, 25, 
        51, 54, 58, 61, 65, 68, 69, 72, 77, 81, 84, 88, 91, 92, 95, 99, 
        125, 128, 132, 135, 139, 142, 146]
    #resnums = [61]
    predefined_resnums = [('A', r) for r in resnums]

    use_enriched = True
    use_abple=True

    rmsd = 0.6

    vdm_cg_aa_atommap_dict_a = {
        ('ph_0'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CZ',
            'correspond_names' : ['CG', 'CD1', 'CD2', 'CZ'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_1'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CZ',
            'correspond_names' : ['CG', 'CD2', 'CD1', 'CZ'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_2'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CG',
            'correspond_names' : ['CD1', 'CG', 'CE1', 'CE2'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_3'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CG',
            'correspond_names' : ['CD1', 'CE1', 'CG', 'CE2'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_4'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CZ',
            'correspond_names' : ['CE1', 'CD1', 'CZ', 'CD2'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_5'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CZ',
            'correspond_names' : ['CE1', 'CZ', 'CD1', 'CD2'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_6'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CZ',
            'correspond_names' : ['CZ', 'CE1', 'CE2', 'CG'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_7'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CZ',
            'correspond_names' : ['CZ', 'CE2', 'CE1', 'CG'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_8'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CZ',
            'correspond_names' : ['CE2', 'CZ', 'CD2', 'CD1'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_9'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CZ',
            'correspond_names' : ['CE2', 'CD2', 'CZ', 'CD1'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_10'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CG',
            'correspond_names' : ['CD2', 'CE2', 'CG', 'CE1'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_11'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'represent_name' : 'CG',
            'correspond_names' : ['CD2', 'CG', 'CE2', 'CE1'],
            'aas' : 'F',
            'filter_hb' : False,
            'filter_cc' : True
        }
    }

    vdm_cg_aa_atommap_dict_b = {
        ('conh2_0'):{
            'cg' : 'conh2',
            'lgd_sel' : ['O3', 'C12', 'N2', 'C11'],
            'correspond_resname' : 'ASN',
            'represent_name' : 'CG',
            'correspond_names' : ['CB', 'CG', 'OD1', 'ND2'],
            'aas' : 'STYH', #'HDE' is also provide good hb but we want to get rid of confusing the metal binding.
            'filter_hb' : True,
            'filter_cc' : False
        },  
    }


def _calc_lig_sel_dist(lig, key_a, key_b, para):
    dists = []
    lig_sel_count = len(para.vdm_cg_aa_atommap_dict_a[key_a]['lgd_sel'])
    for i in range(lig_sel_count):
        a = para.vdm_cg_aa_atommap_dict_a[key_a]['lgd_sel'][i]
        b = para.vdm_cg_aa_atommap_dict_b[key_b]['lgd_sel'][i]
        dist = pr.calcDistance(lig.select('name ' + a), lig.select('name ' + b))[0]
        dists.append(dist)
    lig_sel_coords =[]
    for i in range(lig_sel_count):
        a = para.vdm_cg_aa_atommap_dict_a[key_a]['lgd_sel'][i]
        lig_sel_coords.append(lig.select('name ' + a)[0].getCoords())
    for i in range(lig_sel_count):
        b = para.vdm_cg_aa_atommap_dict_b[key_b]['lgd_sel'][i]
        lig_sel_coords.append(lig.select('name ' + b)[0].getCoords())
    return dists, np.array(lig_sel_coords)

def _vdm_bb_clashing(target, chidres, df_vdm_filter, dist_cut = 2.5):
    '''
    
    '''
    df_vdm_filter_heavy = df_vdm_filter[(df_vdm_filter['resnum'] == 10)
            & ~(df_vdm_filter['name'].isin(['N', 'C', 'CA', 'O']))
            & ~(df_vdm_filter['atom_type_label'].isin(['h_pol', 'h_alkyl', 'h_aro']))]

    bb_coords = target.select('heavy and not ( chid ' + chidres[0] + ' and resnum ' + str(chidres[1]) + ')').getCoords()

    nn = NearestNeighbors(radius= dist_cut).fit(bb_coords)
    x_in_y = nn.radius_neighbors_graph(df_vdm_filter_heavy[['c_x', 'c_y', 'c_z']]).toarray().astype(bool)
    clashing = np.any(x_in_y, axis=1)

    labels_clashing = df_vdm_filter_heavy[['CG', 'rota', 'probe_name']][clashing]

    return labels_clashing

def _get_labels_and_vdm_coords(target, chidres, vdm_cg_aa_atommap_dict, cg_key, abple, aa, para, path_to_database):

    df_vdm = search_lig_indep.load_old_vdm(path_to_database, vdm_cg_aa_atommap_dict[cg_key]['cg'], aa)

    df_vdm_filter = search_lig_indep.filter_db(df_vdm, use_enriched = para.use_enriched, use_abple=para.use_abple, abple=abple)        
    labels = df_vdm_filter[['CG', 'rota', 'probe_name']]

    #Get vdm sc clashing filtered labels index.
    labels_clashing = _vdm_bb_clashing(target, chidres, df_vdm_filter, dist_cut = 2.5)
    labels_clashing = labels_clashing.drop_duplicates().astype(str).T.agg(','.join)

    labels_clashfilter_ind = ~labels.astype(str).T.agg(','.join).isin(labels_clashing)
    df_vdm_filter_filter = df_vdm_filter[labels_clashfilter_ind]

    pos = target.select('chid ' + chidres[0] + ' and resnum ' + str(chidres[1]))
    tf = pr.calcTransformation(constant.ideal_ala_coords, pos.select('name N CA C').getCoords())

    df_vdm_filter_filter.loc[:, ['c_x', 'c_y', 'c_z']] = tf.apply(df_vdm_filter_filter[['c_x', 'c_y', 'c_z']].to_numpy())

    labels, vdm_coords = search_lig_indep.get_vdm_labels_coords_4old_vdm_db(df_vdm_filter_filter, vdm_cg_aa_atommap_dict[cg_key]['correspond_resname'], vdm_cg_aa_atommap_dict[cg_key]['represent_name'], vdm_cg_aa_atommap_dict[cg_key]['correspond_names'])
    
    return labels, vdm_coords, df_vdm_filter_filter


def _get_pair_vdms_adj_matrix(dists, labels_a, vdm_coords_a, labels_b, vdm_coords_b, para):
    '''
    
    '''
    adj_matrix_pair = lil_matrix((len(labels_a), len(labels_b)), dtype=bool)

    for i in range(len(dists)):
        dist = dists[i]
        coords_a = vdm_coords_a[:, i*3:i*3+3]
        coords_b = vdm_coords_b[:, i*3:i*3+3]
        nn_low = NearestNeighbors(radius= dist - para.rmsd).fit(coords_b)
        a_in_b_low = nn_low.radius_neighbors_graph(coords_a).toarray().astype(bool)

        nn_high= NearestNeighbors(radius= dist + para.rmsd).fit(coords_b)
        a_in_b_high = nn_high.radius_neighbors_graph(coords_a).toarray().astype(bool)
        if i == 0:
            adj_matrix_pair = lil_matrix(a_in_b_high > a_in_b_low, dtype=bool)
            if adj_matrix_pair.sum() == 0:
                return adj_matrix_pair
        else:
            adj_matrix_pair = adj_matrix_pair.multiply(a_in_b_high > a_in_b_low)
            if adj_matrix_pair.sum() == 0:
                return adj_matrix_pair
    return adj_matrix_pair.tolil()


def _lig_on_pair_vdms(lig, para, lig_sel_coords, adj_matrix_pair, labels_a, vdm_coords_a, labels_b, vdm_coords_b):
    lig_la_lbs = []
    for r in range(adj_matrix_pair.shape[0]):
        for c in adj_matrix_pair.rows[r]:
            coords = []
            coords.extend(vdm_coords_a[r].reshape(-1, 3))
            coords.extend(vdm_coords_b[c].reshape(-1, 3))
            _rmsd = pr.calcRMSD(lig_sel_coords, np.array(coords))
            if _rmsd > para.rmsd:
                continue
            _lig = lig.copy()
            tf = pr.calcTransformation(lig_sel_coords, np.array(coords))
            tf.apply(_lig)
            lig_la_lbs.append((labels_a.iloc[r], labels_b.iloc[c], _lig))
    return lig_la_lbs


def _search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b):
    #>>> For the selected key, get the dists of a pair of ligand cgs.
    dists, lig_sel_coords = _calc_lig_sel_dist(lig, key_a, key_b, para)

    for _aa_a in para.vdm_cg_aa_atommap_dict_a[key_a]['aas']:
        aa_a = constant.inv_one_letter_code[_aa_a]
        labels_a, vdm_coords_a, df_vdm_filter_filter_a = _get_labels_and_vdm_coords(target, chidres_a, para.vdm_cg_aa_atommap_dict_a, key_a, abple_a, aa_a, para, path_to_database)

        for _aa_b in para.vdm_cg_aa_atommap_dict_b[key_b]['aas']:
            aa_b = constant.inv_one_letter_code[_aa_b]
            labels_b, vdm_coords_b, df_vdm_filter_filter_b = _get_labels_and_vdm_coords(target, chidres_b, para.vdm_cg_aa_atommap_dict_b, key_b, abple_b, aa_b, para, path_to_database)

            adj_matrix_pair = _get_pair_vdms_adj_matrix(dists, labels_a, vdm_coords_a, labels_b, vdm_coords_b, para)
            
            # with open(outdir + '_summary.txt', 'a') as f:
            #     f.write('----------------------------------------------------------\n')
            #     f.write('pos: ({}, {}) key_a: {}, key_b: {}, aa_a:{}, aa_b {}, matrix sum: {}\n'.format(chidres_a, chidres_b, key_a, key_b, aa_a, aa_b, adj_matrix_pair.sum()))

            if adj_matrix_pair.sum() == 0:
                continue
            lig_la_lbs = _lig_on_pair_vdms(lig, para, lig_sel_coords, adj_matrix_pair, labels_a, vdm_coords_a, labels_b, vdm_coords_b)
            print('pos: ({}, {}) key_a: {}, key_b: {}, aa_a:{}, aa_b {}, matrix sum: {}, lig_la_lbs: {}\n'.format(chidres_a, chidres_b, key_a, key_b, aa_a, aa_b, adj_matrix_pair.sum(), len(lig_la_lbs)))
            
            if len(lig_la_lbs)>0:
                title = '_'.join([key_a, key_b, chidres_a[0], str(chidres_a[1]), chidres_b[0], str(chidres_b[0])]) + '_summary.txt'
                with open(outdir + title, 'a') as f:
                    for la, lb, lig in lig_la_lbs:
                        la_str = '_'.join([str(x) for x in la])
                        lb_str = '_'.join([str(x) for x in lb])
                        f.write(target.getTitle() + '\t' + lig.getTitle() + '\t' + la_str + '\t' + lb_str + '\n')
                        pr.writePDB(outdir + 'Lig_' + la_str + lb_str, lig)



def _select_chidres_keys(target, lig, para, path_to_database):
    abples, phipsi = utils.seq_get_ABPLE(target)

    select_chidres_keys = []
    for i, j in itertools.combinations(range(len(para.predefined_resnums)), 2):
        chidres_a = para.predefined_resnums[i]
        chidres_b = para.predefined_resnums[j]

        pos_a = target.select('chid ' + chidres_a[0] + ' and resnum ' + str(chidres_a[1]))
        pos_b = target.select('chid ' + chidres_b[0] + ' and resnum ' + str(chidres_b[1]))

        far = pr.calcDistance(pos_a.select('name C'), pos_b.select('name C'))[0]
        if far > 18:
            continue
        
        abple_a = abples[pos_a.getResindices()[0]]
        abple_b = abples[pos_b.getResindices()[0]]

        for key_a in para.vdm_cg_aa_atommap_dict_a.keys():
            for key_b in para.vdm_cg_aa_atommap_dict_b.keys():
                select_chidres_keys.append((key_a, key_b, chidres_a, chidres_b, abple_a, abple_b))
        
    return select_chidres_keys

def run_local():
    workdir = '/mnt/e/DesignData/Metalloenzyme/'

    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

    lig = pr.parsePDB(workdir + 'ligs/meo_50g_amber14eht_md_out/50g_md_0.pdb')

    target = pr.parsePDB(workdir + 'targets/01_f63440_nick_ala.pdb')

    para = Para()

    select_chidres_keys = _select_chidres_keys(target, lig, para, path_to_database)
    for key_a, key_b, chidres_a, chidres_b, abple_a, abple_b in select_chidres_keys:
        _search_select_pair_vdm(workdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)
    return

def run_wynton():
    workdir = '/wynton/home/degradolab/lonelu/DesignData/Metalloenzyme/'

    path_to_database='/wynton/home/degradolab/lonelu/DesignData/Database/vdMs/'

    lig = pr.parsePDB(workdir + 'meo_50g_amber14eht_md_out/50g_md_0.pdb')

    target = pr.parsePDB(workdir + 'targets/01_f63440_nick_ala.pdb')

    outdir = workdir + 'results/' 
    os.makedirs(outdir)

    para = Para()

    select_chidres_keys = _select_chidres_keys(target, lig, para, path_to_database)

    ind = int(sys.argv[1]) -1

    key_a, key_b, chidres_a, chidres_b, abple_a, abple_b = select_chidres_keys[ind]
    _search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)
    return

def run_wynton_multifile():
    workdir = '/wynton/home/degradolab/lonelu/DesignData/Metalloenzyme/'

    path_to_database='/wynton/home/degradolab/lonelu/DesignData/Database/vdMs/'

    ligs = [pr.parsePDB(workdir + 'meo_50g_amber14eht_md_out/' + x) for x in os.listdir(workdir + 'meo_50g_amber14eht_md_out/') if '.pdb' in x]

    targets = [pr.parsePDB(workdir + 'targets/' + x) for x in os.listdir(workdir + 'targets/') if '.pdb' in x]
    
    para = Para()

    ind = int(sys.argv[1]) -1
    select_chidres_keys = _select_chidres_keys(targets[0], ligs[0], para, path_to_database)
    key_a, key_b, chidres_a, chidres_b, abple_a, abple_b = select_chidres_keys[ind]
    for target in targets:
        for lig in ligs:           
            outdir = workdir + 'results_' + target.getTitle() + '_' + lig.getTitle() + '/'
            os.makedirs(outdir, exist_ok = True)
            _search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)

    return

if __name__=='__main__':
    run_wynton_multifile()

def _test_case():
    '''
    The test case include a real example to develop the functions.
    '''
    workdir = '/mnt/e/DesignData/Metalloenzyme/'

    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

    lig = pr.parsePDB(workdir + 'ligs/meo_50g_amber14eht_md_out/50g_md_0.pdb')

    target = pr.parsePDB(workdir + 'targets/01_f63440_nick_ala.pdb')

    para = Para()

    abples, phipsi = utils.seq_get_ABPLE(target)

    i = 1
    j = 3
    chidres_a = para.predefined_resnums[i]
    chidres_b = para.predefined_resnums[j]

    pos_a = target.select('chid ' + chidres_a[0] + ' and resnum ' + str(chidres_a[1]))
    pos_b = target.select('chid ' + chidres_b[0] + ' and resnum ' + str(chidres_b[1]))

    abple_a = abples[pos_a.getResindices()[0]]
    abple_b = abples[pos_b.getResindices()[0]]


    key_a = 'ph_0'
    key_b = 'conh2_0'
    dists, lig_sel_coords = _calc_lig_sel_dist(lig, key_a, key_b, para)


    aa_a = 'PHE'
    aa_b = 'SER'
    labels_a, vdm_coords_a = _get_labels_and_vdm_coords(target, chidres_a, para.vdm_cg_aa_atommap_dict_a, key_a, abple_a, aa_a, pos_a, para)
    labels_b, vdm_coords_b = _get_labels_and_vdm_coords(target, chidres_b, para.vdm_cg_aa_atommap_dict_b, key_b, abple_b, aa_b, pos_b, para)

    adj_matrix_pair = _get_pair_vdms_adj_matrix(dists, labels_a, vdm_coords_a, labels_b, vdm_coords_b, para)

    print(adj_matrix_pair.sum() == 0)

    lig_la_lbs = _lig_on_pair_vdms(lig, lig_sel_coords, adj_matrix_pair, labels_a, vdm_coords_a, labels_b, vdm_coords_b)

    if len(lig_la_lbs) > 0:
        print(str(i) + ' ' + str(j))
        print(aa_a + ' ' + aa_b)
    return


