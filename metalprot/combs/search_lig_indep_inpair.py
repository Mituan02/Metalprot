
import sys
import os
import pandas as pd
import numpy as np
import prody as pr
from sklearn.neighbors import NearestNeighbors
import itertools
from datetime import datetime
from scipy.sparse import csr_matrix, lil_matrix

from ..basic import constant, utils
from ..combs import gvdm_helper, search_lig_indep


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

    df_vdm = gvdm_helper.load_old_vdm(path_to_database, vdm_cg_aa_atommap_dict[cg_key]['cg'], aa)

    df_vdm_filter = gvdm_helper.filter_db(df_vdm, use_enriched = para.use_enriched, use_abple=para.use_abple, abple=abple)        
    labels = df_vdm_filter[['CG', 'rota', 'probe_name']]

    #Get vdm sc clashing filtered labels index.
    labels_clashing = _vdm_bb_clashing(target, chidres, df_vdm_filter, dist_cut = 2.5)
    labels_clashing = labels_clashing.drop_duplicates().astype(str).T.agg(','.join)

    labels_clashfilter_ind = ~labels.astype(str).T.agg(','.join).isin(labels_clashing)
    df_vdm_filter_filter = df_vdm_filter[labels_clashfilter_ind]

    pos = target.select('chid ' + chidres[0] + ' and resnum ' + str(chidres[1]))
    tf = pr.calcTransformation(constant.ideal_ala_coords, pos.select('name N CA C').getCoords())

    df_vdm_filter_filter.loc[:, ['c_x', 'c_y', 'c_z']] = tf.apply(df_vdm_filter_filter[['c_x', 'c_y', 'c_z']].to_numpy())

    labels, vdm_coords = gvdm_helper.get_vdm_labels_coords_4old_vdm_db(df_vdm_filter_filter, vdm_cg_aa_atommap_dict[cg_key]['correspond_resname'], vdm_cg_aa_atommap_dict[cg_key]['represent_name'], vdm_cg_aa_atommap_dict[cg_key]['correspond_names'])
    
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


def search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b):
    #>>> For the selected key, get the dists of a pair of ligand cgs.
    dists, lig_sel_coords = _calc_lig_sel_dist(lig, key_a, key_b, para)
    #print(dists)
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


def _test_case(para):
    '''
    The test case include a real example to develop the functions.
    '''
    workdir = '/mnt/e/DesignData/Metalloenzyme/'

    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

    lig = pr.parsePDB(workdir + 'ligs/meo_50g_amber14eht_md_out/50g_md_0.pdb')

    target = pr.parsePDB(workdir + 'targets/01_f63440_nick_ala.pdb')

    #para = Para()

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