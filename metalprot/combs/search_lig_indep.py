'''
The basic functions to search/eval ligs or 2ndshell.
'''

import itertools
import os
import pdb
import prody as pr
import pickle
import pandas as pd
import numpy as np
from metalprot.basic import constant, utils, prody_ext
from metalprot.combs import gvdm_helper, position_ligand
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import lil_matrix
import shutil


def prepare_rosetta_target(outdir, target_path, designed_site = [('A', 15), ('A', 19), ('A', 27)]):
    '''
    Or searching rosetta 1st round designed target backbone. 
    All aas are mutated into gly except the metal binding ones.
    '''
    target = pr.parsePDB(target_path)

    resindices, chidres2ind, ind2chidres = prody_ext.transfer2resindices(target, designed_site)

    pdb_gly = prody_ext.target_to_all_gly_ala(target, target.getTitle() + '_gly.pdb', resindices, aa= 'GLY', keep_no_protein = False)
    
    pr.writePDB(outdir + pdb_gly.getTitle() + '.pdb', pdb_gly)

    return pdb_gly, chidres2ind


def lgd_sel_coord(lgd, lgd_sel):
    '''
    The lgd_sel is a list of selected atom.
    '''
    coords = []
    for lsa in lgd_sel:
        coords.extend(lgd.select('name ' + lsa).getCoords().flatten())
    return coords


def get_ligand_coords(filtered_ligands, lgd_sel):
    # Get all ligand coords.
    ligand_coords = []
    for i in range(len(filtered_ligands)):
        lgd = filtered_ligands[i]
        coords = lgd_sel_coord(lgd, lgd_sel)
        ligand_coords.append(coords)
    ligand_coords = np.array(ligand_coords)

    return ligand_coords
    

def get_nearest_vdms_rmsd(vdm_coords, ligand_coords, radius = 1):
    # Nearest Neighbor
    nbr = NearestNeighbors(radius=radius).fit(vdm_coords)
    dists, inds = nbr.radius_neighbors(ligand_coords)

    return dists, inds


def clash_filter_protein(target, all_coords, dist_cut = 2.5):
    '''
    # clashing: any atom close to the bb let's say dist <= dist_cut

    The bb here should contain 'N CA C O of any aa, the sc of binding aa, Metal'
    '''
    
    bb_coords = target.select('heavy').getCoords()
    nn = NearestNeighbors(radius= dist_cut).fit(bb_coords)
    x_in_y = nn.radius_neighbors_graph(all_coords).toarray().astype(bool)
    clashing = np.any(x_in_y, axis=1)

    #np.sum(clashing) #The clashing is too less to be considered.
    return clashing


def clash_filter_proteins(targets, dist_cut = 2.2):
    '''
    # clashing: any atom close to the bb of a pair pdbs let's say dist <= dist_cut
    return True if clashing. 
    '''
    for i, j in itertools.permutations(range(len(targets))):
        t1 = targets[i]
        t2 = targets[j]
        nbrs = NearestNeighbors(radius= dist_cut).fit(t1.select('heavy').getCoords())
        adj_matrix = nbrs.radius_neighbors_graph(t2.select('heavy').getCoords()).astype(bool)

        if np.sum(adj_matrix) >0:
            return True
    #np.sum(clashing) #The clashing is too less to be considered.
    return False


def clash_filter_protein_single(target, chidres, vdm_ag, dist_cut = 2.7):
    '''
    # clashing: any atom close to the bb let's say dist <= dist_cut

    The bb here should contain 'N CA C O of any aa, the sc of binding aa, Metal'
    '''
    target_coords = target.select('heavy and not ( chid ' + chidres[0] + ' and resnum ' + str(chidres[1]) + ')').getCoords()
    vdm_sel = vdm_ag.select('heavy and sc and chain X and resnum 10') #For Gly, it can be None.
    if not vdm_sel:
        return True
    vdm_coords = vdm_sel.getCoords()
    nbrs = NearestNeighbors(radius= dist_cut).fit(target_coords)
    adj_matrix = nbrs.radius_neighbors_graph(vdm_coords).astype(bool)

    if np.sum(adj_matrix) >0:
        return True
    return False


def clash_filter_lig(lig, vdm, dist_cut = 2.5):
    '''
    Check if the vdm clashing the ligand
    '''
    lig_coords = lig.select('heavy').getCoords()
    vdm_sel = vdm.select('heavy and sc and chain X and resnum 10') #For Gly, it can be None.
    if not vdm_sel:
        return True
    vdm_coords = vdm_sel.getCoords()
    nbrs = NearestNeighbors(radius= dist_cut).fit(lig_coords)
    adj_matrix = nbrs.radius_neighbors_graph(vdm_coords).astype(bool)

    if np.sum(adj_matrix) >0:
        return True
    return False


def search_lig_at_cg_aa_resnum(target, chidres, pos, abple, ligs, vdm_cg_aa_atommap_dict, cg_id, df_vdm, ideal_ala_coords, rmsd = 0.7, filter_hb = False, filter_cc = False):
    
    results = []

    _ligs, tf_rev, tf = gvdm_helper.ligscoords_2_ideal_ala(pos, ligs, ideal_ala_coords)

    #<<< Debug
    # outdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/search_2ndshell_result/'
    # _target = target.copy()
    # tf.apply(_target)
    # pr.writePDB(outdir + '_target.pdb', _target)
    # pr.writePDB(outdir + '_lig.pdb', _ligs[0])
    #>>>

    #print(input_dict[cg_id]['lgd_sel'])
    #print(_ligs)
    ligand_coords = get_ligand_coords(_ligs, vdm_cg_aa_atommap_dict[cg_id]['lgd_sel'])

    labels, vdm_coords = gvdm_helper.get_vdm_labels_coords_4old_vdm_db(df_vdm, vdm_cg_aa_atommap_dict[cg_id]['correspond_resname'], vdm_cg_aa_atommap_dict[cg_id]['represent_name'], vdm_cg_aa_atommap_dict[cg_id]['correspond_names'])

    #<<< Debgu
    # for i in range(labels.shape[0]):
    #     x = labels.iloc[i]
    #     v = df_vdm[(df_vdm['CG'] == x['CG']) & (df_vdm['rota'] == x['rota']) & (df_vdm['probe_name'] == x['probe_name'])]
    #     title = '-'.join([str(_s) for _s in (x['CG'], x['rota'], x['probe_name'])])
    #     ag = gvdm_helper.df2ag(v, title)
    #     pr.writePDB(outdir + title + '.pdb', ag)
    #>>>

    if vdm_coords.shape[0] <=0 or vdm_coords.shape[0] != labels.shape[0]:
        print('cg_id not working: {}'.format(cg_id))
        return results
    
    num_cg_atoms = len(vdm_cg_aa_atommap_dict[cg_id]['lgd_sel'])
    radius = np.sqrt(num_cg_atoms) * rmsd
    dists, inds = get_nearest_vdms_rmsd(vdm_coords, ligand_coords, radius = radius)

    #u is lig id, v is vdM id.
    for u in range(len(inds)): 
        for z in range(len(inds[u])):                
            ind = inds[u][z]
            rmsd = dists[u][z]/np.sqrt(num_cg_atoms)
            try:
                x = labels.iloc[ind]
                v = df_vdm[(df_vdm['CG'] == x['CG']) & (df_vdm['rota'] == x['rota']) & (df_vdm['probe_name'] == x['probe_name'])]

                if (not filter_hb and not filter_cc) or (filter_hb and v[['contact_hb']].any()[0]) or (filter_cc and v[['contact_cc']].any()[0]):
                    title = '-'.join([str(_s) for _s in (x['CG'], x['rota'], x['probe_name'])])
                    ag = gvdm_helper.df2ag(v, title)
                    tf_rev.apply(ag)

                    if clash_filter_protein_single(target, chidres, ag) or clash_filter_lig(ligs[u], ag):
                        print('filtered clash_filter_protein_single, clash_filter_lig')
                        continue
                    results.append((cg_id, u, rmsd, ag, ind, v['C_score_ABPLE_' + abple].values[0], None, v[['contact_hb']].any()[0], v[['contact_cc']].any()[0]))
            except:
                print(cg_id)
                print(x)

    return results




