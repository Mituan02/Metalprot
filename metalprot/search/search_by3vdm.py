'''
The AAdAAdAA database is used for the search. 
The AAdAAdAA database contains binding cores with combined 3 single vdMs.
This is the searching part for the ShaB deep learning metal position prediction.

The searching is taking the advantage of the Cluster algorithm.
First the searching coords is appened to the database coords (db_coords).
Then cluster.make_pairwise_rmsd_mat() is applied. 
Then the pairwise_rmsd_mat is used to extract db_pdbs.

'''

import os
import numpy as np
import prody
import pickle
from ..basic import cluster


def load_db3vdm(db_path):
    '''
    Load the AAdAAdAA database.
    '''
    db_pdbs, db_coords = None
    return db_pdbs, db_coords


def get_coords(target, win_comb, metal_coord):
    '''
    Given win_comb, extract the bb coords.
    '''
    coords = [target.select('chid ' + w[0] + ' and resnum ' + str(w[2]) + 'name N C CA O').getCoords() for w in win_comb]
    coords.append(metal_coord)
    return  np.array(coords, dtype ='float32').reshape((1,-1))


def cluster_search(db_coord_all, inds, rmsd):
    '''
    Cluster based searching method.
    The target bb searching coords are appended into the database coords.
    The cluster make_adj_mat is calculated.
    '''
    #>>> Prepare cluster
    clu = cluster.Cluster()
    clu.rmsd_cutoff = rmsd
    clu.pdb_coords = db_coord_all

    #>>> search rmsd mat
    clu.make_pairwise_rmsd_mat()  
    clu.make_adj_mat()

    indices = np.arange(clu.adj_mat.shape[0])

    all_mems = {}
    for cent in inds:
        row = clu.adj_mat.getrow(cent)
        tf = row.toarray().astype(bool)[0]
        mems = indices[tf]
        all_mems[cent] = mems
    return all_mems


def run_search_by3vdm(ss, db_path, win_combs, metal_coords):
    '''
    ss: Search_vdM
    '''
    #>>>Load db
    db_pdbs, db_coords = load_db3vdm(db_path)

    #>>> Extract target bb coords.
    win_coords = []
    for i in range(len(win_combs)):
        win_comb = win_combs[i]
        metal_coord = metal_coords[i]
        _coord = get_coords(ss.target, win_comb, metal_coord)
        win_coords.append(_coord)
    db_coord_all =  np.concatenate((db_coords, win_coords))
    inds = list(range(len(db_coords), len(db_coord_all)))

    all_mems = cluster_search(db_coords, db_coord_all, ss.rmsd)

    for i in range(len(inds)):
        ind = inds[i]
        win_comb = win_combs[i]
        mems = all_mems[ind]
        
    return


