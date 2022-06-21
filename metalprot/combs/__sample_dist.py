'''
Combs searching via a distance based method.
'''


import pandas as pd
import numpy as np
import prody as pr
from sklearn.neighbors import NearestNeighbors
import itertools
from datetime import datetime
from scipy.sparse import csr_matrix


def prepare_ligand_info(ligand, sel_cg_names):
    '''
    Prepare ligand info.
    ligand: prody atomgroup object.
    '''
    ligand_cg_dist_dict = {}
    for i, j in itertools.combinations(range(sel_cg_names.shape[0]), 2):
        ligand_cg_dist_dict[(i, j)] = []
        for ind in range(sel_cg_names.shape[1]):
            ni = sel_cg_names[i][ind]
            nj = sel_cg_names[j][ind]
            dist = pr.calcDistance(ligand.select('name ' + ni), ligand.select('name ' + nj))
            ligand_cg_dist_dict[(i, j)].append(dist[0])

    return ligand_cg_dist_dict


def set_cg_contact_atom_neighbors(cg_dict, ligand_cg_vs_atom_dict, cg_count = 3):
    '''
    For each ligand, we select 3 contacting atom, for each atom we extract the all_coord.
    '''
    print('Setting contact atom neighbors ...')

    all_coord_dict = {}
    df_contact_dict = {}
    for i in range(cg_count):

        frames = []
        for k in cg_dict.keys():
            ligand_cg_vs_atom = ligand_cg_vs_atom_dict[k]
            resnames = [ligand_cg_vs_atom[i][0]]
            df = cg_dict[k]
            if k == 'bb_cco':
                resnames = ['GLY', 'ALA', 'PRO']
                df_contact = df[
                    (df['chain'] == 'Y')  
                    & ( (df['resname'] == resnames[0]) | (df['resname'] == resnames[1]) | (df['resname'] == resnames[2]) )
                     & (df['name'] == ligand_cg_vs_atom[i][2]) 
                ]
            else:
                df_contact = df[
                    (df['chain'] == 'Y')  
                    & (df['resname'] == resnames[0]) 
                    & (df['name'] == ligand_cg_vs_atom[i][2]) 
                ]
            frames.append(df_contact)

        df_contacts = pd.concat(frames)
        df_contact_dict[i] = df_contacts
        all_coords = df_contacts[['c_x', 'c_y', 'c_z']].values
        all_coord_dict[i] = all_coords

    count = len(all_coord_dict[0])
    return all_coord_dict, df_contact_dict, count



def clash_filter(template, all_coords, dist_cut = 2.8):
    # clashing: any atom close to the bb let's say dist <= dist_cut
    bb_coords = template.dataframe[['c_x', 'c_y', 'c_z']].values
    nn = NearestNeighbors(radius= dist_cut).fit(bb_coords)
    x_in_y = nn.radius_neighbors_graph(all_coords).toarray().astype(bool)
    clashing = np.any(x_in_y, axis=1)

    #np.sum(clashing) #The clashing is too less to be considered.
    return clashing



def prepare_adj_matrix(keys, all_coord_dict, ligand_cg_dist_dict, dist_tol = 0.2):
    '''
    For each ligand, we select 3 contacting atom, for each atom we extract the all_coord.
    
    Input:
    keys = {(0, 1), (0, 2), (1, 2) etc}
    ligand_cg_dist_dict = {
        (0, 1):[2.2, 5.2, 6.2]
        (0, 2):[2.2, 5.2, 6.2]
        (1, 2):[2.2, 5.2, 6.2]
        }
    
    Return: 
    adj_matrix_dict {
        0: {(0, 1): matrix, (0, 2): matrix, (1, 2): matrix},
        1: {(0, 1): matrix, (0, 2): matrix, (1, 2): matrix},
        2: {(0, 1): matrix, (0, 2): matrix, (1, 2): matrix}
    }

    '''
    adj_matrix_dict = {}
    for i in range(3):
        all_coords = all_coord_dict[i]
        adj_matrix_d = {}
        for key in keys:   
            dist = ligand_cg_dist_dict[key][i]

            dista = dist + dist_tol
            distb = dist - dist_tol
            nbr_a = NearestNeighbors(radius= dista).fit(all_coords)
            adj_matrix_a = nbr_a.radius_neighbors_graph(all_coords).astype(bool)
            nbr_b = NearestNeighbors(radius= distb).fit(all_coords)
            adj_matrix_b = nbr_b.radius_neighbors_graph(all_coords).astype(bool)
            adj_matrix = (adj_matrix_a > adj_matrix_b) + (adj_matrix_a < adj_matrix_b)
            
            adj_matrix_d[key] = adj_matrix
            
        adj_matrix_dict[i]=adj_matrix_d

    return adj_matrix_dict


def filter_mask(keys, count, df_contact_dict, ligand_cg_atm_dict):
    '''
    Based on the connection between each cg.

    '''
    label_matrixs = []
    label_matrix_dict = {}
    for key in keys:
        labels = np.zeros(count, dtype = bool)
        #TO DO: something wrong here.
        if df_contacts['name'].iloc[i] in ligand_cg_atm_dict[0]:
            labels[i] = True

        label_matrix = np.broadcast_to(labels, (count, count))
        label_matrixs.append(label_matrix)

    for i, j in keys:
        label_matrix_dict[key] = label_matrixs[i].T * label_matrixs[j]

    return label_matrix_dict


def filter_position_mask(count, df_contact_dict, cg_count = 3):
    '''
    Filter connections between vdm on the same bb position.
    '''
    
    poss = np.array(df_contact_dict[0]['seg_chain_resnum'])
    ty_ind = 0
    ty_dict = {} 
    poss_ind = np.zeros(count, dtype = int)
    for i in range(count):
        pos_curr = poss[i]
        if pos_curr not in ty_dict.keys():
            ty_ind += 1
            ty_dict[pos_curr] = ty_ind
        poss_ind[i] = ty_dict[pos_curr]
    print(poss_ind)
    pos_matrix = np.zeros((count, count), dtype = bool)
    for i in range(count):
        pos_matrix[i] = poss_ind == poss_ind[i]
    return ~pos_matrix


def modify_adj_matrix(keys, adj_matrix_dict, pos_matrix, label_matrix_dict = None):
    '''
    adj_matrix_dict has key ((0, 1), 0)
    '''
    m_adj_matrix_dict = {}
    for k in adj_matrix_dict.keys():
        m_adj_matrix_dict[k] = {}
        for key in keys:
            adj_matrix = adj_matrix_dict[k][key]
            m_adj_matrix = adj_matrix.multiply(pos_matrix)

            if label_matrix_dict:
                label_matrix = label_matrix_dict[key]
                m_adj_matrix = m_adj_matrix.multiply(label_matrix)

            m_adj_matrix_dict[k][key] = m_adj_matrix.toarray()

    return m_adj_matrix_dict


def calc_adj_matrix_paths(m_adj_matrix_d, count, path_order_matrix, num_iter =3):
    '''
    Get the paths. Each path represent a solution candicate.
    '''
    paths = []
    for r in range(count):
        if not path_order_matrix[0][r]:
            continue
        comb = [r]
        calc_adj_matrix_paths_helper(m_adj_matrix_d, path_order_matrix, comb, num_iter, 0, paths)          
    return paths


def calc_adj_matrix_paths_helper(m_adj_matrix_d, path_order_matrix, comb, num_iter, ir, paths):
    if len(comb) >= num_iter:
        paths.append(comb)
        return
    r = comb[ir]
    r_nexts = np.where(m_adj_matrix_d[(ir, ir+1)][r])[0]
    if len(r_nexts) <= 0:
        return 
    for r_next in r_nexts:
        if not path_order_matrix[ir+1][r_next]:
            continue
        _comb = comb.copy()
        _comb.append(r_next)
        exist = True
        for i in range(len(_comb) - 2):
            if not m_adj_matrix_d[i, ir+1][_comb[i], r_next]:
                exist = False
        if exist:
            calc_adj_matrix_paths_helper(m_adj_matrix_d, path_order_matrix, _comb, num_iter, ir +1, paths)
    return 


def get_path_order_matrix(df_contact_dict, count, cg_count = 3):
    '''
    The paths have orders corespond to the order of chemical groups selected.
    '''
    path_order_matrix = np.zeros((3, count), dtype = bool)
    df = df_contact_dict[0]
    cg_types = df['CG_type']

    for j in range(count):
        if cg_types.iloc[j] == 'hie':
            path_order_matrix[0,j] = True
        if cg_types.iloc[j]  == 'hid':
            path_order_matrix[1,j] = True
        if cg_types.iloc[j]  == 'bb_cco':
            path_order_matrix[2,j] = True
    return path_order_matrix
            

def run_calc_paths(m_adj_matrix_dict, count, path_order_matrix):
    '''
    For each vdM, we select 3 contact atoms and for each ct_atom, we calc path based on the m_aj_matrix_d.
    '''
    path_dict = {}
    for k in m_adj_matrix_dict.keys():
        m_adj_matrix_d = m_adj_matrix_dict[k]
        paths = calc_adj_matrix_paths(m_adj_matrix_d, count, path_order_matrix, num_iter =3)

        for path in paths:
            if tuple(path) in path_dict.keys():
                path_dict[tuple(path)] +=1
            else: 
                path_dict[tuple(path)] =1

    filtered_paths = [k for k in path_dict.keys() if path_dict[k] == 3]
    return filtered_paths, path_dict


def superimpose_ligand_2_path(ligand, ligand_cgs_ag, filtered_paths, df_contact_dict, rmsd_cut = 0.5):
    pdbs = []
    for path in filtered_paths:
        _ligand_cgs_ag = ligand_cgs_ag.copy()

        #Extract the coord for path
        coords = [] 
        
        for p in path:  
            for i in range(3): 
                coords.append(df_contact_dict[i].iloc[p][['c_x', 'c_y', 'c_z']].values)

        transformation = pr.calcTransformation(_ligand_cgs_ag.getCoords(), np.array(coords, dtype=float))
        transformation.apply(_ligand_cgs_ag)
        rmsd = pr.calcRMSD(_ligand_cgs_ag.getCoords(), np.array(coords, dtype=float))


        _ligand = ligand.copy()           
        transformation.apply(_ligand)
        
        #TO DO: Clash filter in 'combs2'
        pdbs.append((_ligand, path, rmsd))
    
    return pdbs



