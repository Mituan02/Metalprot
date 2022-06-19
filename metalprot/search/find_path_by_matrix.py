import numpy as np
import itertools
from scipy.sparse import lil_matrix
from sklearn.neighbors import NearestNeighbors, radius_neighbors_graph

from ..basic import utils
from ..basic.constant import one_letter_code


'''
test_matrix = np.array(
    [[0, 0, 1, 1],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [0, 0, 0, 0]],
    dtype = bool
)

test_paths = calc_adj_matrix_paths(test_matrix, num_iter)

paths = calc_adj_matrix_paths(m_adj_matrix)

for i in range(len(ss.vdms)):
    if '5od1' in ss.vdms[i].query.getTitle():
        print(i)
        print(ss.vdms[i].query.getTitle())
index = [3776, 4387*1+2865, 4387*2+2192]
index == [3776, 7252, 10966]

all_vdms = []
for i in range(len(wins)):
    all_vdms.extend(ss.vdms)
all_vdms[7252].query.getTitle()


adj_matrix[3776, 7252]
adj_matrix[3776, 10966]
adj_matrix[7252, 10966]

win_mask[3776, 7252]
win_mask[3776, 10966]
win_mask[7252, 10966]


paths = calc_adj_matrix_paths(m_adj_matrix)

'''

def calc_adj_matrix_paths(m_adj_matrix, num_iter =3):
    paths = []
    for r in range(m_adj_matrix.shape[0]):
        inds = m_adj_matrix.rows[r]
        if len(inds) < num_iter-1:
            continue
        for _comb in itertools.combinations(inds, num_iter-1):
            comb = [r]
            comb.extend(_comb)
            valid = calc_adj_matrix_paths_helper(m_adj_matrix, comb, num_iter, 1)
            if valid:
                #print(comb)
                paths.append(comb)

    return paths


def calc_adj_matrix_paths_helper( m_adj_matrix, comb, num_iter, iter):
    if iter >= num_iter -1:
        return True
    r_curr = comb[iter]
    for r_next in comb[iter+1:]:
        if not m_adj_matrix[r_curr, r_next]:
            return False
    return calc_adj_matrix_paths_helper(m_adj_matrix, comb, num_iter, iter +1)


def neighbor_generate_nngraph(ss):
    '''
    Instead of doing this in a pairwise way as the function 'search.neighbor_generate_pair_dict'.
    Here we calc one nearest neighbor object and graph.

    ss is the search.Search_vdM object.
    '''
    wins = sorted(list(ss.neighbor_query_dict.keys()))
    metal_vdm_size = len(ss.all_metal_vdm.get_metal_mem_coords())

    all_coords = []
    win_labels = []
    vdm_inds = []
    for inx in range(len(wins)):
        wx = wins[inx] 
        win_labels.extend([wx]*metal_vdm_size)
        n_x = ss.neighbor_query_dict[wx].get_metal_mem_coords()
        all_coords.extend(n_x)
        vdm_inds.extend(range(metal_vdm_size))

    bb_coords = ss.target.select('name N C CA O').getCoords()

    mc_coords = []
    for inx in range(len(wins)):
        wx = wins[inx]
        n_x = ss.neighbor_query_dict[wx].get_metalcontact_mem_coords()
        mc_coords.extend(n_x)

    # sc_coords = []
    # sc_coord_inds = []
    # for inx in range(len(wins)):
    #     wx = wins[inx]
    #     n_x, n_x_id = ss.neighbor_query_dict[wx].get_sc_mem_coords_and_ids()
    #     sc_coords.extend(n_x)
    #     sc_coord_inds.extend(n_x_id)


    #>>> create mask (The method is not used.)
    # mask = generate_filter_mask(ss, wins, win_labels, metal_vdm_size, adj_matrix_bb)
    # #calc modified adj matrix
    # m_adj_matrix = adj_matrix.multiply(mask)
    print('generate_mask_labels.')
    mask_labels = generate_mask_labels(ss, wins, metal_vdm_size, bb_coords, all_coords)
    print('filter_adj_matrix by label.')
    m_adj_matrix = filter_adj_matrix(ss, metal_vdm_size, all_coords, mask_labels)
    print('filter_adj_matrix_by_mc_min.')
    #m_adj_matrix = filter_adj_matrix_by_mc(ss, mc_coords, m_adj_matrix)
    m_adj_matrix = filter_adj_matrix_by_mc_min(ss, wins, metal_vdm_size, mc_coords, m_adj_matrix)
    print('filter_adj_matrix_by_mc_max.')
    m_adj_matrix = filter_adj_matrix_by_mc_max(ss, wins, metal_vdm_size, mc_coords, m_adj_matrix)
    print(m_adj_matrix.shape)
    return m_adj_matrix.tolil(), win_labels, vdm_inds

    
def generate_mask_labels(ss, wins, metal_vdm_size, bb_coords, all_coords):
    '''
    One issue for the mask method is that the matrix is sparse, there will be a lot of unnecessary calculation.
    For example, filter phi psi, if a vdm is never be in any neighbor, there is no need to calc it.
    '''

    #>>> metal_clashing with bb
    nbrs_bb = NearestNeighbors(radius= 3.5).fit(all_coords)
    adj_matrix_bb = nbrs_bb.radius_neighbors_graph(bb_coords).astype(bool)
    #print(adj_matrix_bb.shape)

    mask_labels = np.ones(len(wins)*metal_vdm_size, dtype=bool)

    # Metal bb Clashing filter
    for i in range(adj_matrix_bb.shape[0]):
        mask_labels *= ~(adj_matrix_bb[i].toarray().reshape(len(mask_labels),))

    #aa origin filter. 
    if ss.validateOriginStruct:
        v_aa = np.array([v.aa_type for v in ss.vdms])
        ress = [one_letter_code[ss.target.select('name CA and resindex ' + str(wx)).getResnames()[0]] for wx in wins]
        aa_labels = np.zeros(len(wins)*metal_vdm_size, dtype=bool)
        for inx in range(len(wins)):
            aa_labels[inx*metal_vdm_size:(inx+1)*metal_vdm_size] = v_aa == ress[inx]
        mask_labels *= aa_labels

    # filter vdM by score:
    if ss.search_filter.filter_vdm_score:
        v_scores = np.array([v.score for v in ss.vdms])
        vdm_score_labels = np.zeros(len(wins)*metal_vdm_size, dtype=bool)
        for inx in range(len(wins)):
            vdm_score_labels[inx*metal_vdm_size:(inx+1)*metal_vdm_size] = v_scores >= ss.search_filter.min_vdm_score
        mask_labels *= vdm_score_labels

    # filter vdM by count:
    if ss.search_filter.filter_vdm_count:
        v_count = np.array([v.clu_num for v in ss.vdms])
        vdm_count_labels = np.zeros(len(wins)*metal_vdm_size, dtype=bool)
        for inx in range(len(wins)):
            vdm_count_labels[inx*metal_vdm_size:(inx+1)*metal_vdm_size] = v_count >= ss.search_filter.min_vdm_clu_num
        mask_labels *= vdm_count_labels

    # abple filter
    if ss.search_filter.filter_abple:
        v_abples = np.array([v.abple for v in ss.vdms])
        apxs = [ss.target_abple[wx] for wx in wins]
        abple_labels = np.zeros(len(wins)*metal_vdm_size, dtype=bool)
        for inx in range(len(wins)):
            abple_labels[inx*metal_vdm_size:(inx+1)*metal_vdm_size] = v_abples == apxs[inx]
        mask_labels *= abple_labels

    #Filter unwanted amino acids. if ss.allowed_aas = {'H', 'D'}, then {'E', 'S'} will be eliminated.
    if not ss.validateOriginStruct and len(ss.allowed_aas) > 0 and len(ss.allowed_aas) < 4:
        v_aa = np.array([v.aa_type for v in ss.vdms])
        aa_allow_labels = np.zeros(len(wins)*metal_vdm_size, dtype=bool)
        for inx in range(len(wins)):
            aa_allow_labels[inx*metal_vdm_size:(inx+1)*metal_vdm_size] = np.array([v in ss.allowed_aas for v in v_aa])
        mask_labels *= aa_allow_labels

    #phi psi filter
    if ss.search_filter.filter_phipsi:
        #TO DO: filter phi psi need to be changed to be able to broadcast.
        v_phis = [v.phi for v in ss.vdms]
        v_psis = [v.psi for v in ss.vdms]

        phis = [ss.phipsi[wx][0] for wx in wins]
        psis = [ss.phipsi[wx][1] for wx in wins]
        phipsi_labels = np.zeros(len(wins)*metal_vdm_size, dtype=bool)
        for inx in range(len(wins)):  
            for i in range(metal_vdm_size):
                phi_ok = utils.filter_phipsi(phis[inx], v_phis[i], ss.search_filter.max_phipsi_val)
                psi_ok = utils.filter_phipsi(psis[inx], v_psis[i], ss.search_filter.max_phipsi_val)
                if phi_ok and psi_ok:
                    phipsi_labels[inx*metal_vdm_size + i] = True
        mask_labels *= phipsi_labels

    return mask_labels


def filter_adj_matrix(ss, metal_vdm_size, all_coords, mask_labels):
    '''
    One issue for the mask method is that the matrix is sparse, there will be a lot of unnecessary calculation.
    For example, filter phi psi, if a vdm is never be in any neighbor, there is no need to calc it.
    '''     
    #>>> calc radius_neighbors_graph for metal-metal.
    #nbrs = NearestNeighbors(radius= ss.metal_metal_dist).fit(all_coords)
    #adj_matrix = nbrs.radius_neighbors_graph(all_coords).astype(bool)
    adj_matrix = radius_neighbors_graph(all_coords, radius= ss.metal_metal_dist).astype(bool)

    #print(adj_matrix.shape)
    #>>> Final output matrix
    m_adj_matrix = lil_matrix(adj_matrix.shape, dtype=bool)
    for r in range(adj_matrix.shape[0]):
        if not mask_labels[r]:
            continue 

        for ind in range(adj_matrix.indptr[r], adj_matrix.indptr[r+1]):
            c = adj_matrix.indices[ind]
            
            #>>> vdm on the same position don't connect.     
            if c < r or (not mask_labels[c]) or r//metal_vdm_size == c//metal_vdm_size:
                continue

            m_adj_matrix[r, c] = True

    return m_adj_matrix


def filter_adj_matrix_by_mc_min(ss, wins, metal_vdm_size, mc_coords, m_adj_matrix):
    '''
    The method is used to filter metalContactingAtom-metalContactingAtom distance low boundry.
    '''
    #>>> The method used below to calculate the whole adj_matrix_mc is not very efficient. 
    #>>> Considering the radius is large and mc_coords from the same aa position will have a lot True value.
    #adj_matrix_mc_min = radius_neighbors_graph(mc_coords, radius= ss.search_filter.pair_aa_aa_dist_range[0], mode = 'connectivity', includes_self = False).astype(bool)
    
    #>>> The new method is to calculate the adj_matrix in a separate way.
    for i in range(len(wins)):
        for j in range(len(wins)):
            if i ==j:
                continue
            xs = mc_coords[i*metal_vdm_size: (i+1)*metal_vdm_size]
            ys = mc_coords[j*metal_vdm_size: (j+1)*metal_vdm_size]
            nbrs_mc = NearestNeighbors(radius= ss.search_filter.pair_aa_aa_dist_range[0]).fit(ys)
            adj_matrix_mc = nbrs_mc.radius_neighbors_graph(xs).astype(bool)
            m_adj_matrix[i*metal_vdm_size: (i+1)*metal_vdm_size, j*metal_vdm_size: (j+1)*metal_vdm_size] = m_adj_matrix[i*metal_vdm_size: (i+1)*metal_vdm_size, j*metal_vdm_size: (j+1)*metal_vdm_size] > adj_matrix_mc

    return m_adj_matrix


def filter_adj_matrix_by_mc_max(ss, wins, metal_vdm_size, mc_coords, m_adj_matrix):
    '''
    The method is used to filter metalContactingAtom-metalContactingAtom distance low boundry.
    '''
    for i in range(len(wins)):
        for j in range(len(wins)):
            if i ==j:
                continue
            xs = mc_coords[i*metal_vdm_size: (i+1)*metal_vdm_size]
            ys = mc_coords[j*metal_vdm_size: (j+1)*metal_vdm_size]
            nbrs_mc = NearestNeighbors(radius= ss.search_filter.pair_aa_aa_dist_range[1]).fit(ys)
            adj_matrix_mc = nbrs_mc.radius_neighbors_graph(xs).astype(bool).tolil()
            m_adj_matrix[i*metal_vdm_size: (i+1)*metal_vdm_size, j*metal_vdm_size: (j+1)*metal_vdm_size] = m_adj_matrix[i*metal_vdm_size: (i+1)*metal_vdm_size, j*metal_vdm_size: (j+1)*metal_vdm_size].multiply(adj_matrix_mc) 

    return m_adj_matrix


#>>> The funciton is deprecated.
def generate_filter_mask(ss, wins, win_labels, metal_vdm_size, adj_matrix_bb, adj_matrix_mc):
    '''
    One issue for the mask method is that the matrix is sparse, there will be a lot of unnecessary calculation.
    For example, filter phi psi, if a vdm is never be in any neighbor, there is no need to calc it.
    '''
        
    wins = sorted(list(ss.neighbor_query_dict.keys()))
    metal_vdm_size = len(ss.all_metal_vdm.get_metal_mem_coords())    
    
    win_mask = np.ones((len(win_labels), len(win_labels)), dtype=bool)
    win_mask = np.triu(win_mask)

    # win filter: vdm on the same position don't connect.
    for inx in range(len(wins)):
        win_mask[inx*metal_vdm_size:(inx+1)*metal_vdm_size, inx*metal_vdm_size:(inx+1)*metal_vdm_size] = 0

    # Metal bb Clashing filter.
    metal_clashing_vec = np.ones(len(win_labels), dtype=bool)
    for i in range(adj_matrix_bb.shape[0]):
        metal_clashing_vec *= ~(adj_matrix_bb[i].toarray().reshape(len(win_labels),))
    labels_m = np.broadcast_to(metal_clashing_vec, (len(wins)*metal_vdm_size, len(wins)*metal_vdm_size))
    win_mask *= labels_m.T
    win_mask *= labels_m

    # Metal contacting atom clashing filter.
    mc_clashing_vec = np.ones(len(win_labels), dtype=bool)
    for i in range(adj_matrix_mc.shape[0]):
        mc_clashing_vec *= ~(mc_clashing_vec[i].toarray().reshape(len(win_labels),))
    labels_mc = np.broadcast_to(mc_clashing_vec, (len(wins)*metal_vdm_size, len(wins)*metal_vdm_size))
    win_mask *= labels_mc.T
    win_mask *= labels_mc

    # Modify mask with aa filter. 
    if ss.validateOriginStruct:
        v_aa = np.array([v.aa_type for v in ss.vdms])
        ress = [one_letter_code[ss.target.select('name CA and resindex ' + str(wx)).getResnames()[0]] for wx in wins]
        labels = np.zeros(len(wins)*metal_vdm_size, dtype=bool)
        for inx in range(len(wins)):
            labels[inx*metal_vdm_size:(inx+1)*metal_vdm_size] = v_aa == ress[inx]
        labels_m = np.broadcast_to(labels, (len(wins)*metal_vdm_size, len(wins)*metal_vdm_size))
        win_mask *= labels_m.T
        win_mask *= labels_m

    if ss.search_filter.filter_abple:
        v_abples = np.array([v.abple for v in ss.vdms])
        apxs = [ss.target_abple[wx] for wx in wins]
        labels = np.zeros(len(wins)*metal_vdm_size, dtype=bool)
        for inx in range(len(wins)):
            labels[inx*metal_vdm_size:(inx+1)*metal_vdm_size] = v_abples == apxs[inx]
        labels_m = np.broadcast_to(labels, (len(wins)*metal_vdm_size, len(wins)*metal_vdm_size))
        win_mask *= labels_m.T
        win_mask *= labels_m      

    #Filter unwanted amino acids. if ss.allowed_aas = {'H', 'D'}, then {'E', 'S'} will be eliminated.
    if not ss.validateOriginStruct and len(ss.allowed_aas) > 0 and len(ss.allowed_aas) < 4:
        v_aa = np.array([v.aa_type for v in ss.vdms])
        labels = np.zeros(len(wins)*metal_vdm_size, dtype=bool)
        for inx in range(len(wins)):
            labels[inx*metal_vdm_size:(inx+1)*metal_vdm_size] = np.array([v in ss.allowed_aas for v in v_aa])
        labels_m = np.broadcast_to(labels, (len(wins)*metal_vdm_size, len(wins)*metal_vdm_size))
        win_mask *= labels_m.T
        win_mask *= labels_m 

    if ss.search_filter.filter_phipsi:
        #TO DO: filter phi psi need to be changed to be able to broadcast.
        v_phis = [v.phi for v in ss.vdms]
        v_psis = [v.psi for v in ss.vdms]

        phis = [ss.phipsi[wx][0] for wx in wins]
        psis = [ss.phipsi[wx][1] for wx in wins]
        labels = np.zeros(len(wins)*metal_vdm_size, dtype=bool)
        for inx in range(len(wins)):  
            for i in range(metal_vdm_size):
                #TO DO: how to ignore unnecessary psiphi
                #if any(win_mask[inx*metal_vdm_size + i, ]) and any(adj_matrix[inx*metal_vdm_size + i, ]):
                phi_ok = utils.filter_phipsi(phis[inx], v_phis[i], ss.search_filter.max_phipsi_val)
                psi_ok = utils.filter_phipsi(psis[inx], v_psis[i], ss.search_filter.max_phipsi_val)
                if phi_ok and psi_ok:
                    labels[inx*metal_vdm_size + i] = True
        labels_m = np.broadcast_to(labels, (len(wins)*metal_vdm_size, len(wins)*metal_vdm_size))
        win_mask *= labels_m.T
        win_mask *= labels_m  
            
    return win_mask
