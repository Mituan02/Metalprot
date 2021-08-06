import numpy as np
from scipy.spatial.distance import cdist
import itertools

def get_contact_map(target, win_filter = None):
    '''
    calculate contact map for 2aa_sep database.
    return the ordered distance array and resindex array.
    '''
    xyzs = []
    for c in target.select('protein and name CA').getCoords():
        xyzs.append(c)
    xyzs = np.vstack(xyzs)  
    dists = cdist(xyzs, xyzs)

    dist_array = []
    id_array = []

    for i in range(len(xyzs)):
        for j in range(i+1, len(xyzs)):
            if win_filter:
                if i not in win_filter or j not in win_filter:
                    continue
            dist_array.append(dists[i, j])  
            id_array.append((i, j))
    dist_array, id_array = zip(*sorted(zip(dist_array, id_array)))
    return dist_array, id_array, dists


def check_pair_distance_satisfy(wx, wy, dists):
    '''
    wx, wy could be int or tuple
    '''
    wins = []
    if isinstance(wx, tuple):
        wins.extend(wx)
    else:
        wins.append(wx)

    if isinstance(wy, tuple):
        wins.extend(wy)
    else:
        wins.append(wy)
    
    for _wx, _wy in itertools.combinations(wins, 2):
        _x, _y = (_wx, _wy) if _wx < _wy else (_wy, _wx) 
        if dists[_x, _y] > 15:
            return False, wins
    return True, wins

class Node:
    def __init__(self, id, all_pos_len):
        self.id = id
        self.all_pos_len = all_pos_len

        #list of bool that represent the connectivity
        #self.left = [False for i in range(all_pos_len)]
        self.right = [False for i in range(all_pos_len)]


class Graph:
    def __init__(self, inds, num_iter = 3):
        self.inds = inds
        self.ind_len = len(inds)

        self.num_iter = num_iter
        self.nodes = []
        for i in range(num_iter):
            nodes = []
            for j in range(self.ind_len):
                nodes.append(Node((i, j), self.ind_len))
            self.nodes.append(nodes)

        self.all_paths = []

    def calc_pair_connectivity(self, pair):
        '''
        For all the pair that satisfy the connection requirement. 
        Change the Node connectivity.
        '''
        x, y = pair
        for i in range(self.num_iter):
            self.nodes[i][x].right[y] = True
        return 

    def calc_pair_connectivity_all(self, all_pair):
        for pair in all_pair:
            self.calc_pair_connectivity(pair)
        return

    def get_paths(self):
        '''
        After calc connectivity, calc path with dynamic programming. 
        '''
        c = 0
        for r in range(self.ind_len):
            self.get_path_helper(c, r, [r])
        return

    def get_path_helper(self, c, r, temp):
        '''
        Dynamic programming. 
        '''
        #if c == self.num_iter -1:
        if len(temp) == self.num_iter:
            self.all_paths.append([t for t in temp])
            return
        #print('c ' + str(c) + ' r ' + str(r))
        rs = [i for i, x in enumerate(self.nodes[c][r].right) if x]
        if len(rs) == 0:
            return
        for _r in rs:
            #print('col ' + str(c) + ' temp ' + '_'.join([str(t) for t in temp]))
            _temp = temp.copy()
            _temp.append(_r)
            self.get_path_helper(c+1, _r, [t for t in _temp])

        return


#How to handle pair-wise combinations.
def combination_calc(inds, pair_set, num_iter = 3):
    '''
    inds: List of positions for combination calcluation.
    previously use of 'itertools.combinations(resinds, 3)', which is not efficient.
    Should use a graph based method with dynamic programming to remove useless combinations.
    '''
    graph = Graph(inds, num_iter)
    graph.calc_pair_connectivity_all(pair_set)
    graph.get_paths()

    return graph.all_paths


