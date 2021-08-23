import numpy as np
from scipy.spatial.distance import cdist
import itertools
import prody as pr
from . import constant

def get_ABPLE(resn, phi, psi):
    try:
        psi = int(np.ceil(psi / 10.0)) * 10
        phi = int(np.ceil(phi / 10.0)) * 10
        if psi == -180:
            psi = -170
        if phi == -180:
            phi = -170
        return constant.abple_dict[resn][psi][phi]
    except ValueError:
        return 'n'

def seq_get_ABPLE(target):
    abples = []
    phipsi = []
    for resn in target.iterResidues():
        try:
            phi = pr.calcPhi(resn)
            psi = pr.calcPsi(resn)
            phipsi.append((phi, psi))

            ap = get_ABPLE(resn.getResname(), phi, psi)
            abples.append(ap)
        except:
            phipsi.append((0, 0))
            abples.append('n')
    #print(abples)
    return abples, phipsi
    

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

def filter_phipsi(angle, query_angle, filter_angle_val = 30):
    '''
    phi psi ranges is between [-180, 180]
    calc the range of phi or psi angle with a filter range.
    Ex. phi = 170, filter range = +-30, [140, 200] as 200 > 180. The true range is [(140, 180), (0, 20)]
    The filter_angle_val must be < 180
    '''
    filter_ranges = []

    low = angle - filter_angle_val

    high = angle + filter_angle_val

    if low < -180:
        filter_ranges.append((-180, high))
        filter_ranges.append((low + 360, 180))

    elif high > 180:
        filter_ranges.append((low, 180))
        filter_ranges.append((-180, high - 360))
    
    else:
        filter_ranges.append((low, high))

    for r in filter_ranges:
        if query_angle >= r[0] and query_angle <= r[1]:
            return True
    return False

class CombGraph:
    def __init__(self, inds, pair_set, num_iter = 3):
        self.inds = inds
        self.ind_len = len(inds)

        self.num_iter = num_iter

        self.pair_set = pair_set
        self.all_paths = []


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

        rs = []
        for i in range(r + 1, self.ind_len):
            if (r, i) in self.pair_set:
                rs.append(i)
        if len(rs) == 0:
            return
        for _r in rs:           
            satisfy = True
            for ci in range(len(temp)-1):
                tv = temp[ci]
                if not (tv, _r) in self.pair_set:
                    satisfy = False
            if satisfy:      
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

    graph = CombGraph(inds, pair_set, num_iter)
    graph.get_paths()

    return graph.all_paths


