import os
from typing import Dict
import numpy as np
import itertools
from numpy.lib.function_base import extract
import prody as pr
from prody.measure.transform import superpose
from prody.proteins.pdbfile import writePDB
from scipy.spatial.distance import cdist, dice
import datetime
from .ligand_database import clu_info
from .extract_vdm import get_vdm_mem
from .quco import Query, Comb, pair_wist_geometry
from . import core
from . import hull
from . import utils

from sklearn.neighbors import NearestNeighbors


class Node:
    def __init__(self, wid, rid, win_inds, all_pos_len):
        self.wid = wid
        self.rid = rid
        self.all_pos_len = all_pos_len

        self.win_dict = []
        for w in win_inds:
            cnn = [False for i in range(all_pos_len)]
            self.win_dict.append(cnn)

class Graph:
    def __init__(self, wins, all_pos_len):
        self.wins = wins
        self.win_inds = list(range(len(wins)))
        self.all_pos_len = all_pos_len

        self.nodes = []
        for i in range(len(wins)):
            nodes = []
            for j in range(self.all_pos_len):
                nodes.append(Node((i, j), self.all_pos_len))
        self.all_paths = []


    def calc_pair_connectivity(self, neighbor_pair_dict):

        for c, c2 in itertools.permutations(self.win_inds):
            wx = self.wins[c]
            wy = self.wins[c2]    
            for pair in neighbor_pair_dict[(wx, wy)]:
                for r in range(self.all_pos_len):              
                    for y in pair[wy]:

                        self.nodes[c][r].win_dict[c2][y] = True

        return

    
    def get_paths(self):

        c = 0
        for r in range(self.all_pos_len):
            self.get_path_helper(c, r, [r])
        return


    def get_path_helper(self, c, r, temp):
        '''
        Dynamic programming. 
        '''
        #if c == self.num_iter -1:
        if len(temp) == len(self.wins):
            self.all_paths.append([t for t in temp])
            return
        #print('c ' + str(c) + ' r ' + str(r))
        rs = [i for i, x in enumerate(self.nodes[c][r].win_dict[c+1]) if x]
        if len(rs) == 0:
            return      
        for _r in rs:
            #print('col ' + str(c) + ' temp ' + '_'.join([str(t) for t in temp]))
            for ci in range(len(temp)-1):
                tv = temp[ci]
                if not self.nodes[c+1][_r].win_dict[ci][tv]:
                    return
            _temp = temp.copy()
            _temp.append(_r)
            self.get_path_helper(c+1, _r, [t for t in _temp])

        return

            



class Search_vdm:
    '''
    The function to search comb
    '''
    def __init__(self, target_pdb, workdir, querys, all_metal_query, num_iter, rmsd = 0.25, win_filtered = None, 
    contact_querys = None, secondshell_querys = None, validateOriginStruct = False):

        if workdir:
            _workdir = os.path.realpath(workdir)
            if not os.path.exists(_workdir):
                os.mkdir(_workdir)
        else:
            _workdir = os.getcwd() + '/output_' + datetime.now().strftime('%Y-%m-%d-%H-%M-%S')          
            os.mkdir(_workdir)
        self.workdir = _workdir + '/'

        self.target = pr.parsePDB(target_pdb)

        self.num_iter = num_iter

        self.dist_array, self.id_array, self.dists = utils.get_contact_map(self.target, win_filtered)

        self.rmsd = rmsd

        self.win_filtered = win_filtered

        self.validateOriginStruct = validateOriginStruct

        #neighbor searching strategy---------- 
        self.querys = querys             
        self.all_metal_query = all_metal_query #The query with all_metal_coord_ag
        self.cluster_dict = {} # {metal_id 1234: (HIS cluster 0)} 

        self.neighbor_query_dict = dict() # {93: [the only centroid query with all metal coords]}
        self.neighbor_pair_dict = dict() # {(33, 37): [xs-33 near 37 coords]}
        self.neighbor_comb_dict = dict() # {(33, 37, 42): } Please check neighbor_win2comb()
        
        #contact-----------------------
        self.contact_querys = contact_querys

        #secondshell-----------------------
        self.secondshell_querys = secondshell_querys
        #end---------------------------- 

    #region Neighbor Search
    '''
    A new searching method based on sklearn NerestNeighbor radius_neighbors.
    The method is similarly used in COMBS by Nick Polizzi.
    One difference is that, we need to consider the overlap of metal coords from >3 aa sidechains. 
    While for general ligand, one chemical group of aa bind one chemical group of ligand. 
    '''
    def run_neighbor_search(self):
        '''
        All functions need to run the neighbor search.
        '''
        print('run_neighbor_search')

        #TO DO: where should I apply filters: win filter, query_metal filter, phipsi, etc.
        self.neighbor_generate_query_dict()

        self.neighbor_generate_pair_dict()

        self.neighbor_search_wins()



        return


    def neighbor_generate_query_dict(self):
        '''
        return self.neighbor_query_dict = dict() # {93: [the only centroid query with all metal coords]}

        '''
        print('hull_generate_query_dict')
        wins = []
        if self.win_filter:
            wins.extend([w for w in self.win_filter])
        else:
            t = self.target.select('name CA').getResindices()
            wins.extend(([w for w in t]))

        for w in wins:     
            cquery = self.supperimpose_target_bb(align_sel='name N CA C')
            self.neighbor_query_dict[w] = cquery
        return

        
    def supperimpose_target_bb(self, win, align_sel='name N CA C'):
        '''
        Copy the all_metal_query to a new object.
        Transform the copied all_metal_query to the target win. 
        '''
        _query = self.all_metal_query.copy()
        target_sel = self.target.select('resindex ' + str(win))

        if len(_query.query.select(align_sel)) != len(target_sel.select(align_sel)):
            return None
        
        transform = pr.calcTransformation(_query.query.select(align_sel), target_sel.select(align_sel))
        transform.apply(_query.query)
        transform.apply(_query.hull_ag)       

        return _query


    def neighbor_generate_pair_dict(self):
        '''
        To generate pair win dict, we will first filter out the impossible win pair by distance.
        Apply utils.check_pair_distance_satisfy()

        '''

        print('neighbor_generate_pair_dict')
        
        wins = list(self.hull_query_dict.keys())

        for inx in range(len(wins)):
            for iny in range(inx + 1, len(wins)):
                wx = wins[inx]
                wy = wins[iny]  
                dist_ok, ws = utils.check_pair_distance_satisfy(wx, wy, self.dists)
                if not dist_ok:
                    continue
                # if self.query_all_metal and not check_hull_satisfy(self.target, wx, wy, self.query_all_metal_x, self.query_all_metal_y):
                #     continue
                n_x = self.neighbor_query_dict[wx].get_hull_points()
                n_y = self.neighbor_query_dict[wy].get_hull_points()

                x_in_y, x_has_y = self.calc_pairwise_neighbor(n_x, n_y, self.rmsd)
                y_in_x, y_has_x = self.calc_pairwise_neighbor(n_y, n_x, self.rmsd)
                if x_has_y and y_has_x:
                    self.neighbor_pair_dict[(wx, wy)] = (x_in_y)
                    self.neighbor_pair_dict[(wy, wx)] = (y_in_x)
        return


    def calc_pairwise_neighbor(self, n_x, n_y, rmsd = 0.25):
        '''
        Use sklean NearestNeighbors
        '''

        neigh_y = NearestNeighbors(radius= rmsd) 
        neigh_y.fit(n_y)

        x_in_y = neigh_y.radius_neighbors(n_x)
        x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])

        return x_in_y, x_has_y

    
    def neighbor_search_wins(self):
        '''
        The combinations of positions are extracted from all possible positions with pairs.
        '''
        wins = sorted(set(self.hull_query_dict.keys()))

        print('All wins with overlap {}'.format(wins))

        win_combs = itertools.combinations(wins, self.num_iter)
        for win_comb in win_combs:
            print(win_comb)
            
            self.neighbor_win2comb(win_comb)

        return


    def neighbor_construct_comb(self, win_comb):
        '''
        win_comb: [0, 1, 2, 3]

        cluster_dict: {(0, 0, 0, 0): {0:[1, 3, 4], 1:[2, 3, 4], 2: [2, 6, 7], 3:[1, 2, 3]}}
        The key is win, the value is list of boolean, each represent one metal coord exist in all other wins' metal.
        
        self.neighbor_comb_dict: {(0, 1, 2, 3): cluster_dict}

        As we calculate all metals at the same time. 
        We need to extract vdm representations.
        '''
        print('neighbor_construct comb')

        for x, y in itertools.permutations(win_comb, 2):
            if (x, y) not in self.neighbor_pair_dict.keys():
                return None
        
        graph = Graph(win_comb, len(self.querys))

        graph.calc_pair_connectivity(self.neighbor_pair_dict)

        graph.get_paths()

        clu_dict = {}

        for path in graph.all_paths:

            clu_key = tuple([self.cluster_dict[p] for p in path])
            
            if clu_key in clu_dict:
                for i in range(len(win_comb)):
                    clu_dict[clu_key][win_comb[i]].add(path[i])
            else:
                clu_dict[clu_key] = {}
                for i in range(len(win_comb)):
                    clu_dict[clu_key][win_comb[i]] = set(path[i])

        if len(clu_dict) > 0:
            self.neighbor_comb_dict[tuple(win_comb)] = clu_dict

        return

    
    def neighbor_calc_comb_score(self):
        '''
        The summed vdM score could not reflect the designability.
        Here is a new score method with weight added.
        '''
        print('hull_calc_comb_score')
        for key in self.neighbor_comb_dict.keys():
            score = 0
            total = 0
            total_all = 0
            for i in range(len(key[0])):
                win = key[0][i]
                clu = key[1][i] 
                query = self.hull_query_dict[win][clu]
                score += len(query.candidates_metal_points)/len(query.get_hull_points())
                total += query.clu_num
                total_all += query.clu_total_num
            score = np.log(score*(total/total_all))
            self.hull_score_dict[key] = score
        return        


            
            

