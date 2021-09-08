import contextlib
from math import comb, log
import os
from typing import Dict
from metalprot.apps import transformation
from numba.cuda import target
import numpy as np
import itertools
from numpy.lib.function_base import extract
import prody as pr
from prody.measure.transform import superpose
from prody.proteins.pdbfile import writePDB
import scipy.spatial
from scipy.spatial.distance import cdist, dice
import datetime
from .ligand_database import clu_info
from .extract_vdm import get_vdm_mem
from . import quco
from . import core
from . import hull
from . import utils
from . import ligand_database

from sklearn.neighbors import NearestNeighbors
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool

class Graph:
    def __init__(self, wins, all_pos_len):
        self.wins = wins
        self.win_inds = list(range(len(wins)))
        self.all_pos_len = all_pos_len

        self.pair_dict = {}
        for c, c2 in itertools.permutations(self.win_inds, 2):
            self.pair_dict[(c, c2)] = set()

        self.all_paths = []


    def calc_pair_connectivity(self, neighbor_pair_dict):

        for c, c2 in itertools.permutations(self.win_inds, 2):
            wx = self.wins[c]
            wy = self.wins[c2]    
            for r in range(self.all_pos_len):
                connect = neighbor_pair_dict[(wx, wy)][r]
                #print('The len of connect in {} is {}'.format(r, len(connect)))
                if len(connect) <= 0:
                    continue     
                for y in connect:
                    self.pair_dict[(c, c2)].add((r, y))
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
        if len(temp) == len(self.wins):
            self.all_paths.append([t for t in temp])
            return
        #print('c ' + str(c) + ' r ' + str(r))
        rs = []
        for i in range(self.all_pos_len):
            if (r, i) in self.pair_dict[(c, c+1)]:
                rs.append(i)

        if len(rs) == 0:
            return      
        for _r in rs:
            #print('col ' + str(c) + ' temp ' + '_'.join([str(t) for t in temp]))
            satisfy = True
            for ci in range(len(temp)-1):
                tv = temp[ci]
                if not (tv, _r) in self.pair_dict[(ci, c+1)]:
                    satisfy = False
                    #break
            if satisfy:
                _temp = temp.copy()
                _temp.append(_r)
                self.get_path_helper(c+1, _r, [t for t in _temp])

        return


class CombInfo:
    def __init__(self):
        self.totals = None
        self.scores = None

        self.geometry = None
        self.aa_aa_pair = None
        self.metal_aa_pair = None
        self.angle_pair = None

        self.query_dict = {}
        self.volume = 0
        self.volPerMetal = 0
        self.diameter = 0
        self.centroid_dict = {}
        self.fracScore = -10.00
        self.multiScore = -10.00  #Calc multiScore (By Bill: -ln(Na/SumNa * Nb/SumNb * Nc/SumNc))

        #evaluation property for evaluation search
        self.eval_mins = None
        self.eval_min_vdMs = None
        self.eval_is_origin = False

        #After search filter property
        self.pair_aa_aa_dist_ok = 0 #0: unchecked. -1: condition unsatisfied; 1: condition satisfied.
        self.pair_angle_ok = 0
        self.pair_metal_aa_dist_ok = 0


    def calc_geometry(self):
        all_coords = []
        metal_coords = []  

        for key in self.query_dict.keys():
            coords = []
            for _query in self.query_dict[key]:                                  
                coords.append(_query.get_contact_coord())
                metal_coords.append(_query.get_metal_coord())
            all_coords.append(pr.calcCenter(hull.transfer2pdb(coords)))         
        all_coords.append(pr.calcCenter(hull.transfer2pdb(metal_coords)))

        self.geometry = hull.transfer2pdb(all_coords, ['NI' if i == len(all_coords)-1 else 'N' for i in range(len(all_coords))])
        self.aa_aa_pair, self.metal_aa_pair, self.angle_pair  = quco.pair_wise_geometry(self.geometry)
        return

    def after_search_condition_satisfied(self, pair_angle_range = None, pair_aa_aa_dist_range = None, pair_metal_aa_dist_range = None):
        '''
        range = (75, 125) for Zn.
        if all pairwise angle is between the range. The geometry is satisfied.
        '''
        if pair_angle_range:
            for an in self.angle_pair:
                if an < pair_angle_range[0] or an > pair_angle_range[1]:
                    self.pair_angle_ok = -1
                    return False
                else:
                    self.pair_angle_ok = 1
        if pair_aa_aa_dist_range:           
            for ad in self.aa_aa_pair:
                if ad < pair_aa_aa_dist_range[0] or ad > pair_aa_aa_dist_range[1]:
                    self.pair_aa_aa_dist_ok = -1
                    return False
                else:
                    self.pair_aa_aa_dist_ok = 1

        if pair_metal_aa_dist_range:
            for amd in self.metal_aa_pair:
                if amd < pair_metal_aa_dist_range[0] or amd > pair_metal_aa_dist_range[1]:
                    self.pair_aa_metal_dist_ok = -1
                    return False
                else:
                    self.pair_aa_metal_dist_ok = 1

        return True

class Search_vdM:
    '''
    The function to search comb is based on nearest neighbor function of sklearn.
    '''
    def __init__(self, target_pdb, workdir, querys, id_cluster_dict, cluster_centroid_dict, all_metal_query, cluster_centroid_origin_dict = None, num_iters = [3], 
    rmsd = 0.25, win_filtered = None, contact_querys = None, secondshell_querys = None, 
    validateOriginStruct = False, filter_abple = False, filter_phipsi = False, filter_phipsi_val = 20, 
    after_search_filter = False, pair_angle_range = None, pair_aa_aa_dist_range = None, pair_metal_aa_dist_range = None,
    parallel = False):

        if workdir:
            _workdir = os.path.realpath(workdir)
            if not os.path.exists(_workdir):
                os.mkdir(_workdir)
        else:
            _workdir = os.getcwd() + '/output_' + datetime.now().strftime('%Y-%m-%d-%H-%M-%S')          
            os.mkdir(_workdir)
        self.workdir = _workdir + '/'

        self.target = pr.parsePDB(target_pdb)

        self.num_iters = num_iters

        self.dist_array, self.id_array, self.dists = utils.get_contact_map(self.target, win_filtered)

        self.target_abple, self.phipsi = utils.seq_get_ABPLE(self.target)

        self.rmsd = rmsd

        self.win_filtered = win_filtered


        #neighbor in search filter--------------
        self.validateOriginStruct = validateOriginStruct

        self.filter_abple = filter_abple

        self.filter_phipsi = filter_phipsi

        self.filter_phipsi_val = filter_phipsi_val

        #neighbor after search filter-----------
        self.after_search_filter = after_search_filter
        self.pair_angle_range = pair_angle_range
        self.pair_aa_aa_dist_range =  pair_aa_aa_dist_range
        self.pair_metal_aa_dist_range = pair_metal_aa_dist_range

        #neighbor parallel mechanism------------- 
        self.parallel = parallel

        #neighbor searching strategy------------- 
        self.querys = querys #[query]
        self.all_metal_query = all_metal_query #The query with all_metal_coord_ag
        self.id_cluster_dict = id_cluster_dict # {metal_id 1234: (HIS, 0)} 
        self.cluster_centroid_dict = cluster_centroid_dict #{(HIS, 0): centroid}
        self.cluster_centroid_origin_dict = cluster_centroid_origin_dict #{(HIS, 0): centroid} Similar to cluster_centroid_dict, except keep the origin cluster coords.


        self.neighbor_query_dict = dict() # {93: [the only centroid query with all metal coords]}
        self.neighbor_pair_dict = dict() # {(33, 37): [xs-33 near 37 coords]}
        self.neighbor_comb_dict = dict() 
        # { (wins, ids), (comb, combinfo)}
        # {((0, 1, 2, 3)(0, 0, 0, 0)): {(0:[1, 3, 4], 1:[2, 3, 4], 2: [2, 6, 7], 3:[1, 2, 3]), combinfo}} Please check neighbor_win2comb()

        #contact atoms for geometry vdM score----
        self.contact_querys = contact_querys

        #secondshell-----------------------------
        self.secondshell_querys = secondshell_querys
        
        #For developing output purpose-----------
        self.log = ''
        #end-------------------------------------


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

        if self.parallel:
            self.neighbor_search_wins_pool()
        else:
            self.neighbor_search_wins()

        self.neighbor_write_represents()

        self.neighbor_write_summary()

        return


    def neighbor_generate_query_dict(self):
        '''
        return self.neighbor_query_dict = dict() # {93: [the only centroid query with all metal coords]}

        '''
        print('neighbor_generate_query_dict')
        wins = []
        if self.win_filtered:
            wins.extend([w for w in self.win_filtered])
        else:
            t = self.target.select('name CA').getResindices()
            wins.extend(([w for w in t]))

        for w in wins: 
            if self.validateOriginStruct:
                #here only filter aa, note the overlap that HIS still supperimpose to GLU.
                if self.target.select('resindex ' + str(w) + ' name CA').getResnames()[0] not in ['HIS', 'GLU', 'ASP', 'CYS']:
                    continue

            _query = self.all_metal_query.copy()
            x = self.supperimpose_target_bb(_query, w, align_sel='name N CA C')
            if x:
                self.neighbor_query_dict[w] = _query
        return

        
    def supperimpose_target_bb(self, _query, win, align_sel='name N CA C'):
        '''
        Copy the all_metal_query to a new object.
        Transform the copied all_metal_query to the target win. 
        '''
        #_query = self.all_metal_query.copy()

        target_sel = 'resindex ' + str(win) + ' and ' + align_sel
        query_sel = 'resindex ' + str(_query.contact_resind) + ' and '+ align_sel

        if len(_query.query.select(query_sel)) != len(self.target.select(target_sel)):
            print('supperimpose_target_bb not happening')
            return False
        
        transform = pr.calcTransformation(_query.query.select(query_sel), self.target.select(target_sel))
        transform.apply(_query.query)
        if _query.hull_ag:
            transform.apply(_query.hull_ag)       

        return True

    def supperimpose_centroid(self, _query, centroid, align_sel='heavy'):
        '''
        supperimpose_centroid
        '''

        if len(_query.query.select(align_sel)) != len(centroid.query.select(align_sel)):
            print('supperimpose_target_bb not happening')
            return False
        
        transform = pr.calcTransformation(_query.query.select(align_sel), centroid.query.select(align_sel))
        transform.apply(_query.query)

        return True


    def neighbor_generate_pair_dict(self):
        '''
        To generate pair win dict, we will first filter out the impossible win pair by distance.
        Apply utils.check_pair_distance_satisfy()
        '''

        print('neighbor_generate_pair_dict')
        
        wins = sorted(list(self.neighbor_query_dict.keys()))
        print(wins)
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
                # print('-------------------')
                # print(wx)
                # print(wy)
                x_in_y, x_has_y = self.calc_pairwise_neighbor(n_x, n_y, self.rmsd)
                y_in_x, y_has_x = self.calc_pairwise_neighbor(n_y, n_x, self.rmsd)

                # print(x_has_y)
                # print(y_has_x)

                x_in_y, x_has_y = self.neighbor_filter_new(wx, wy, x_in_y)
                y_in_x, y_has_x = self.neighbor_filter_new(wy, wx, y_in_x)

                # print(x_has_y)
                # print(y_has_x)

                if x_has_y and y_has_x:
                    self.neighbor_pair_dict[(wx, wy)] = x_in_y
                    self.neighbor_pair_dict[(wy, wx)] = y_in_x
        return


    def calc_pairwise_neighbor(self, n_x, n_y, rmsd = 0.25):
        '''
        Use sklean NearestNeighbors
        '''

        neigh_y = NearestNeighbors(radius= rmsd) 
        neigh_y.fit(n_y)

        x_in_y = neigh_y.radius_neighbors(n_x)
        x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])

        return x_in_y[1], x_has_y

    
    def neighbor_filter_new(self, wx, wy, x_in_y):
        '''
        filter impossible pairwise connect.
        1. validateOriginStruct
        2. ABPLE filter. 
        3. phipsi angle filter. (future)
        '''

        x_in_y_mask = [[True for j in range(len(x_in_y[i]))] for i in range(len(x_in_y))]

        for i in range(len(x_in_y)):
            if len(x_in_y[i]) <= 0:
                continue

            if self.validateOriginStruct:
                resx = self.target.select('name CA and resindex ' + str(wx)).getResnames()[0]
                resy = self.target.select('name CA and resindex ' + str(wy)).getResnames()[0]

                if not self.id_cluster_dict[i][0] == resx:
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))] 
                    continue

            if self.filter_abple:
                apx = self.target_abple[wx]
                apy = self.target_abple[wy]
       
                if not self.querys[i].abple == apx:
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))] 
                    continue

            
            if self.filter_phipsi:
                phix, psix = self.phipsi[wx]
                phiy, psiy = self.phipsi[wy]

                if (not utils.filter_phipsi(phix, self.querys[i].phi, self.filter_phipsi_val)) or not (utils.filter_phipsi(psix, self.querys[i].psi, self.filter_phipsi_val)):
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))]
                    continue

                
            for j in range(len(x_in_y[i])):
                j_ind = x_in_y[i][j]

                if self.validateOriginStruct:
                    resy = self.target.select('name CA and resindex ' + str(wy)).getResnames()[0]
                    if not self.id_cluster_dict[j_ind][0] == resy:
                        x_in_y_mask[i][j] = False
                        continue

                if self.filter_abple:
                    apy = self.target_abple[wy]
                    if not self.querys[j_ind].abple == apy:
                        x_in_y_mask[i][j] = False
                        continue

                if self.filter_phipsi:
                    phiy, psiy = self.phipsi[wy]
                    if (not utils.filter_phipsi(phiy, self.querys[j_ind].phi, self.filter_phipsi_val)) or (not utils.filter_phipsi(psiy, self.querys[j_ind].psi, self.filter_phipsi_val)) :
                        x_in_y_mask[i][j] = False
                        continue
        
        x_in_y_filter = [[x_in_y[i][j] for j in range(len(x_in_y[i])) if x_in_y_mask[i][j]] for i in range(len(x_in_y))]

        x_has_y = any([True if len(a) >0 else False for a in x_in_y_filter])

        return x_in_y_filter, x_has_y


    def neighbor_filter(self, wx, wy, x_in_y):
        '''
        filter impossible pairwise connect.
        1. validateOriginStruct
        2. ABPLE filter. 
        3. phipsi angle filter. (future)
        '''

        x_in_y_filter = [[j for j in x_in_y[i]] for i in range(len(x_in_y))]
        if self.validateOriginStruct:
            resx = self.target.select('name CA and resindex ' + str(wx)).getResnames()[0]
            resy = self.target.select('name CA and resindex ' + str(wy)).getResnames()[0]
            
            for i in range(len(x_in_y)):
                if len(x_in_y[i]) <= 0:
                    continue
                if not self.id_cluster_dict[i][0] == resx:
                    x_in_y_filter[i].clear()
                    continue
                
                for j in range(len(x_in_y[i])):
                    j_ind = x_in_y[i][j]
                    if not self.id_cluster_dict[j_ind][0] == resy:
                        x_in_y_filter[i].remove(j_ind)


        if self.filter_abple:
            apx = self.target_abple[wx]
            apy = self.target_abple[wy]

            for i in range(len(x_in_y)):
                if len(x_in_y[i]) <= 0:
                    continue

                if not self.querys[i].abple == apx:
                    x_in_y_filter[i].clear()
                    continue

                for j in range(len(x_in_y[i])):
                    j_ind = x_in_y[i][j]
                    if j_ind not in x_in_y_filter[i]:
                        continue
                    if not self.querys[j_ind].abple == apy:
                        x_in_y_filter[i].remove(j_ind)


        if self.filter_phipsi:
            phix, psix = self.phipsi[wx]
            phiy, psiy = self.phipsi[wy]

            for i in range(len(x_in_y)):
                if len(x_in_y[i]) <= 0:
                    continue                

                if (not utils.filter_phipsi(phix, self.querys[i].phi, self.filter_phipsi_val)) or not (utils.filter_phipsi(psix, self.querys[i].psi, self.filter_phipsi_val)):
                    x_in_y_filter[i].clear()
                    continue

                for j in range(len(x_in_y[i])):
                    j_ind = x_in_y[i][j]
                    if j_ind not in x_in_y_filter[i]:
                        continue

                    if (not utils.filter_phipsi(phiy, self.querys[j_ind].phi, self.filter_phipsi_val)) or (not utils.filter_phipsi(psiy, self.querys[j_ind].psi, self.filter_phipsi_val)) :
                        x_in_y_filter[i].remove(j_ind)
                        continue
        

        ### Debug print.
        # for i in range(len(x_in_y)):
        #     if len(x_in_y[i]) <= 0:
        #         continue
        #     if self.id_cluster_dict[i][0] == 'HIS' and 'CYS' in [self.id_cluster_dict[j_ind][0] for j_ind in x_in_y[i]] and 'HIS' in [self.id_cluster_dict[j_ind][0] for j_ind in x_in_y[i]]:
        #         #print(x_in_y[i])          
        #         print([self.id_cluster_dict[j_ind][0] for j_ind in x_in_y[i]])
        #         #print(resx)
        #         #print(resy)
        #         print('x: ' + self.id_cluster_dict[i][0])
        #         #print(x_in_y_filter[i])
        #         print([self.id_cluster_dict[j_ind][0] for j_ind in x_in_y_filter[i]])
        #         print('---------')
                  
        x_has_y = any([True if len(a) >0 else False for a in x_in_y_filter])
        return x_in_y_filter, x_has_y


    def neighbor_get_win_combs(self):
        '''
        The combinations of positions are extracted from all possible positions with pairs.
        '''
        print('pair_dict_keys: {}'.format(self.neighbor_pair_dict.keys()))
        wins = sorted(set(self.neighbor_query_dict.keys()))
        print('All wins with overlap {}'.format(wins))

        win_ind_dict = {}
        for i in range(len(wins)):
            win_ind_dict[wins[i]] = i

        pair_set = set()

        for x, y in self.neighbor_pair_dict.keys():
            i = win_ind_dict[x]
            j = win_ind_dict[y]
            pair_set.add((i, j))

        paths = []
        for n in self.num_iters: 
            paths.extend(utils.combination_calc(list(range(len(wins))), pair_set, n))

        win_combs = []
        for path in paths:
            win_comb = [wins[p] for p in path]
            win_combs.append(win_comb)
        return win_combs

    
    def neighbor_search_wins(self):
        '''
        The combinations of positions are extracted from all possible positions with pairs.
        '''
        win_combs = self.neighbor_get_win_combs()

        for win_comb in win_combs:
            print(win_comb)      
            comb_dict = self.neighbor_construct_comb(win_comb)
            self.neighbor_write_win(comb_dict)
            self.neighbor_comb_dict.update(comb_dict)
        return


    def neighbor_search_wins_pool(self):
        '''
        multithread.
        '''
        win_combs = self.neighbor_get_win_combs()
        print('multithread search win_combs: {}'.format(len(win_combs)))
        if len(win_combs) <= 0:
            return
        num_cores = int(mp.cpu_count() - 1)
        print('pool: {}'.format(num_cores))
        # pool = mp.Pool(num_cores)
        # results = [pool.apply_async(self.neighbor_construct_comb, args=win_comb) for win_comb in win_combs]
        # results = [p.get() for p in results]
        pool = ThreadPool(num_cores)
        results = pool.map(self.neighbor_construct_comb, win_combs)
        pool.close()
        pool.join()
        for r in results:
            self.neighbor_write_win(r)
            self.neighbor_comb_dict.update(r)
        return


    def neighbor_construct_comb(self, win_comb):
        '''
        win_comb: [0, 1, 2, 3]

        cluster_dict: {(0, 0, 0, 0): {0:[1, 3, 4], 1:[2, 3, 4], 2: [2, 6, 7], 3:[1, 2, 3]}}
        The key is win, the value is list of index, each represent one metal coord exist in all other wins' metal.
        
        self.neighbor_comb_dict
        # { (wins, ids), (comb, combinfo)}
        # {((0, 1, 2, 3)(0, 0, 0, 0)): {(0:[1, 3, 4], 1:[2, 3, 4], 2: [2, 6, 7], 3:[1, 2, 3]), combinfo}}


        As we calculate all metals at the same time. 
        We need to extract vdm representations.
        '''
        print('neighbor_construct comb: {}'.format(win_comb))

        comb_dict = dict()

        for x, y in itertools.permutations(win_comb, 2):
            if (x, y) not in self.neighbor_pair_dict.keys():
                return None
        
        graph = Graph(win_comb, len(self.querys))

        graph.calc_pair_connectivity(self.neighbor_pair_dict)

        graph.get_paths()

        print('graph.paths len {}'.format(len(graph.all_paths)))

        #TO DO: Here is a temp method to solve extream solutions. Mostly happened in 4 CYS binding cores.
        if len(graph.all_paths) > 100000:
            print('Too many paths to be considered so far.')
            graph.all_paths = graph.all_paths[0:100001]

        clu_dict = {}

        # path represent the id of each metal vdM.
        for path in graph.all_paths:
            
            clu_key = tuple([self.id_cluster_dict[p] for p in path])

            ### Debug
            # if not (path[0] == 7208 and path[1] == 7864 and path[2] == 8539):
            #     continue

            # if clu_key == (('HIS', 7), ('HIS', 17), ('HIS', 49)):
            #     print(path)

            
            if clu_key in clu_dict:
                for i in range(len(win_comb)):
                    clu_dict[clu_key][win_comb[i]].add(path[i])
            else:
                clu_dict[clu_key] = {}
                for i in range(len(win_comb)):
                    clu_dict[clu_key][win_comb[i]] = set()
                    clu_dict[clu_key][win_comb[i]].add(path[i])

        if len(clu_dict) > 0:
            for clu_key in clu_dict.keys():
                comb = clu_dict[clu_key]
                comb_dict[(tuple(win_comb), clu_key)] = (comb, CombInfo())

        self.neighbor_extract_query(comb_dict)
        self.neighbor_calc_comb_score(comb_dict)
        self.neighbor_calc_geometry(comb_dict)

        return comb_dict


    def neighbor_extract_query(self, comb_dict):
        '''
        for each (comb, combinfo), we extract candidate querys and add it in combinfo.
        '''
        print('neighbor_extract_query')
        for wins, clu_key in comb_dict.keys():
            for i in range(len(wins)):
                win = wins[i]
                comb_dict[(wins, clu_key)][1].query_dict[win] = []

                clu = clu_key[i]
                centroid = self.cluster_centroid_dict[clu].copy()
                self.supperimpose_target_bb(centroid, win)
                comb_dict[(wins, clu_key)][1].centroid_dict[win] = centroid

                for id in comb_dict[(wins, clu_key)][0][win]:
                    _query = self.querys[id].copy()
                    
                    self.supperimpose_target_bb(_query, win)
                    #self.supperimpose_centroid(_query, centroid, align_sel='heavy') # Test different supperimposition effect of search.
                    comb_dict[(wins, clu_key)][1].query_dict[win].append(_query)

        return


    def neighbor_calc_comb_score(self, comb_dict):
        '''
        The summed vdM score could not reflect the designability.
        Here is a new score method with weight added.

        In the future, the vdM score can be pre-calcluated as COMBS did.
        '''
        print('neighbor_calc_comb_score')

        #For multiScore.
        aa_num_dict = {}
        for key in self.id_cluster_dict.keys():
            aa = self.id_cluster_dict[key][0]
            if aa not in aa_num_dict.keys():
                aa_num_dict[aa] = 0
            else:
                aa_num_dict[aa] += 1
        print(aa_num_dict)

        for wins, clu_key in comb_dict.keys():       
            comb = comb_dict[(wins, clu_key)][0] 

            # From here we will calculate the score for each comb. 
            totals = [len(comb[wins[i]]) for i in range(len(wins))]
            #total_clu = sum([self.cluster_centroid_dict[clu_key[i]].total_clu for i in range(len(wins))])
            scores = [self.cluster_centroid_dict[clu_key[i]].score for i in range(len(wins))]
            comb_dict[(wins, clu_key)][1].totals = totals
            comb_dict[(wins, clu_key)][1].scores = scores

            #Calc fraction score
            fracScore = sum([self.cluster_centroid_dict[clu_key[i]].clu_num for i in range(len(wins))])/sum([self.cluster_centroid_dict[clu_key[i]].clu_total_num for i in range(len(wins))])
            weight = 0
            for i in range(len(wins)):
                c = self.cluster_centroid_dict[clu_key[i]]
                weight += sum(totals)/c.clu_total_num
            comb_dict[(wins, clu_key)][1].fracScore = fracScore*weight*1000

            #Calc multiScore (By Bill: -ln(Na/SumNa * Nb/SumNb * Nc/SumNc))
            aa_numbs = [aa_num_dict[c[0]] for c in clu_key]
            comb_dict[(wins, clu_key)][1].multiScore = -log(np.prod(totals)/np.prod(aa_numbs))

        return
    

    def neighbor_calc_geometry(self, comb_dict):
        '''
        The overlap has several members for each query. 
        The geometry is a centroid contact atom of each query's candidates.
        Check CombInfo.calc_geometry()
        '''
        print('neighbor_calc_geometry')
        for wins, clu_key in comb_dict.keys():
            comb_dict[(wins, clu_key)][1].calc_geometry()         
        return


    def neighbor_write_represents(self):
        print('neighbor_write_represents')
        for key in self.neighbor_comb_dict.keys():  
            info = self.neighbor_comb_dict[key][1]
            if self.after_search_filter:
                if not info.after_search_condition_satisfied(self.pair_angle_range, self.pair_aa_aa_dist_range, self.pair_metal_aa_dist_range):
                    continue

            outdir = self.workdir + 'represents/'
            if not os.path.exists(outdir):
                os.mkdir(outdir)        

            tag = 'win_' + '-'.join([str(k) for k in key[0]]) + '_clu_' + '-'.join(k[0] + '-' + str(k[1]) for k in key[1]) 
            for w in key[0]:
                c = self.neighbor_comb_dict[key][1].query_dict[w][0]
                pdb_path = outdir + tag + '_w_' + str(w) + '_apble_' + c.abple + '_' + c.query.getTitle() + '.pdb'
                pr.writePDB(pdb_path, c.query)                  

        return

    def neighbor_write_win(self, comb_dict):
        '''
        Write output.
        Too many output, May need optimization. 
        '''
        print('neighbor_write')
        for key in comb_dict.keys():  
            info = comb_dict[key][1]
            if self.after_search_filter:
                if not info.after_search_condition_satisfied(self.pair_angle_range, self.pair_aa_aa_dist_range, self.pair_metal_aa_dist_range):
                    continue
            outpath = 'win_' + '-'.join([str(k) for k in key[0]]) + '/'
            outdir = self.workdir + outpath
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            tag = 'win_' + '-'.join([str(k) for k in key[0]]) + '_clu_' + '-'.join(k[0] + '-' + str(k[1]) for k in key[1]) 

            # Write geometry       
            pr.writePDB(outdir + tag +'_geometry.pdb', comb_dict[key][1].geometry) 
            
            #Write 
            all_overlap_metal_coords = []
            for w in key[0]:
                candidates = comb_dict[key][1].query_dict[w]
                metal_coords = []

                max_out = 100 #To reduce number of output.
                for c in candidates:
                    if max_out >0:
                        pdb_path = outdir + tag + '_w_' + str(w) + '_apble_' + c.abple + '_' + c.query.getTitle() + '.pdb'
                        pr.writePDB(pdb_path, c.query)                  
                        max_out-=1
                        metal_coords.append(c.get_metal_coord())
                        all_overlap_metal_coords.extend(metal_coords)
                    else:
                        break

                hull.write2pymol(metal_coords, outdir, tag + '_w_' + str(w) +'_overlap_points.pdb')

            # if len(all_overlap_metal_coords) > 3:
            #     vol = scipy.spatial.ConvexHull(all_overlap_metal_coords)
            #     volume = vol.volume
            # else:
            #     volume = 0
            # self.neighbor_comb_dict[key][1].volume = volume
            # self.neighbor_comb_dict[key][1].volPerMetal = volume/len(all_overlap_metal_coords)
            # hdist = cdist(all_overlap_metal_coords, all_overlap_metal_coords, metric='euclidean')
            # self.neighbor_comb_dict[key][1].diameter = hdist.argmax()
            
            #Write Centroid and all metal coords in the cluster
            for w in key[0]:
                centroid = comb_dict[key][1].centroid_dict[w]
                pdb_path = outdir + tag + '_centroid_' + centroid.query.getTitle() + '.pdb'
                pr.writePDB(pdb_path, centroid.query)
                clu_allmetal_coords = centroid.get_hull_points()
                hull.write2pymol(clu_allmetal_coords, outdir, tag + '_w_' + str(w) +'_points.pdb')  
        return   

    def neighbor_write(self):
        '''
        Write output all comb_dict results together.
        '''
        print('neighbor_write')
        for key in self.neighbor_comb_dict.keys():  
            comb_dict = self.neighbor_comb_dict[key]
            self.neighbor_write_win(comb_dict)
        return    


    def neighbor_write_summary(self, eval = False):
        '''
        Write a tab dilimited file.
        '''
        print('neighbor_write_summary')
        with open(self.workdir + '_summary.tsv', 'w') as f:
            f.write('Wins\tClusterIDs\tproteinABPLEs\tCentroidABPLEs\tproteinPhiPsi\tCentroidPhiPsi\tvolume\tvol2metal\tdiameter\tTotalVdMScore\tFracScore\tMultiScore\taa_aa_dists\tmetal_aa_dists\tPair_angles\toverlap#\toverlaps#\tvdm_scores\ttotal_clu#\tclu_nums')
            if eval:
                f.write('\teval_min_rmsd\teval_min_vdMs\teval_phi\teval_psi\teval_abple\teval_is_origin')
            f.write('\n')
            for key in self.neighbor_comb_dict.keys(): 
                info = self.neighbor_comb_dict[key][1]
                if self.after_search_filter:
                    if not info.after_search_condition_satisfied(self.pair_angle_range, self.pair_aa_aa_dist_range, self.pair_metal_aa_dist_range):
                        continue

                #centroids = [c.query.getTitle() for c in info.centroid_dict.values()]
                vdm_scores = [c for c in info.scores]
                overlaps = [c for c in info.totals]
                clu_nums = [c.clu_num for c in info.centroid_dict.values()]
                
                f.write('_'.join([str(x) for x in key[0]]) + '\t')
                f.write('_'.join([x[0] + '-' + str(x[1]) for x in key[1]]) + '\t')

                f.write('_'.join([self.target_abple[x] for x in key[0]]) + '\t')
                f.write('_'.join([c.abple for c in info.centroid_dict.values()]) + '\t')
                f.write('_'.join([str((round(self.phipsi[x][0],2), round(self.phipsi[x][1],2))) for x in key[0]]) + '\t')
                f.write('_'.join([str((round(c.phi, 2), round(c.psi, 2))) for c in info.centroid_dict.values()]) + '\t')

                f.write(str(round(info.volume,2))+ '\t')
                f.write(str(round(info.volPerMetal,2))+ '\t')
                f.write(str(round(info.diameter,2))+ '\t')

                f.write(str(round(sum(info.scores), 2)) + '\t')
                f.write(str(round(info.fracScore, 2)) + '\t')
                f.write(str(round(info.multiScore, 2)) + '\t')

                
                f.write('||'.join([str(round(d, 2)) for d in info.aa_aa_pair])  + '\t')
                f.write('||'.join([str(round(d, 2)) for d in info.metal_aa_pair])  + '\t')
                f.write('||'.join([str(round(a, 2)) for a in info.angle_pair])  + '\t')

                f.write(str(sum(overlaps)) + '\t')
                f.write('||'.join([str(s) for s in overlaps]) + '\t')

                f.write('||'.join([str(round(s, 2)) for s in vdm_scores]) + '\t')
                f.write(str(sum(clu_nums)) + '\t')
                f.write('||'.join([str(c) for c in clu_nums]) + '\t')

                if eval:
                    f.write('||'.join([str(round(m, 2)) for m in info.eval_mins]) + '\t')
                    f.write('||'.join([v.query.getTitle() for v in info.eval_min_vdMs]) + '\t')
                    f.write('||'.join([str(round(v.phi,2)) for v in info.eval_min_vdMs]) + '\t')
                    f.write('||'.join([str(round(v.psi,2)) for v in info.eval_min_vdMs]) + '\t')
                    f.write('||'.join([v.abple for v in info.eval_min_vdMs]) + '\t')
                    f.write(str(info.eval_is_origin) + '\t')

                f.write('\n')


        if len(self.log) > 0:
            with open(self.workdir + '_log.txt', 'w') as f:
                f.write(self.log)
            
        return 


    def eval_search(self):
        '''
        Given a metal binding protein, extract the binding core. For each contact aa, find the original or best vdM as target. 
        For all the 'target' vdMs, use nearest neighbor search to check overlap and combinfo.
        '''
        wins, combs = self.eval_get_comb()

        uni_wins = set()
        [uni_wins.add(w) for win in wins for w in win]
        self.win_filtered = sorted(uni_wins)

        self.neighbor_generate_query_dict()

        #Extract closest vdMs and the cluster infomation. 
        self.eval_extract_closest_vdMs(wins, combs)

        #The normal search process focus on the current wins.
        self.neighbor_generate_pair_dict()

        for win_comb in wins:
            print(win_comb)      
            comb_dict = self.neighbor_construct_comb(win_comb)
            if not comb_dict: continue
            self.neighbor_write_win()
            self.neighbor_comb_dict.update(comb_dict)     

        #self.neighbor_write()

        #Evaluate search result.
        self.eval_search_results(wins, combs)

        self.neighbor_write_summary(eval=True)

        return 
                

    def eval_get_comb(self):
        '''
        Extract vdMs from the target protein.
        '''
        aas = ['HIS', 'GLU', 'ASP', 'CYS']

        wins = []
        combs = []
        _cores = ligand_database.get_metal_core_seq(self.target, quco.metal_sel, extend = 4)
        cores = [core.Core(c[1]) for c in _cores]

        if len(_cores) <= 0:
            self.log += 'No core exist (old vdM).\n'

        for c in cores:
            for aa in c.contact_aas.getResnames():
                if aa not in aas:
                    self.log += 'core contain other aa.\n'
                    continue
            comb = [c._generate_AA_phipsi_Metal(w) for w in c.contact_aa_resinds]
            if None in comb:
                self.log += 'core extract None.\n'
                continue
            wins.append(c.contact_aa_resinds)
            combs.append(comb)
        return wins, combs


    def _eval_extract_closest_vdM(self, w, v):
        best_v = None
        best_id = -1
        min_rmsd = 10

        if not v:
            print('The vdM at {} from the protein bb is not successfully extracted.'.format(w))
            return best_v, best_id, min_rmsd

        coord = [v.select(quco.metal_sel)[0].getCoords()]

        ns = self.neighbor_query_dict[w].get_hull_points()

        x_in_y, x_has_y = self.calc_pairwise_neighbor(coord, ns, 1.5)


        for ind in x_in_y[0]:

            if len(self.querys[ind].query.select('heavy')) != len(v.select('heavy')):
                #print('supperimpose_target_bb not happening')
                continue

            test_v = self.querys[ind].copy()
    
            transform = pr.calcTransformation(test_v.query.select('heavy'), v.select('heavy'))
            transform.apply(test_v.query)
            rmsd = pr.calcRMSD(v.select('heavy'), test_v.query.select('heavy'))

            if rmsd < min_rmsd:
                best_v = test_v
                best_id = ind
                min_rmsd = rmsd

        return best_v, best_id, min_rmsd


    def write_closest_vdM_clu_points(self, best_v, w, evaldir, tag):
        clu_key = best_v.get_cluster_key()

        origin_centroid = self.cluster_centroid_origin_dict[clu_key].copy()
        self.supperimpose_target_bb(origin_centroid, w)

        #pr.writePDB(evaldir + tag + 'origin_' + origin_centroid.query.getTitle(), origin_centroid.query)
        clu_origin_allmetal_coords = origin_centroid.get_hull_points()
        hull.write2pymol(clu_origin_allmetal_coords, evaldir, tag + '_origin_w_' + str(w) +'_points.pdb') 


        centroid = self.cluster_centroid_dict[clu_key].copy()
        self.supperimpose_target_bb(centroid, w)

        pr.writePDB(evaldir + tag + centroid.query.getTitle(), centroid.query)
        clu_allmetal_coords = centroid.get_hull_points()
        hull.write2pymol(clu_allmetal_coords, evaldir, tag + '_w_' + str(w) +'_points.pdb') 

        origin_best_v = best_v.copy()
        transform = pr.calcTransformation(origin_best_v.query.select('heavy'), centroid.query.select('heavy'))
        transform.apply(origin_best_v.query)
        pr.writePDB(evaldir + tag + 'origin_' + origin_best_v.query.getTitle(), origin_best_v.query)

        return


    def eval_extract_closest_vdMs(self, wins, combs):
        '''
        First supperimpose the all_metal_query. 
        Then use nearest neighbor to get candidates by calculate the metal distance with dist 0.25.
        Then superimpose and obtain the one with min rmsd. 
        '''

        for i in range(len(wins)):
            evaldir = self.workdir + 'closest_win_' + '_'.join([str(w) for w in wins[i]])
            os.makedirs(evaldir, exist_ok=True)

            for j in range(len(wins[i])):
                w = wins[i][j]
                v = combs[i][j]
                best_v, best_id, min_rmsd = self._eval_extract_closest_vdM(w, v)
                
                if not v:
                    continue
                pr.writePDB(evaldir + '/win_' + str(w) + '_' + v.getTitle(), v)
                if best_v:
                    clu_id = self.id_cluster_dict[best_id]
                    tag = '/win_' + str(w) + '_clu_' + '_'.join([str(ci) for ci in clu_id]) + '_rmsd_' + str(round(min_rmsd, 3)) + '_'                   
                    print(tag + ' : best_id ' + str(best_id))
                    pr.writePDB(evaldir + tag + best_v.query.getTitle(), best_v.query)
                    self.write_closest_vdM_clu_points(best_v, w, evaldir, tag)

        return

    
    def eval_extract_comb_closest_vdMs(self, w, v, combinfo):
        '''
        '''
        best_v = None
        min_rmsd = 10

        if not v:
            return best_v, min_rmsd

        for test_v in combinfo.query_dict[w]:
            if len(test_v.query.select('heavy')) != len(v.select('heavy')):
                #print('supperimpose_target_bb not happening')
                continue
            
            transform = pr.calcTransformation(test_v.query.select('heavy'), v.select('heavy'))
            transform.apply(test_v.query)
            rmsd = pr.calcRMSD(v.select('heavy'), test_v.query.select('heavy'))

            if rmsd < min_rmsd:
                best_v = test_v
                min_rmsd = rmsd

        return best_v, min_rmsd


    def eval_search_results(self, wins, combs):
        '''
        After the nearest neighbor search, find the closest one from each neighbor_comb_dict.
        '''
        for i in range(len(wins)):
            evaldir = self.workdir + 'eval_win_' + '_'.join([str(w) for w in wins[i]])
            os.makedirs(evaldir, exist_ok=True)
            
            for key in self.neighbor_comb_dict.keys():
                if not key[0] == tuple([w for w in wins[i]]): 
                    continue
                
                clu_id = key[1]
                combinfo = self.neighbor_comb_dict[key][1]

                min_rmsds = []
                best_vs = []
                for j in range(len(wins[i])):
                    w = wins[i][j]
                    v = combs[i][j]                    

                    best_v, min_rmsd = self.eval_extract_comb_closest_vdMs(w, v, combinfo)
                    best_vs.append(best_v)
                    min_rmsds.append(min_rmsd)

                    if not v:
                        continue
                    pr.writePDB(evaldir + '/win_' + str(w) + '_' + v.getTitle(), v)
                    if best_v:
                        tag = '/clu_' + '_'.join([str(ci) for cid in clu_id for ci in cid]) + '_rmsd_' + str(round(min_rmsd, 3)) + '_win_' + str(w) + '_clu_' + '_'.join([str(ci) for ci in clu_id[j]]) + '_'                   
                        pr.writePDB(evaldir + tag + best_v.query.getTitle(), best_v.query)   
                self.neighbor_comb_dict[key][1].eval_mins = min_rmsds
                self.neighbor_comb_dict[key][1].eval_min_vdMs = best_vs
                if all([m < 0.05 for m in min_rmsds]):
                    self.neighbor_comb_dict[key][1].eval_is_origin = True
        return