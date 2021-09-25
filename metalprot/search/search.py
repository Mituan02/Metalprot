import contextlib
from math import comb, log
import os
from typing import Dict
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

from ..basic import quco
from ..basic import hull
from ..basic import utils
from ..basic.filter import Search_filter
from ..database import core
from ..database import database_extract
from .graph import Graph
from .comb_info import CombInfo

from sklearn.neighbors import NearestNeighbors
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool


def supperimpose_target_bb(_target, _query, win, align_sel='name N CA C'):
    '''
    Copy the all_metal_query to a new object.
    Transform the copied all_metal_query to the target win. 
    '''

    target_sel = 'resindex ' + str(win) + ' and ' + align_sel
    query_sel = 'resindex ' + str(_query.contact_resind) + ' and '+ align_sel
    # print('a-----------------------------------')
    # print(target_sel)
    # print(query_sel)
    # print(len(_query.query))
    # print(len(_target))
    q = _query.query.select(query_sel)
    t = _target.select(target_sel)
    if len(q) != len(t):
        print('supperimpose-target-bb not happening')
        return False
    # print('b-----------------------------------')
    # print(target_sel)
    # print(query_sel)
    # print(len(q))
    # print(len(t))
    transform = pr.calcTransformation(q, t)
    transform.apply(_query.query)
    if _query.hull_ag:
        transform.apply(_query.hull_ag)  

    return True


def supperimpose_centroid(_query, centroid, align_sel='heavy'):
    '''
    supperimpose_centroid
    '''

    if len(_query.query.select(align_sel)) != len(centroid.query.select(align_sel)):
        print('supperimpose-centroid not happening')
        return False
    
    transform = pr.calcTransformation(_query.query.select(align_sel), centroid.query.select(align_sel))
    transform.apply(_query.query)

    return True

class Search_vdM:
    '''
    The function to search comb is based on nearest neighbor function of sklearn.
    '''
    def __init__(self, target_pdb, workdir, querys, id_cluster_dict, cluster_centroid_dict, all_metal_query, cluster_centroid_origin_dict = None, num_iters = [3], 
    rmsd = 0.25, win_filtered = None, contact_querys = None, secondshell_querys = None, 
    validateOriginStruct = False, search_filter = None, parallel = False, selfcenter_rmsd = 0.45):

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

        self.search_filter = search_filter

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

        #For_scoring-----------------------------
        self.aa_num_dict = None
        
        #For developing output purpose-----------
        self.log = ''

        #For Search_selfcenter-------------------
        self.selfcenter_rmsd = selfcenter_rmsd
        self.best_aa_comb_dict = {} # To store&Write the best comb for each combinations of wins. 
        #----------------------------------------
        self.setup()
        #end-------------------------------------


    def setup(self):

        #For multiScore.
        aa_num_dict = {}
        for key in self.id_cluster_dict.keys():
            aa = self.id_cluster_dict[key][0]
            if aa not in aa_num_dict.keys():
                aa_num_dict[aa] = 0
            else:
                aa_num_dict[aa] += 1
        print(aa_num_dict)

        self.aa_num_dict = aa_num_dict
        return


    def para2string(self):
        parameters = "Search parameters: \n"
        parameters += 'target: ' + self.target.getTitle() + ' \n'
        parameters += 'num_iters: ' + str(self.num_iters) + ' \n'
        parameters += 'rmsd: ' + str(self.rmsd) + ' \n'
        parameters += 'win_filtered: ' + str(self.win_filtered) + ' \n'
        parameters += 'validateOriginStruct: ' + str(self.validateOriginStruct) + ' \n'
        parameters += 'parallel: ' + str(self.parallel) + ' \n'
        parameters += 'selfcenter_rmsd: ' + str(self.selfcenter_rmsd) + ' \n'
        return parameters

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

        self.neighbor_write_summary(self.workdir, self.neighbor_comb_dict)

        self.neighbor_write_log()

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
            x = supperimpose_target_bb(self.target, _query, w, align_sel='name N CA C')
            if x:
                self.neighbor_query_dict[w] = _query
        return



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
        3. phipsi angle filter. 
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

            if self.search_filter.filter_abple:
                apx = self.target_abple[wx]
                apy = self.target_abple[wy]
       
                if not self.querys[i].abple == apx:
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))] 
                    continue

            
            if self.search_filter.filter_phipsi:
                phix, psix = self.phipsi[wx]
                phiy, psiy = self.phipsi[wy]

                if (not utils.filter_phipsi(phix, self.querys[i].phi, self.search_filter.max_phipsi_val)) or not (utils.filter_phipsi(psix, self.querys[i].psi, self.search_filter.max_phipsi_val)):
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))]
                    continue
           

            if self.search_filter.filter_vdm_score:
                if self.querys[i].score < self.search_filter.min_vdm_score:
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))]
                    continue
            

            if self.search_filter.filter_vdm_count:
                if self.querys[i].clu_num < self.search_filter.min_vdm_clu_num:
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))]
                    continue
                           

            for j in range(len(x_in_y[i])):
                j_ind = x_in_y[i][j]

                if self.validateOriginStruct:
                    resy = self.target.select('name CA and resindex ' + str(wy)).getResnames()[0]
                    if not self.id_cluster_dict[j_ind][0] == resy:
                        x_in_y_mask[i][j] = False
                        continue

                if self.search_filter.filter_abple:
                    apy = self.target_abple[wy]
                    if not self.querys[j_ind].abple == apy:
                        x_in_y_mask[i][j] = False
                        continue

                if self.search_filter.filter_phipsi:
                    phiy, psiy = self.phipsi[wy]
                    if (not utils.filter_phipsi(phiy, self.querys[j_ind].phi, self.search_filter.max_phipsi_val)) or (not utils.filter_phipsi(psiy, self.querys[j_ind].psi, self.search_filter.max_phipsi_val)) :
                        x_in_y_mask[i][j] = False
                        continue

                if self.search_filter.filter_vdm_score:
                    if self.querys[j_ind].score < self.search_filter.min_vdm_score:
                        x_in_y_mask[i][j] = False 
                        continue
            

                if self.search_filter.filter_vdm_count:
                    if self.querys[j_ind].clu_num < self.search_filter.min_vdm_clu_num:
                        x_in_y_mask[i][j] = False 
                        continue
        
        x_in_y_filter = [[x_in_y[i][j] for j in range(len(x_in_y[i])) if x_in_y_mask[i][j]] for i in range(len(x_in_y))]

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
            comb_dict = self.neighbor_run_comb(win_comb)
            if not comb_dict:
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
        results = pool.map(self.neighbor_run_comb, win_combs)
        pool.close()
        pool.join()
        for r in results: 
            if not r:
                self.neighbor_comb_dict.update(r)
        return


    def neighbor_run_comb(self, win_comb):
        #try:
        comb_dict = self.neighbor_construct_comb(win_comb)
        if len([comb_dict.keys()]) <= 0:
            return comb_dict
        _target = self.target.copy()
        self.neighbor_extract_query(_target, comb_dict)
        self.neighbor_calc_comb_score(comb_dict)
        self.neighbor_calc_geometry(comb_dict)
        self.neighbor_aftersearch_filt(_target, comb_dict)
        self.neighbor_write_win(comb_dict)
        return comb_dict
        # except:
        #     self.log += 'Error in win_comb: ' + '-'.join([str(w) for w in win_comb]) + '\n'
        #     return {}

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
        
        graph = Graph(win_comb, [len(self.querys) for i in range(len(win_comb))])

        graph.calc_pair_connectivity(self.neighbor_pair_dict)

        graph.get_paths()

        print('graph.paths len {}'.format(len(graph.all_paths)))

        #TO DO: Here is a temp method to solve extream solutions. Mostly happened in 4 CYS binding cores.
        if len(graph.all_paths) > 10000:
            print('Too many paths to be considered so far.')
            graph.all_paths = graph.all_paths[0:10001]

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
                combinfo = CombInfo() 
                combinfo.comb = clu_dict[clu_key]
                comb_dict[(tuple(win_comb), clu_key)] = combinfo

        return comb_dict


    def neighbor_extract_query(self, _target, comb_dict):
        '''
        for each (comb, combinfo), we extract candidate querys and add it in combinfo.
        '''
        print('neighbor-extract-query')
        for wins, clu_key in comb_dict.keys():
            for i in range(len(wins)):
                win = wins[i]
                comb_dict[(wins, clu_key)].query_dict[win] = []

                clu = clu_key[i]
                centroid = self.cluster_centroid_dict[clu].copy()
                
                try:
                    supperimpose_target_bb(_target, centroid, win)
                except:
                    print('x-------------------------')
                    print(wins)
                    print(clu_key)
                    print(len(centroid.query))
                    print(len(_target))
                    print(centroid.query.getTitle())
                comb_dict[(wins, clu_key)].centroid_dict[win] = centroid

                for id in comb_dict[(wins, clu_key)].comb[win]:
                    _query = self.querys[id].copy()
                    try:
                        supperimpose_target_bb(_target, _query, win)
                    except:
                        print('y-------------------------')
                        print(wins)
                        print(clu_key)
                        print(len(centroid.query))
                        print(len(_target))
                        print(_query.query.getTitle())
                    comb_dict[(wins, clu_key)].query_dict[win].append(_query)

        return


    def neighbor_calc_comb_score(self, comb_dict):
        '''
        The summed vdM score could not reflect the designability.
        Here is a new score method with weight added.

        In the future, the vdM score can be pre-calcluated as COMBS did.
        '''
        print('neighbor_calc_comb_score')

        for key in comb_dict.keys():       
            wins = key[0]
            clu_key = key[1]
            
            # From here we will calculate the score for each comb. 
            totals = [len(qs) for qs in comb_dict[key].query_dict.values()]
            try:
                #This is for the Search_selfcenter.
                if comb_dict[key].overlap_dict:
                    totals = [len(v) for v in comb_dict[key].overlap_dict.values()]
            except:
                print(totals) 
            #total_clu = sum([self.cluster_centroid_dict[clu_key[i]].total_clu for i in range(len(wins))])
            scores = [self.cluster_centroid_dict[clu_key[i]].score for i in range(len(wins))]
            comb_dict[key].totals = totals
            comb_dict[key].scores = scores

            #Calc fraction score
            fracScore = sum([self.cluster_centroid_dict[clu_key[i]].clu_num for i in range(len(wins))])/sum([self.cluster_centroid_dict[clu_key[i]].clu_total_num for i in range(len(wins))])
            weight = 0
            for i in range(len(wins)):
                c = self.cluster_centroid_dict[clu_key[i]]
                weight += sum(totals)/c.clu_total_num
            comb_dict[key].fracScore = fracScore*weight

            #Calc multiScore (By Bill: -ln(Na/SumNa * Nb/SumNb * Nc/SumNc))
            aa_numbs = [self.aa_num_dict[c[0]] for c in clu_key]
            comb_dict[key].multiScore = -log(np.prod(totals)/np.prod(aa_numbs))

        return
    

    def neighbor_calc_geometry(self, comb_dict):
        '''
        The overlap has several members for each query. 
        The geometry is a centroid contact atom of each query's candidates.
        Check CombInfo.calc_geometry()
        '''
        print('neighbor_calc_geometry')
        for wins, clu_key in comb_dict.keys():
            comb_dict[(wins, clu_key)].calc_geometry()         
        return


    def neighbor_aftersearch_filt(self, _target, comb_dict):
        '''
        After get the comb_dict, filter the pair angle, pair dists; filter the clash.
        remove the filtered key-value.
        '''
        for key in comb_dict.keys():  
            info = comb_dict[key]
            if self.search_filter.after_search_filter:
                if not info.after_search_condition_satisfied(self.search_filter.pair_angle_range, self.search_filter.pair_aa_aa_dist_range, self.search_filter.pair_metal_aa_dist_range):
                    comb_dict[key].after_search_filtered = True
                    
            if self.search_filter.filter_qt_clash:
                wins = [w for w in info.query_dict.keys()]
                vdms = [info.query_dict[w][0] for w in wins]

                if Search_filter.vdm_clash(vdms, _target, unsupperimposed=False, wins=wins):
                    comb_dict[key].vdm_no_clash = -1
                    comb_dict[key].after_search_filtered = True  
                else:
                    comb_dict[key].vdm_no_clash = 1       
        return


    def neighbor_write_represents(self):
        print('neighbor_write_represents')
        for key in self.neighbor_comb_dict.keys():  

            outdir = self.workdir + 'represents/'
            if not os.path.exists(outdir):
                os.mkdir(outdir)        

            tag = 'win_' + '-'.join([str(k) for k in key[0]]) + '_clu_' + '-'.join(k[0] + '-' + str(k[1]) for k in key[1]) 
            for w in key[0]:
                c = self.neighbor_comb_dict[key].query_dict[w][0]
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
            info = comb_dict[key]
            if not self.search_filter.write_filtered_result and info.after_search_filtered:
                continue

            outpath = 'win_' + '-'.join([str(k) for k in key[0]]) + '/'
            outdir = self.workdir + outpath
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            tag = 'win_' + '-'.join([str(k) for k in key[0]]) + '_clu_' + '-'.join(k[0] + '-' + str(k[1]) for k in key[1]) 

            # Write geometry       
            pr.writePDB(outdir + tag +'_geometry.pdb', comb_dict[key].geometry) 
            
            #Write 
            all_overlap_metal_coords = []
            for w in key[0]:
                candidates = comb_dict[key].query_dict[w]
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
            # self.neighbor_comb_dict[key].volume = volume
            # self.neighbor_comb_dict[key].volPerMetal = volume/len(all_overlap_metal_coords)
            # hdist = cdist(all_overlap_metal_coords, all_overlap_metal_coords, metric='euclidean')
            # self.neighbor_comb_dict[key].diameter = hdist.argmax()
            
            #Write Centroid and all metal coords in the cluster
            for w in key[0]:
                centroid = comb_dict[key].centroid_dict[w]
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


    def neighbor_write_summary(self, outdir, comb_dict, eval = False):
        '''
        Write a tab dilimited file.
        '''
        print('neighbor-write-summary')

        os.makedirs(outdir, exist_ok=True)

        with open(outdir + '_summary.tsv', 'w') as f:
            f.write('Wins\tClusterIDs\tproteinABPLEs\tCentroidABPLEs\tproteinPhiPsi\tCentroidPhiPsi\tvolume\tvol2metal\tdiameter\tTotalVdMScore\tFracScore\tMultiScore\taa_aa_dists\tmetal_aa_dists\tPair_angles\toverlap#\toverlaps#\tvdm_scores\ttotal_clu#\tclu_nums')
            f.write('\tpair_aa_aa_dist_ok\tpair_angle_ok\tpair_metal_aa_dist_ok\tvdm_no_clash')
            if eval:
                f.write('\teval_min_rmsd\teval_min_vdMs\teval_phi\teval_psi\teval_abple\teval_is_origin')
            f.write('\n')
            for key in comb_dict.keys(): 
                info = comb_dict[key]
                if not self.search_filter.write_filtered_result and info.after_search_filtered:
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

                f.write(str(info.pair_aa_aa_dist_ok) + '\t')
                f.write(str(info.pair_angle_ok) + '\t')
                f.write(str(info.pair_metal_aa_dist_ok) + '\t')
                f.write(str(info.vdm_no_clash) + '\t')

                if eval:
                    f.write('||'.join([str(round(m, 2)) for m in info.eval_mins]) + '\t')
                    f.write('||'.join([v.query.getTitle() for v in info.eval_min_vdMs]) + '\t')
                    f.write('||'.join([str(round(v.phi,2)) for v in info.eval_min_vdMs]) + '\t')
                    f.write('||'.join([str(round(v.psi,2)) for v in info.eval_min_vdMs]) + '\t')
                    f.write('||'.join([v.abple for v in info.eval_min_vdMs]) + '\t')
                    f.write(str(info.eval_is_origin) + '\t')

                f.write('\n')          
        return 


    def neighbor_write_log(self):
        if len(self.log) > 0:
            with open(self.workdir + '_log.txt', 'w') as f:
                f.write(self.log)

        with open(self.workdir + '_parameters.txt', 'w') as f:
            f.write(self.para2string())
            f.write(self.search_filter.para2string())



