import os
from typing import Dict
from numba.cuda import target
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
from .quco import Query, Comb, pair_wise_geometry
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
                nodes.append(Node(i, j, self.win_inds, self.all_pos_len))
            self.nodes.append(nodes)
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



class CombInfo:
    def __init__(self):
        self.totals = None
        self.scores = []
        self.geometry = None
        self.query_dict = {}
        self.centroid_dict = {}

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

        return

        

class Search_vdM:
    '''
    The function to search comb
    '''
    def __init__(self, target_pdb, workdir, querys, id_cluster_dict, cluster_centroid_dict, all_metal_query, num_iter, rmsd = 0.25, win_filtered = None, 
    contact_querys = None, secondshell_querys = None, validateOriginStruct = False, filter_abple = False):

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

        self.target_abple = utils.seq_get_ABPLE(self.target)

        self.rmsd = rmsd

        self.win_filtered = win_filtered

        self.validateOriginStruct = validateOriginStruct

        self.filter_abple = filter_abple

        #neighbor searching strategy---------- 
        self.querys = querys             
        self.all_metal_query = all_metal_query #The query with all_metal_coord_ag
        self.id_cluster_dict = id_cluster_dict # {metal_id 1234: (HIS, 0)} 
        self.cluster_centroid_dict = cluster_centroid_dict #{(HIS, 0): centroid}


        self.neighbor_query_dict = dict() # {93: [the only centroid query with all metal coords]}
        self.neighbor_pair_dict = dict() # {(33, 37): [xs-33 near 37 coords]}
        self.neighbor_comb_dict = dict() 
        # { (wins, ids), (comb, combinfo)}
        # {((0, 1, 2, 3)(0, 0, 0, 0)): {(0:[1, 3, 4], 1:[2, 3, 4], 2: [2, 6, 7], 3:[1, 2, 3]), combinfo}} Please check neighbor_win2comb()
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

        self.neighbor_extract_query()

        self.neighbor_calc_comb_score()

        self.neighbor_calc_geometry()

        self.neighbor_write()

        self.neighbor_write_summary()

        return


    def neighbor_generate_query_dict(self):
        '''
        return self.neighbor_query_dict = dict() # {93: [the only centroid query with all metal coords]}

        '''
        print('hull_generate_query_dict')
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


    def neighbor_generate_pair_dict(self):
        '''
        To generate pair win dict, we will first filter out the impossible win pair by distance.
        Apply utils.check_pair_distance_satisfy()
        '''

        print('neighbor_generate_pair_dict')
        
        wins = list(self.neighbor_query_dict.keys())

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

                x_in_y, x_has_y = self.neighbor_filter(wx, wy, x_in_y)
                y_in_x, y_has_x = self.neighbor_filter(wy, wx, y_in_x)

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
            apx = self.target_abple[wx - 1]
            apy = self.target_abple[wy - 1]

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

    
    def neighbor_search_wins(self):
        '''
        The combinations of positions are extracted from all possible positions with pairs.
        '''
        wins = sorted(set(self.neighbor_query_dict.keys()))

        print('All wins with overlap {}'.format(wins))

        win_combs = itertools.combinations(wins, self.num_iter)
        for win_comb in win_combs:
            print(win_comb)
            
            self.neighbor_construct_comb(win_comb)

        return


    def neighbor_construct_comb(self, win_comb):
        '''
        win_comb: [0, 1, 2, 3]

        cluster_dict: {(0, 0, 0, 0): {0:[1, 3, 4], 1:[2, 3, 4], 2: [2, 6, 7], 3:[1, 2, 3]}}
        The key is win, the value is list of boolean, each represent one metal coord exist in all other wins' metal.
        
        self.neighbor_comb_dict
        # { (wins, ids), (comb, combinfo)}
        # {((0, 1, 2, 3)(0, 0, 0, 0)): {(0:[1, 3, 4], 1:[2, 3, 4], 2: [2, 6, 7], 3:[1, 2, 3]), combinfo}}


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

        print('graph.paths len {}'.format(len(graph.all_paths)))

        clu_dict = {}

        # path represent the id of each metal vdM.
        for path in graph.all_paths:
            
            clu_key = tuple([self.id_cluster_dict[p] for p in path])

            ### Debug purpose
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
                self.neighbor_comb_dict[(tuple(win_comb), clu_key)] = (comb, CombInfo())

        return


    def neighbor_extract_query(self):
        '''
        for each (comb, combinfo), we extract candidate querys and add it in combinfo.
        '''
        print('neighbor_calc_geometry')
        for wins, clu_key in self.neighbor_comb_dict.keys():
            for i in range(len(wins)):
                win = wins[i]
                self.neighbor_comb_dict[(wins, clu_key)][1].query_dict[win] = []

                clu = clu_key[i]
                centroid = self.cluster_centroid_dict[clu].copy()
                self.supperimpose_target_bb(centroid, win)
                self.neighbor_comb_dict[(wins, clu_key)][1].centroid_dict[win] = centroid

                for id in self.neighbor_comb_dict[(wins, clu_key)][0][win]:
                    _query = self.querys[id].copy()
                    self.supperimpose_target_bb(_query, win)
                    self.neighbor_comb_dict[(wins, clu_key)][1].query_dict[win].append(_query)

        return


    def neighbor_calc_comb_score(self):
        '''
        The summed vdM score could not reflect the designability.
        Here is a new score method with weight added.

        In the future, the vdM score can be pre-calcluated as COMBS did.
        '''
        print('hull_calc_comb_score')
        for wins, clu_key in self.neighbor_comb_dict.keys():       
            comb = self.neighbor_comb_dict[(wins, clu_key)][0] 

            # From here we will calculate the score for each comb. 
            totals = [len(comb[wins[i]]) for i in range(len(wins))]
            #total_clu = sum([self.cluster_centroid_dict[clu_key[i]].total_clu for i in range(len(wins))])
            score = 0
            self.neighbor_comb_dict[(wins, clu_key)][1].totals = totals
            self.neighbor_comb_dict[(wins, clu_key)][1].scores.append(score) 

        return

        

    def neighbor_calc_geometry(self):
        '''
        The overlap has several members for each query. 
        The geometry is a centroid contact atom of each query's candidates.
        Check CombInfo.calc_geometry()
        '''
        print('neighbor_calc_geometry')
        for wins, clu_key in self.neighbor_comb_dict.keys():
            self.neighbor_comb_dict[(wins, clu_key)][1].calc_geometry()         
        return
            

    def neighbor_write(self):
        '''
        Write output.
        Too many output, May need optimization. 
        '''
        print('neighbor_write')
        for key in self.neighbor_comb_dict.keys():       
            outpath = 'win_' + '-'.join([str(k) for k in key[0]]) + '/'
            outdir = self.workdir + outpath
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            tag = 'win_' + '-'.join([str(k) for k in key[0]]) + '_clu_' + '-'.join(str(k) for k in key[1]) 

            # Write geometry       
            pr.writePDB(outdir + tag +'_geometry.pdb', self.neighbor_comb_dict[key][1].geometry) 
            
            #Write pdbs
            for w in key[0]:
                candidates = self.neighbor_comb_dict[key][1].query_dict[w]
                metal_coords = []

                for c in candidates:
                    pdb_path = outdir + tag + '_' + c.query.getTitle()
                    pr.writePDB(pdb_path, c.query)

                    metal_coords.append(c.get_metal_coord())
                hull.write2pymol(metal_coords, outdir, tag + '_win_' + str(w) +'_points.pdb')

            #Write Centroid and all metal coords in the cluster
            for w in key[0]:
                centroid = self.neighbor_comb_dict[key][1].centroid_dict[w]
                pdb_path = outdir + tag + '_centroid_' + centroid.query.getTitle()
                pr.writePDB(pdb_path, centroid.query)
                clu_allmetal_coords = centroid.get_hull_points()
                hull.write2pymol(clu_allmetal_coords, outdir, tag + '_win_' + str(w) +'_points.pdb')  
        return    


    def neighbor_write_summary(self):
        '''
        Write a tab dilimited file.
        '''
        print('neighbor_write_summary')
        with open(self.workdir + '_summary.tsv', 'w') as f:
            f.write('Wins\tClusterIDs\tTotalScore\taa_aa_dists\tmetal_aa_dists\tPair_angles\toverlap#\toverlaps#\tvdm_scores\ttotal_clu#\tclu_nums\tCentroids\n')
            for key in self.neighbor_comb_dict.keys(): 
                info = self.neighbor_comb_dict[key][1]
                centroids = [c.query.getTitle() for c in info.centroid_dict.values()]
                vdm_scores = [c for c in info.scores]
                overlaps = [c for c in info.totals]
                clu_nums = [c.clu_num for c in info.centroid_dict.values()]
                
                f.write('_'.join([str(x) for x in key[0]]) + '\t')
                f.write('_'.join([str(x) for x in key[1]]) + '\t')
                f.write(str(round(sum(info.scores), 2)) + '\t')

                aa_aa_pair, metal_aa_pair, angle_pair  = pair_wise_geometry(info.geometry)
                f.write('||'.join([str(round(d, 2)) for d in aa_aa_pair])  + '\t')
                f.write('||'.join([str(round(d, 2)) for d in metal_aa_pair])  + '\t')
                f.write('||'.join([str(round(a, 2)) for a in angle_pair])  + '\t')

                f.write(str(sum(overlaps)) + '\t')
                f.write('||'.join([str(s) for s in overlaps]) + '\t')

                f.write('||'.join([str(round(s, 2)) for s in vdm_scores]) + '\t')
                f.write(str(sum(clu_nums)) + '\t')
                f.write('||'.join([str(c) for c in clu_nums]) + '\t')
                f.write('||'.join(centroids) + '\n')
        return 

