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

from sklearn import neighbors

from ..basic import quco
from ..basic import hull
from ..basic import utils
from ..database import core
from ..database import database_extract

from sklearn.neighbors import NearestNeighbors
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool

from .search import Search_vdM
from .graph import Graph
from .comb_info import CombInfo

class Search_selfcenter(Search_vdM):
    '''
    Inheritated from Search_vdM. 
    Except here is searching the selfcenter vdM database. 
    '''

    def get_query_id_dict(self):
        #TO DO: Temp useage for query id extraction. the id will be in quco.query in the future.
        query_id_dict = {}
        for i in range(len(self.querys)):
            key = self.querys[i].query.getTitle()
            query_id_dict[key] = i
        return query_id_dict

    def run_neighbor_search(self):
        '''
        All functions need to run the neighbor search.
        '''
        print('run-selfcenter-neighbor-search')

        #TO DO: where should I apply filters: win filter, query_metal filter, phipsi, etc.
        self.neighbor_generate_query_dict()

        self.neighbor_generate_pair_dict()

        if self.parallel:
            self.neighbor_search_wins_pool()
        else:
            self.neighbor_search_wins()

        self.neighbor_write_summary(self.workdir, self.best_aa_comb_dict)

        self.neighbor_write_log()

        return


    def neighbor_run_comb(self, win_comb):
        # try:
        print('selfcenter-run at: ' + ','.join([str(w) for w in win_comb]))
        comb_dict = self.selfcenter_construct_comb(win_comb)
        if len([comb_dict.keys()]) <= 0:
            return comb_dict
            
        _target = self.target.copy()
        
        self.neighbor_extract_query(_target, comb_dict)

        self.neighbor_calc_geometry(comb_dict)
        self.comb_overlap(comb_dict)
        self.neighbor_calc_comb_score(comb_dict)

        self.neighbor_aftersearch_filt(_target, comb_dict)
        comb_dict = self.selfcenter_redu(comb_dict)
        
        if len([comb_dict.keys()]) <= 0:
            return comb_dict

        self.neighbor_write_win(comb_dict)
        
        self.selfcenter_get_write_represents(comb_dict)

        outpath = 'win_' + '-'.join([str(w) for w in win_comb]) + '/'
        outdir = self.workdir + outpath  
        
        if not self.search_filter.write_filtered_result:
            if len([key for key in comb_dict.keys() if not comb_dict[key].after_search_filtered]) > 0:
                self.neighbor_write_summary(outdir, comb_dict)
        else:
            self.neighbor_write_summary(outdir, comb_dict)

        return comb_dict
        # except:

        #     self.log += 'Error in win_comb: ' + '-'.join([str(w) for w in win_comb]) + '\n'
        #     return {}

    def selfcenter_construct_comb(self, win_comb):
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

        # path represent the id of each metal vdM.
        for path in graph.all_paths:
            
            clu_key = tuple([self.id_cluster_dict[p] for p in path])

            comb = dict()
            for i in range(len(win_comb)):
                comb[win_comb[i]] = [path[i]]

            combinfo = CombInfo()
            combinfo.comb = comb 
            comb_dict[(tuple(win_comb), clu_key)] = combinfo

        return comb_dict


    def comb_overlap(self, comb_dict):
        '''
        For each path, calc the overlap.
        '''
        for key in comb_dict.keys():

            win_comb = key[0]

            pair_dict = {}

            len_s = []
            for w in win_comb:
                ns = comb_dict[key].centroid_dict[w].get_hull_points()
                len_s.append(len(ns))

            for wx, wy in itertools.combinations(win_comb, 2):

                n_x = comb_dict[key].centroid_dict[wx].get_hull_points()
                n_y = comb_dict[key].centroid_dict[wy].get_hull_points()
                # print('-------------------')
                # print(wx)
                # print(wy)
                x_in_y, x_has_y = self.calc_pairwise_neighbor(n_x, n_y, self.selfcenter_rmsd)
                y_in_x, y_has_x = self.calc_pairwise_neighbor(n_y, n_x, self.selfcenter_rmsd)

                #if x_has_y and y_has_x: 
                pair_dict[(wx, wy)] = x_in_y
                pair_dict[(wy, wx)] = y_in_x

            graph = Graph(win_comb, len_s)

            graph.calc_pair_connectivity(pair_dict)

            graph.get_paths()

            overlap_dict = {} # {win:[overlap ids]} need to get {win:[query ids]} 

            for path in graph.all_paths:
                for i in range(len(win_comb)):
                    w = win_comb[i]
                    value = list(comb_dict[key].centroid_dict[w].selfcenter_cluster_queryid)[path[i]]
                    #value = path[i]
                    if w in overlap_dict.keys():
                        overlap_dict[w].add(value)
                    else:
                        overlap_dict[w] = set()
                        overlap_dict[w].add(value)

            comb_dict[key].overlap_dict = overlap_dict 

        return

    def selfcenter_redu(self, comb_dict):
        '''
        Here we try to remove any solution have seen in a better solution. The comb is in the overlap of a better comb.
        '''
        comb_dict_filter = {}

        query_id_dict = self.get_query_id_dict()

        comb_dict_sorted = {k: v for k, v in sorted(comb_dict.items(), key=lambda item: sum(item[1].totals), reverse = True)}

        #sum(list(comb_dict_sorted.values())[0].totals)

        if len(list(comb_dict_sorted.keys())) <= 0:
            return comb_dict_filter

        best_keys = [list(comb_dict_sorted.keys())[0]]

        for key in comb_dict_sorted:

            if key in best_keys:
                continue
            
            info = comb_dict_sorted[key]

            seen = []
            for bkey in best_keys:
                binfo = comb_dict_sorted[bkey]

                ### If the filter result is different, there is no need to check the members.
                if info.after_search_filtered != binfo.after_search_filtered:
                    seen.append(False)
                    continue

                seen_here = []
                for w in key[0]:
                    title = info.centroid_dict[w].query.getTitle()
                    id = query_id_dict[title]
                    if id not in binfo.overlap_dict[w]:
                        seen_here.append(False)
                    else:
                        seen_here.append(True)
                if all(seen_here):
                    seen.append(True)
                    break
                else:
                    seen.append(False)
                
            if not any(seen):
                best_keys.append(key)

        
        for key in best_keys:
            comb_dict_filter[key] = comb_dict[key]

        return comb_dict_filter


    def selfcenter_write_win(self, comb_dict):
        '''
        Write output.
        Too many output, May need optimization. 
        '''
        print('selfcenter search neighbor_write')
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
            
            #Write overlap
            for w in key[0]:
                candidate_ids = comb_dict[key].overlap_dict[w]
                print(len(candidate_ids))
                metal_coords = []

                max_out = 2 #To reduce number of output.
                for cid in candidate_ids:   
                    cquery = self.querys[cid].copy()
                    pdb_path = outdir + tag + cquery.query.getTitle() + '_w_' + str(w) + '_apble_' + cquery.abple + '.pdb'
                 
                    try:
                        pr.calcTransformation(cquery.query.select('heavy'), comb_dict[key].centroid_dict[w].query.select('heavy')).apply(cquery.query)
                        if max_out >0:
                            pr.writePDB(pdb_path, cquery.query)  
                            max_out-=1 

                    except:
                        print('Transformation failed:' + pdb_path)
                    
                    metal_coords.append(cquery.get_metal_coord())


                # centroid = comb_dict[key].centroid_dict[w]
                # clu_allmetal_coords = centroid.get_hull_points()
                # for cid in candidate_ids:
                #     #print(len(clu_allmetal_coords[cid]))
                #     metal_coords.append(clu_allmetal_coords[cid])

                hull.write2pymol(metal_coords, outdir, tag + '_w_' + str(w) +'_overlap_points.pdb')
            
            #Write Centroid and all metal coords in the cluster
            for w in key[0]:
                centroid = comb_dict[key].centroid_dict[w]
                pdb_path = outdir + tag + '_centroid_' + centroid.query.getTitle() + '.pdb'
                pr.writePDB(pdb_path, centroid.query)
                clu_allmetal_coords = centroid.get_hull_points()
                hull.write2pymol(clu_allmetal_coords, outdir, tag + '_w_' + str(w) +'_points.pdb')  
        return   
        

    def selfcenter_get_write_represents(self, comb_dict):
        '''
        Here is different from the parent class.
        Here we only want to write the 'best' comb in each win_comb with same aa_comb.
        '''
        print('selfcenter-write-represents.')
        outdir = self.workdir + 'represents/'
        if not os.path.exists(outdir):
            os.mkdir(outdir)      

        for key in comb_dict.keys():  
            info = comb_dict[key]
            if not self.search_filter.write_filtered_result and info.after_search_filtered:
                continue
            wins = key[0]
            aas = tuple([c[0] for c in key[1]])
            if (wins, aas) in self.best_aa_comb_dict.keys():
                if sum(self.best_aa_comb_dict[(wins, aas)].totals) < sum(comb_dict[key].totals):
                    self.best_aa_comb_dict[(wins, aas)] = comb_dict[key]
            else:
                self.best_aa_comb_dict[(wins, aas)] = comb_dict[key]
            
        
        for key in self.best_aa_comb_dict.keys():

            tag = 'win_' + '-'.join([str(k) for k in key[0]]) + '_clu_' + '-'.join(k[0] + '-' + str(k[1]) for k in key[1]) 
            for w in key[0]:
                c = self.best_aa_comb_dict[key].query_dict[w][0]
                pdb_path = outdir + tag + '_w_' + str(w) + '_apble_' + c.abple + '_' + c.query.getTitle() + '.pdb'
                pr.writePDB(pdb_path, c.query)                  

        return