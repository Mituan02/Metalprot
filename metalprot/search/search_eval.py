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

from .search import Graph, CombInfo, Search_vdM


class Search_eval(Search_vdM):

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
            comb_dict = self.neighbor_run_comb(win_comb)
            if not comb_dict: continue
            self.neighbor_comb_dict.update(comb_dict)     

        #self.neighbor_write()

        #Evaluate search result.
        self.eval_search_results(wins, combs)

        self.neighbor_write_summary(self.workdir, self.neighbor_comb_dict, eval=True)

        return 
                

    def eval_get_comb(self):
        '''
        Extract vdMs from the target protein.
        '''
        aas = ['HIS', 'GLU', 'ASP', 'CYS']

        wins = []
        combs = []
        _cores = database_extract.get_metal_core_seq(self.target, quco.metal_sel, extend = 4)
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
                combinfo = self.neighbor_comb_dict[key]

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
                self.neighbor_comb_dict[key].eval_mins = min_rmsds
                self.neighbor_comb_dict[key].eval_min_vdMs = best_vs
                if all([m < 0.05 for m in min_rmsds]):
                    self.neighbor_comb_dict[key].eval_is_origin = True
        return