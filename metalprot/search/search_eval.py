import os
import numpy as np
from numpy.core.numeric import count_nonzero
import prody as pr

from ..basic import vdmer
from ..basic import hull
from ..basic import utils
from ..database import core
from ..database import database_extract
from ..database import database_vdMAtomOrder

from sklearn.neighbors import NearestNeighbors
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool

from .search import Search_vdM, supperimpose_target_bb, calc_pairwise_neighbor
from .search_selfcenter import Search_selfcenter
from .graph import Graph
from .comb_info import CombInfo

class Search_eval(Search_selfcenter):


    def run_eval_selfcenter_search(self):
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
        win_combs, vdM_combs = self.eval_extract_closest_vdMs(wins, combs)

        ### No search but summary the result
        nature_comb_dict = self.eval_selfcenter_construct(win_combs[0], vdM_combs[0])

        if len(list(nature_comb_dict)) <= 0:
            return 

        #Specify the geometry for selfcenter search.
        nature_key = list(nature_comb_dict.keys())[0]
        self.geo_struct = nature_comb_dict[nature_key].geometry.copy()
        self.allowed_aa_combinations = [[v.aa_type for v in nature_comb_dict[nature_key].centroid_dict.values()]]


        if self.eval_mmdist:
            # Test metal-metal-dist
            for mmd in range(15, 50, 10):
                self.metal_metal_dist = mmd/100
                self.neighbor_comb_dict.clear()
                self.best_aa_comb_dict.clear()
                self.run_search_selfcenter()
                
                ### Evaluate search result.
                self.eval_search_results(wins, combs)
                self.find_best_for_nature_sel()
                self.neighbor_write_summary(self.workdir, self.neighbor_comb_dict, name = '_summary_' + str(mmd)  + '_' + self.target.getTitle() + '_' + self.time_tag + '.tsv', eval=True)
                self.neighbor_write_summary(self.workdir, self.best_aa_comb_dict, name = '_best_summary_' + str(mmd)  + '_' + self.target.getTitle() + '_' + self.time_tag + '.tsv', eval=True)

        if not self.eval_density:
            self.run_search_selfcenter()
            ### Evaluate search result.
            self.eval_search_results(wins, combs)
            self.find_best_for_nature_sel()
            self.neighbor_write_summary(self.workdir, self.neighbor_comb_dict, name = '_summary_' + self.target.getTitle() + '_' + self.time_tag + '.tsv', eval=True)
                
        self.neighbor_write_log()
        return 


    def eval_selfcenter_construct(self, win_comb, vdM_comb):
        '''
        Use the path to create comb_dict.
        '''
        clu_key = tuple([v.get_cluster_key() for v in vdM_comb])
        path = [v.id for v in vdM_comb]

        comb_dict = {}
        comb = dict()
        for i in range(len(win_comb)):
            comb[win_comb[i]] = [path[i]]
        combinfo = CombInfo()
        combinfo.comb = comb 
        comb_dict[(tuple(win_comb), clu_key)] = combinfo
        
        _target = self.target.copy()
        self.neighbor_extract_query(_target, comb_dict)

        self.neighbor_calc_geometry(comb_dict)
        
        if self.eval_density:
            #TO test different density radius.
            self.log += 'key\tradius\toverlap\tvolume\tdensity\tov1\tov2\tov3\ttotal_clu\tclu1\tclu2\tclu3\tf_total\tf_max\tf_avg\tf_median\n'
            for radius in range(20, 105, 5):
                self.selfcenter_calc_density(comb_dict, radius/100) 

        
        self.selfcenter_calc_density(comb_dict, self.density_radius)

        #self.comb_overlap(comb_dict)
        self.neighbor_calc_comb_score(comb_dict)

        
        # if len([comb_dict.keys()]) <= 0:
        #     return comb_dict

        self.selfcenter_write_win(comb_dict)
        outpath = 'closest_win_' + '_'.join([str(w) for w in win_comb]) + '/'
        outdir = self.workdir + outpath  
        self.neighbor_write_summary(outdir, comb_dict, name = '_closest_summary_' + self.target.getTitle() + '.tsv')

        return comb_dict


    def eval_get_comb(self):
        '''
        Extract vdMs from the target protein.
        '''
        aas = ['HIS', 'GLU', 'ASP', 'CYS']

        wins = []
        combs = []
        _cores = database_extract.get_metal_core_seq(self.target, vdmer.metal_sel, extend = 4)
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
            #Try to shit th oxy of asp and glu to the right order.
            _comb = []
            for c in comb:
                if vdmer.get_contact_atom(c).getResname() == 'GLU' or vdmer.get_contact_atom(c).getResname() == 'ASP':
                    database_vdMAtomOrder.asp_glu_oxy_shift(c)
                _comb.append(c)

            combs.append(_comb)
        return wins, combs


    def _eval_extract_closest_vdM(self, w, v):
        best_v = None
        best_id = -1
        min_rmsd = 10

        if not v:
            print('The vdM at {} from the protein bb is not successfully extracted.'.format(w))
            return best_v, best_id, min_rmsd

        coord = [v.select(vdmer.metal_sel)[0].getCoords()]

        ns = self.neighbor_query_dict[w].get_metal_mem_coords()

        x_in_y, x_has_y = calc_pairwise_neighbor(coord, ns, 1.5)


        for ind in x_in_y[0]:
        #for ind in range(len(self.vdms)):
            #TO DO: there is a bug of the vdM database. 
            #If the order of the atom names is not consistent (such as '-ND1 CD2-' in vdM but '-CD2 ND1-' in bb), it could not target the orginal vdM.
            if len(self.vdms[ind].query.select('heavy')) != len(v.select('heavy')):
                #print('supperimpose_target_bb not happening')
                continue
            #print(self.vdms[ind].query.getTitle())
            test_v = self.vdms[ind].copy()
    
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

        ## Plan to depre
        # if self.cluster_centroid_origin_dict:
        #     origin_centroid = self.cluster_centroid_origin_dict[clu_key].copy()
        #     supperimpose_target_bb(self.target, origin_centroid, w)

        #     #pr.writePDB(evaldir + tag + 'origin_' + origin_centroid.query.getTitle(), origin_centroid.query)
        #     clu_origin_allmetal_coords = origin_centroid.get_metal_mem_coords()
        #     hull.write2pymol(clu_origin_allmetal_coords, evaldir, tag + '_origin_w_' + str(w) +'_points.pdb') 

        centroid_id = self.cluster_centroid_dict[clu_key]
        centroid = self.vdms[centroid_id].copy()
        supperimpose_target_bb(self.target, centroid, w)

        pr.writePDB(evaldir + tag + centroid.query.getTitle(), centroid.query)
        clu_allmetal_coords = centroid.get_metal_mem_coords()
        hull.write2pymol(clu_allmetal_coords, evaldir, tag + '_w_' + str(w) +'_points.pdb') 

        origin_best_v = best_v.copy()
        transform = pr.calcTransformation(origin_best_v.query.select('heavy'), centroid.query.select('heavy'))
        transform.apply(origin_best_v.query)
        pr.writePDB(evaldir + tag + 'origin_' + origin_best_v.query.getTitle(), origin_best_v.query)

        return


    def eval_extract_closest_vdMs(self, wins, combs):
        '''
        First supperimpose the all_metal_vdm. 
        Then use nearest neighbor to get candidates by calculate the metal distance with dist 0.25.
        Then superimpose and obtain the one with min rmsd. 
        '''
        win_combs = []
        vdM_combs = []
        for i in range(len(wins)):
            evaldir = self.workdir + 'closest_win_' + '_'.join([str(w) for w in wins[i]])
            os.makedirs(evaldir, exist_ok=True)
            best_wins = []
            best_vdMs = []
            for j in range(len(wins[i])):
                w = wins[i][j]
                v = combs[i][j]
                best_v, best_id, min_rmsd = self._eval_extract_closest_vdM(w, v)               
                if not v:
                    continue
                pr.writePDB(evaldir + '/win_' + str(w) + '_' + v.getTitle(), v)
                if best_v:
                    best_wins.append(w)
                    best_vdMs.append(best_v)
                    clu_id = self.vdms[best_id].get_cluster_key()
                    tag = '/win_' + str(w) + '_clu_' + '_'.join([str(ci) for ci in clu_id]) + '_rmsd_' + str(round(min_rmsd, 3)) + '_'                   
                    print(tag + ' : best_id ' + str(best_id))
                    pr.writePDB(evaldir + tag + best_v.query.getTitle(), best_v.query)
                    self.write_closest_vdM_clu_points(best_v, w, evaldir, tag)
        win_combs.append(best_wins)
        vdM_combs.append(best_vdMs)
        return win_combs, vdM_combs

    
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

            for key in self.neighbor_comb_dict.keys():
                if not key[0] == tuple([w for w in wins[i]]): 
                    continue

                info = self.neighbor_comb_dict[key]
                if not self.search_filter.write_filtered_result and info.after_search_filtered:
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
  
                self.neighbor_comb_dict[key].eval_mins = min_rmsds
                self.neighbor_comb_dict[key].eval_min_vdMs = best_vs
                if all([m < 0.05 for m in min_rmsds]):
                    self.neighbor_comb_dict[key].eval_is_origin = True
        return


    def find_best_for_nature_sel(self):
        '''
        In the eval search, we can find the nature's selected vdm if the protein is included in the database. 
        The nature's selected vdms may be low overlap vdms, but they may be the member of a better vdm combinations.
        '''

        nature_key = [key for key in self.neighbor_comb_dict.keys() if self.neighbor_comb_dict[key].eval_is_origin == True][0]
        nature_info = self.neighbor_comb_dict[nature_key]

        for key in self.neighbor_comb_dict.keys():
            info = self.neighbor_comb_dict[key]
            if not self.search_filter.write_filtered_result and info.after_search_filtered:
                continue
            seen_here = []
            for w in key[0]:
                id = nature_info.centroid_dict[w].id
                if id not in info.overlap_query_id_dict[w]:
                    seen_here.append(False)
                else:
                    seen_here.append(True)
                
            if all(seen_here):
                self.neighbor_comb_dict[key].eval_contain_origin = True

        return
