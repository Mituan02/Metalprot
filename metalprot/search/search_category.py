'''
The database used for the search should be the category database which is not selfcentered.

The category search only implement the pairwise-neighbor method, 
which tried to find all combination wins first. 
Then for each win_comb, search the metal-metal distance.

The method may be deprecated in the future.
'''
import os
import itertools
from sklearn.neighbors import NearestNeighbors
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool
import prody as pr

from .search import Search_vdM, calc_pairwise_neighbor
from .graph import Graph
from .comb_info import CombInfo
from ..basic import hull
from ..basic import utils

class Search_category(Search_vdM):

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

        self.neighbor_write_summary(self.workdir, self.neighbor_comb_dict, name = '_summary_' + self.target.getTitle() + '_' + self.time_tag + '.tsv')

        self.neighbor_write_log()

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

                n_x = self.neighbor_query_dict[wx].get_metal_mem_coords()
                n_y = self.neighbor_query_dict[wy].get_metal_mem_coords()

                x_in_y, x_has_y = calc_pairwise_neighbor(n_x, n_y, self.metal_metal_dist)

                x_in_y, x_has_y = self.neighbor_filter_new(wx, wy, x_in_y)

                if x_has_y:
                    self.neighbor_pair_dict[(wx, wy)] = x_in_y

        return


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

                if not self.vdms[i].aa_type == resx:
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))] 
                    continue

            if self.search_filter.filter_abple:
                apx = self.target_abple[wx]
                apy = self.target_abple[wy]
       
                if not self.vdms[i].abple == apx:
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))] 
                    continue

            
            if self.search_filter.filter_phipsi:
                phix, psix = self.phipsi[wx]
                phiy, psiy = self.phipsi[wy]

                if (not utils.filter_phipsi(phix, self.vdms[i].phi, self.search_filter.max_phipsi_val)) or (not utils.filter_phipsi(psix, self.vdms[i].psi, self.search_filter.max_phipsi_val)):
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))]
                    continue
           

            if self.search_filter.filter_vdm_score:
                if self.vdms[i].score < self.search_filter.min_vdm_score:
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))]
                    continue
            

            if self.search_filter.filter_vdm_count:
                if self.vdms[i].clu_num < self.search_filter.min_vdm_clu_num:
                    x_in_y_mask[i] = [False for j in range(len(x_in_y[i]))]
                    continue
                           

            for j in range(len(x_in_y[i])):
                j_ind = x_in_y[i][j]

                if self.validateOriginStruct:
                    resy = self.target.select('name CA and resindex ' + str(wy)).getResnames()[0]
                    if not self.vdms[j_ind].aa_type == resy:
                        x_in_y_mask[i][j] = False
                        continue

                if self.search_filter.filter_abple:
                    apy = self.target_abple[wy]
                    if not self.querys[j_ind].abple == apy:
                        x_in_y_mask[i][j] = False
                        continue

                if self.search_filter.filter_phipsi:
                    phiy, psiy = self.phipsi[wy]
                    if (not utils.filter_phipsi(phiy, self.vdms[j_ind].phi, self.search_filter.max_phipsi_val)) or (not utils.filter_phipsi(psiy, self.vdms[j_ind].psi, self.search_filter.max_phipsi_val)) :
                        x_in_y_mask[i][j] = False
                        continue

                if self.search_filter.filter_vdm_score:
                    if self.vdms[j_ind].score < self.search_filter.min_vdm_score:
                        x_in_y_mask[i][j] = False 
                        continue
            

                if self.search_filter.filter_vdm_count:
                    if self.vdms[j_ind].clu_num < self.search_filter.min_vdm_clu_num:
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
                return comb_dict
        
        graph = Graph(win_comb,len(self.vdms))

        graph.calc_pair_connectivity(self.neighbor_pair_dict)

        graph.get_paths()

        print('graph.paths len {}'.format(len(graph.all_paths)))

        #TO DO: Here is a temp method to solve extream solutions. Mostly happened in 4 CYS binding cores.
        if len(graph.all_paths) > 15000:
            print('Too many paths to be considered so far.')
            graph.all_paths = graph.all_paths[0:15001]

        clu_dict = {}

        # path represent the id of each metal vdM.
        for path in graph.all_paths:
            
            clu_key = tuple([self.vdms[p].get_cluster_key() for p in path])
            
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


    def neighbor_construct_comb2(self, win_comb):
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
                return comb_dict
        
        graph = Graph(win_comb, len(self.vdms))

        graph.calc_pair_connectivity(self.neighbor_pair_dict)

        graph.get_paths()

        print('graph.paths len {}'.format(len(graph.all_paths)))

        #TO DO: Here is a temp method to solve extream solutions. Mostly happened in 4 CYS binding cores.
        if len(graph.all_paths) > 15000:
            print('Too many paths to be considered so far.')
            graph.all_paths = graph.all_paths[0:15001]

        # path represent the id of each metal vdM.
        count = 0
        for path in graph.all_paths:
            clu_key = tuple([self.vdms[p].get_cluster_key() for p in path])

            comb = dict()
            for i in range(len(win_comb)):
                comb[win_comb[i]] = [path[i]]

            combinfo = CombInfo()
            combinfo.comb = comb 
            comb_dict[(tuple(win_comb), clu_key, count)] = combinfo           
            count += 1

        return comb_dict


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
            
            #Write Centroid and all metal coords in the cluster
            for w in key[0]:
                centroid = comb_dict[key].centroid_dict[w]
                pdb_path = outdir + tag + '_centroid_' + centroid.query.getTitle() + '.pdb'
                pr.writePDB(pdb_path, centroid.query)

                if self.output_wincomb_overlap:
                    clu_allmetal_coords = centroid.get_metal_mem_coords()
                    hull.write2pymol(clu_allmetal_coords, outdir, tag + '_w_' + str(w) +'_points.pdb') 


                    candidates = comb_dict[key].query_dict[w]
                    metal_coords = []

                    max_out = 100 #To reduce number of output.
                    for c in candidates:
                        if max_out >0:
                            pdb_path = outdir + tag + '_w_' + str(w) + '_' + c.query.getTitle() + '.pdb'
                            pr.writePDB(pdb_path, c.query)                  
                            max_out-=1
                            metal_coords.append(c.get_metal_coord())
                        else:
                            break

                    hull.write2pymol(metal_coords, outdir, tag + '_w_' + str(w) +'_overlap_points.pdb')

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

