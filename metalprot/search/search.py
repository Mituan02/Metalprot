import os
import numpy as np
import itertools
from numpy.core.defchararray import multiply
import prody as pr
import datetime
import math
import string

from ..basic import hull
from ..basic import utils
from ..basic.filter import Search_filter
from .graph import Graph
from .comb_info import CombInfo
from ..basic import constant

from sklearn.neighbors import NearestNeighbors
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool


def supperimpose_target_bb(_target, _vdm, win, align_sel='name N CA C'):
    '''
    Copy the all_metal_vdm to a new object.
    Transform the copied all_metal_vdm to the target win. 
    '''

    target_sel = 'resindex ' + str(win) + ' and ' + align_sel
    query_sel = 'resindex ' + str(_vdm.contact_resind) + ' and '+ align_sel

    q = _vdm.query.select(query_sel)
    t = _target.select(target_sel)
    if len(q) != len(t):
        print('supperimpose-target-bb not happening')
        return False

    transform = pr.calcTransformation(q, t)
    transform.apply(_vdm.query)
    if _vdm.metal_atomgroup:
        transform.apply(_vdm.metal_atomgroup)  

    return True


def supperimpose_centroid(_vdm, centroid, align_sel='heavy'):
    '''
    supperimpose_centroid
    '''

    if len(_vdm.query.select(align_sel)) != len(centroid.query.select(align_sel)):
        print('supperimpose-centroid not happening')
        return False
    
    transform = pr.calcTransformation(_vdm.query.select(align_sel), centroid.query.select(align_sel))
    transform.apply(_vdm.query)
    if _vdm.metal_atomgroup:
        transform.apply(_vdm.metal_atomgroup) 
    return True


def calc_pairwise_neighbor(n_x, n_y, r = 0.25):
    '''
    Use sklean NearestNeighbors. find every x point has how many neighbors from y.
    '''

    neigh_y = NearestNeighbors(radius= r) 
    neigh_y.fit(n_y)

    x_in_y = neigh_y.radius_neighbors(n_x)
    x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])

    return x_in_y[1], x_has_y


def combine_vdm_into_ag(centroid_dict, key, tag, geometry = None, ideal_geometry = None):
    '''
    Merge all vdms from one CombInfo.centroid_dict in to an AtomGroup.
    '''
    ag = pr.AtomGroup(tag)
    coords = []
    chids = []
    names = []
    resnames = []
    resnums = []
    chain_num = 0
    for w in key[0]:
        c = centroid_dict[w].query.select('not name NI MN ZN CO CU MG FE' )
        c.setChids(string.ascii_uppercase[chain_num])
        coords.extend(c.getCoords())
        chids.extend(c.getChids())
        names.extend(c.getNames())
        resnames.extend(c.getResnames())
        resnums.extend(c.getResnums())
        chain_num += 1

    geometry.setChids(string.ascii_uppercase[chain_num])
    _geo = geometry.select('name NI')
    coords.extend(_geo.getCoords())
    chids.extend(_geo.getChids())
    names.extend(_geo.getNames())
    resnames.extend(_geo.getResnames())
    resnums.extend(_geo.getResnums())
    chain_num += 1

    if ideal_geometry:     
        ideal_geometry.setChids(string.ascii_uppercase[chain_num])
        coords.extend(ideal_geometry.getCoords())
        chids.extend(ideal_geometry.getChids())
        names.extend(ideal_geometry.getNames())
        resnames.extend(ideal_geometry.getResnames())
        resnums.extend(ideal_geometry.getResnums())

    ag.setCoords(np.array(coords))
    ag.setChids(chids)
    ag.setNames(names)
    ag.setResnames(resnames)
    ag.setResnums(resnums)
    return ag


class Search_vdM:
    '''
    The function to search comb is based on nearest neighbor function of sklearn.
    '''
    def __init__(self, target_pdb, workdir, vdms, cluster_centroid_dict, all_metal_vdm, 
    num_contact_vdms = [3], metal_metal_dist = 0.45, win_filtered = [], 
    validateOriginStruct = False, search_filter = None, geometry_path= None, parallel = False, density_radius = 0.6,
    secondshell_vdms = None, rmsd_2ndshell = 0.5, allowed_aa_combinations = [],
    output_wincomb_overlap = False):
        self.time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 
        if workdir:
            _workdir = os.path.realpath(workdir) + '_' +  self.time_tag          
        else:
            _workdir = os.getcwd() + '/output_' + self.time_tag              

        os.makedirs(_workdir, exist_ok=True)
        self.workdir = _workdir + '/'
        self.outdir_represent = self.workdir + 'represents/'
        os.makedirs(self.outdir_represent, exist_ok=True)  

        self.target = pr.parsePDB(target_pdb)

        self.num_contact_vdms = num_contact_vdms

        self._resnum_filtered = win_filtered


        self.target_abple, self.phipsi = utils.seq_get_ABPLE(self.target)

        self.metal_metal_dist = metal_metal_dist


        #neighbor in search filter-------------
        self.validateOriginStruct = validateOriginStruct
        self.allowed_aa_combinations = allowed_aa_combinations

        self.search_filter = search_filter
        if geometry_path:
            self.geo_struct = pr.parsePDB(geometry_path)
        else:
            self.geo_struct = constant.tetrahydra_geo

        #neighbor parallel mechanism-----------
        self.parallel = parallel

        #neighbor searching strategy-----------
        self.vdms = vdms #[query]
        self.all_metal_vdm = all_metal_vdm #The query with all_metal_coord_ag
        self.cluster_centroid_dict = cluster_centroid_dict #{(HIS, 0): centroid}
        self.neighbor_query_dict = dict() # {93: [the only centroid query with all metal coords]}

        #NearestNeighbor_pairwise dist strategy. Will be deprecated.
        self.neighbor_pair_dict = dict() # {(33, 37): [xs-33 near 37 coords]}
        self.neighbor_comb_dict = dict() # { (wins, ids), (comb, combinfo)}
        # exp. {((0, 1, 2, 3)(0, 0, 0, 0)): {(0:[1, 3, 4], 1:[2, 3, 4], 2: [2, 6, 7], 3:[1, 2, 3]), combinfo}} Please check neighbor_win2comb()
      
        #NearestNeighbor_graph strategy.


        #secondshell----------------------------
        self.secondshell_vdms = secondshell_vdms
        self.rmsd_2ndshell = rmsd_2ndshell

        #For multi scoring----------------------
        self.aa_num_dict = None
        self.aa_vdm_info_dict = {
            'H':[2662, 285, 92.08, 60],
            'E':[829, 50, 16.32, 11],
            'D':[896, 90, 26.12, 22],
            'C':[3211, 913, 457, 457]
        }


        #For Search_selfcenter-------------------
        self.density_radius = density_radius

        #Output control--------------------------
        self.log = '' #For developing output purpose
        self.best_aa_comb_dict = {} # To store&Write the best comb for each combinations of wins. 
        self.output_wincomb_overlap = output_wincomb_overlap

        #----------------------------------------
        self.setup()
        #end-------------------------------------


    def setup(self):
        '''
        Write parameters and calc some properties.
        '''
        with open(self.workdir + self.target.getTitle() + '_' + self.time_tag + '_parameters.txt', 'w') as f:
            f.write(self.para2string())
            f.write(self.search_filter.para2string())

        # reschain + resnum to resindex
        target_resnums = self.target.select('name CA').getResnums()
        target_chids = self.target.select('name CA').getChids()
        self.target_index_dict = {}
        resnum2ind = {}
        for ind in range(len(target_resnums)):
            if len(np.unique(target_chids)) != 1:
                chid_resnum = target_chids[ind] + '-' + str(target_resnums[ind])
            else:
                chid_resnum = str(target_resnums[ind])
            resnum2ind[chid_resnum] = ind
            self.target_index_dict[ind] = chid_resnum

        self.win_filtered = [resnum2ind[str(rn)] for rn in self._resnum_filtered]

        #allowed_aa_types example: [[H,H,H], [H,H,E], [H,H,D]]
        self.allowed_aas = set()
        self.allowed_aa_combinations_sorted = set()
        if len(self.allowed_aa_combinations) >0:
            for aas in self.allowed_aa_combinations:
                self.allowed_aa_combinations_sorted.add(tuple(sorted(aas)))
                for aa in aas:
                    self.allowed_aas.add(aa)
        print(self.allowed_aas)
        print(self.allowed_aa_combinations_sorted)
        # The contact map is used for the pair-wise neighbor method used before.
        # self.dist_array, self.id_array, self.dists = utils.get_contact_map(self.target, self.win_filtered)

        #Vdm database infomation.
        aa_num_dict = {}
        for v in self.vdms:
            aa = v.aa_type
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
        parameters += 'num_contact_vdms: ' + str(self.num_contact_vdms) + ' \n'
        parameters += 'metal_metal_dist: ' + str(self.metal_metal_dist) + ' \n'
        parameters += 'win_filtered: ' + str(self._resnum_filtered) + ' \n'
        parameters += 'validateOriginStruct: ' + str(self.validateOriginStruct) + ' \n'
        parameters += 'parallel: ' + str(self.parallel) + ' \n'
        parameters += 'density_radius: ' + str(self.density_radius) + ' \n'
        return parameters

    #region Neighbor Search

    def neighbor_generate_query_dict(self):
        '''
        return self.neighbor_query_dict = dict() # {93: [the only centroid query with all metal coords]}

        '''
        print('neighbor_generate_query_dict')
        wins = []
        if len(self.win_filtered) > 0:
            wins.extend([w for w in self.win_filtered])
        else:
            t = self.target.select('name CA').getResindices()
            wins.extend(([w for w in t]))

        for w in wins: 
            if self.validateOriginStruct:
                #here only filter aa, note the overlap that HIS still supperimpose to GLU.
                if self.target.select('resindex ' + str(w) + ' name CA').getResnames()[0] not in ['HIS', 'GLU', 'ASP', 'CYS']:
                    continue

            _vdm = self.all_metal_vdm.copy()
            x = supperimpose_target_bb(self.target, _vdm, w, align_sel='name N CA C')
            if x:
                self.neighbor_query_dict[w] = _vdm
        return


    def neighbor_extract_query(self, _target, comb_dict):
        '''
        for each (comb, combinfo), we extract candidate querys and add it in combinfo.
        '''
        print('neighbor-extract-query')
        for key in comb_dict.keys():
            wins = key[0]
            clu_key = key[1]
            for i in range(len(wins)):
                win = wins[i]
                comb_dict[key].query_dict[win] = []

                clu = clu_key[i]
                centroid_id = self.cluster_centroid_dict[clu]
                centroid = self.vdms[centroid_id].copy()

                supperimpose_target_bb(_target, centroid, win)
                comb_dict[key].centroid_dict[win] = centroid

                for id in comb_dict[key].comb[win]:
                    _query = self.vdms[id].copy()

                    supperimpose_target_bb(_target, _query, win)
                    comb_dict[key].query_dict[win].append(_query)

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

            info = comb_dict[key]
            # From here we will calculate the score for each comb. 
            overlaps = [len(info.query_dict[w]) for w in wins]
            try:
                #This is for the Search_selfcenter.
                if info.overlap_query_id_dict:
                    overlaps = [len(info.overlap_query_id_dict[w]) for w in wins]
            except:
                print('totals: ' + str(overlaps)) 
            
            scores = [self.vdms[self.cluster_centroid_dict[clu_key[i]]].score for i in range(len(wins))]

            core_types = [clu[0] for clu in clu_key]
            clu_medians = [self.aa_vdm_info_dict[c][3] for c in core_types]
            clu_totals = [self.aa_vdm_info_dict[c][0] for c in core_types]
            #cluters = [self.cluster_centroid_dict[clu_key[i]].total_clu for i in range(len(wins))]
            clu_nums = [c.clu_num for c in info.centroid_dict.values()]

            info.volume = 4.0/3.0 * math.pi * (self.density_radius**3)

            comb_dict[key].overlaps = overlaps
            comb_dict[key].scores = scores
            comb_dict[key].cluScore = sum([math.log(clu_nums[i]/clu_medians[i]) for i in range(len(wins))])
            comb_dict[key].overlapScore = np.prod([overlaps[i]*clu_totals[i]/clu_medians[i] for i in range(len(wins))])**(1/(len(wins)))/info.volume
            
        return
    

    def neighbor_calc_geometry(self, comb_dict):
        '''
        The overlap has several members for each query. 
        The geometry is a centroid contact atom of each query's candidates.
        Check CombInfo.calc_geometry()
        '''
        print('neighbor-calc-geometry')
        for key in comb_dict.keys():
            comb_dict[key].calc_geometry()         
        return


    def neighbor_aftersearch_filt(self, _target, comb_dict):
        '''
        After get the comb_dict, filter the pair angle, pair dists; filter the clash.
        remove the filtered key-value.
        '''
        for key in comb_dict.keys():  
            info = comb_dict[key]
            if self.search_filter.after_search_filter_geometry:
                if self.search_filter.filter_based_geometry_structure:
                    #TO DO: geometry structure based filter.
                    ideal_geometry, rmsd = Search_filter.get_min_geo(info.geometry, self.geo_struct)                    
                    if not Search_filter.after_search_geo_strcut_satisfied(info, ideal_geometry, self.search_filter.angle_tol, self.search_filter.aa_aa_tol, self.search_filter.aa_metal_tol):
                        comb_dict[key].after_search_filtered = True
                else:
                    if not Search_filter.after_search_geo_pairwise_satisfied(info, self.search_filter.pair_angle_range, self.search_filter.pair_aa_aa_dist_range, self.search_filter.pair_metal_aa_dist_range):
                        comb_dict[key].after_search_filtered = True
                    
            if self.search_filter.after_search_filter_qt_clash:
                wins = [w for w in info.query_dict.keys()]
                vdms = [info.query_dict[w][0] for w in wins]

                if Search_filter.vdm_clash(vdms, _target, unsupperimposed=False, wins=wins):
                    comb_dict[key].vdm_no_clash = -1
                    comb_dict[key].after_search_filtered = True  
                else:
                    comb_dict[key].vdm_no_clash = 1       
        return


    def neighbor_write_summary(self, outdir, comb_dict, name = '_summary.tsv', eval = False):
        '''
        Write a tab dilimited file.
        '''
        #print('write summary file.')

        os.makedirs(outdir, exist_ok=True)

        with open(outdir + name, 'w') as f:
            f.write('Wins\tClusterIDs\tDensityRadius\tCluScore\tOverlapScore\tOverlapScoreLn\tGeoRmsd\taa_aa_dists\tmetal_aa_dists\tPair_angles\toverlap#\toverlaps#\tclu_nums')
            f.write('\tpair_aa_aa_dist_ok\tpair_angle_ok\tpair_metal_aa_dist_ok\tvdm_no_clash\tproteinABPLEs\tCentroidABPLEs\tproteinPhiPsi\tCentroidPhiPsi')
            if eval:
                f.write('\teval_min_rmsd\teval_min_vdMs\teval_phi\teval_psi\teval_abple\teval_is_origin')
            f.write('\n')
            for key in comb_dict.keys(): 
                info = comb_dict[key]
                if not self.search_filter.write_filtered_result and info.after_search_filtered:
                    continue
                #centroids = [c.query.getTitle() for c in info.centroid_dict.values()]
                vdm_scores = [c for c in info.scores]
                overlaps = [c for c in info.overlaps]
                clu_nums = [c.clu_num for c in info.centroid_dict.values()]
                max_clu_nums = [c.max_clu_num for c in info.centroid_dict.values()]
                
                f.write('_'.join([self.target_index_dict[x] for x in key[0]]) + '\t')
                f.write('_'.join([x[0] + '-' + str(x[1]) for x in key[1]]) + '\t')

                f.write(str(round(self.density_radius, 2))+ '\t')
                f.write(str(round(info.cluScore, 2)) + '\t')
                f.write(str(round(info.overlapScore, 2)) + '\t')
                f.write(str(round(math.log(info.overlapScore), 2)) + '\t')
                f.write(str(round(info.geo_rmsd, 3)) + '\t')
                
                f.write('||'.join([str(round(d, 2)) for d in info.aa_aa_pair])  + '\t')
                f.write('||'.join([str(round(d, 2)) for d in info.metal_aa_pair])  + '\t')
                f.write('||'.join([str(round(a, 2)) for a in info.angle_pair])  + '\t')

                f.write(str(sum(overlaps)) + '\t')
                f.write('||'.join([str(s) for s in overlaps]) + '\t')
                f.write('||'.join([str(c) for c in clu_nums]) + '\t')

                f.write(str(info.pair_aa_aa_dist_ok) + '\t')
                f.write(str(info.pair_angle_ok) + '\t')
                f.write(str(info.pair_metal_aa_dist_ok) + '\t')
                f.write(str(info.vdm_no_clash) + '\t')

                f.write('_'.join([self.target_abple[x] for x in key[0]]) + '\t')
                f.write('_'.join([c.abple for c in info.centroid_dict.values()]) + '\t')
                f.write('_'.join([str((round(self.phipsi[x][0],2), round(self.phipsi[x][1],2))) for x in key[0]]) + '\t')
                f.write('_'.join([str((round(c.phi, 2), round(c.psi, 2))) for c in info.centroid_dict.values()]) + '\t')

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



