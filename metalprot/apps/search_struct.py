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


metal_sel = 'ion or name NI MN ZN CO CU MG FE' 


def _connectivity_filter(arr, inds):
    d = np.diff(arr[inds, :], axis=0)
    tf = (d[:, 0] == 1) & (d[:, 1] == 0) & (d[:, 2] == 0)
    return tf.all()


def connectivity_filter(pdb, window_inds):
    '''
    connectivity_filter copy and modified from qbits.filter
    '''
    N = len(pdb)
    arr = np.zeros((N, 3))
    resnums = pdb.getResnums().astype(np.int16)
    arr[:, 0] = resnums
    chains = pdb.getChids().astype('object')
    arr[:, 1] = np.array(list(map(ord, chains))) 
    segids = pdb.getSegnames().astype('object')
    #segids[segids == ''] = 'a' 
    segids[np.array([len(s)!=1 for s in segids])] = 'a'
    arr[:, 2] = np.array(list(map(ord, segids))) 
    return window_inds[[_connectivity_filter(arr, inds) for inds in window_inds]]


def target_position_filter(window_inds, select_inds):
    filter_window_inds = []
    for inds in window_inds:
        exist = True
        for ind in inds:
            if not ind in select_inds:
                exist = False
        if exist:
            filter_window_inds.append(inds)
    return np.array(filter_window_inds)


def supperimpose_target_bb(query, target, win_filter=None, rmsd_cut = 0.5, query_align_sel = 'name N CA C O', target_align_sel = 'name N CA C O', validateOriginStruct = False):
    '''
    Two possible way:1. Master search; 2. prody calcrmsd.
    query: prody pdb
    target: prody pdb
    '''
    try:
        query_len = len(query.query.select('name CA ' + query_align_sel))
        target_len = len(target.select('protein and name CA'))
    except:
        print(query.query.getTitle())
        return []
    ind = np.arange(query_len)
    window_inds = np.array([ind + i for i in range(target_len- query_len + 1)])

    window_inds = connectivity_filter(target.select('protein and name CA'), window_inds)

    if win_filter:
        window_inds = target_position_filter(window_inds, win_filter)

    new_querys = []
    for win in window_inds:
        new_query = query.copy()
        target_sel = target.select('resindex ' + ' '.join([str(w) for w in win]))
        if len(new_query.query.select(query_align_sel)) != len(target_sel.select(target_align_sel)):
            continue
        if validateOriginStruct and target_sel.select('name CA').getResnames()[0] != new_query.query.select('name CA').getResnames()[0]:
            continue
        #TO DO: The calcTransformation here will change the position of pdb. 
        #This will make the output pdb not align well. Current solved by re align.
        transform = pr.calcTransformation(new_query.query.select(query_align_sel), target_sel.select(target_align_sel))
        transform.apply(new_query.query)
        if new_query.hull_ag:
            transform.apply(new_query.hull_ag)
        rmsd = pr.calcRMSD(target_sel.select(target_align_sel), new_query.query.select(query_align_sel))
        if rmsd < rmsd_cut and not query_target_clash(new_query.query, win, target):
        #if rmsd < rmsd_cut:
            new_query.set_win(win)   
            new_querys.append(new_query)
    return new_querys

def simple_clash(query, query_2nd, clash_dist = 2.0):
    '''
    If the two query has CA within 2, then it is a crash.
    '''
    overlap = [w for w in query.win if w in query_2nd.win]
    if len(overlap) > 0:
        return True
    return False
    '''
    xyzs = []
    try:
        ni_index = query.query.select('ion or name NI MN ZN CO CU MG FE')[0].getIndex()
        all_near = query.query.select('protein and within 2.84 of index ' + str(ni_index))
        inds = all_near.select('nitrogen or oxygen or sulfur').getResindices()
        query_contact_sc = query.query.select('protein and heavy and resindex ' + ' '.join([str(ind) for ind in inds]))
    except:
        print('clashing fail: ' + query.getTitle())
        return True

    for c in query_contact_sc.getCoords():
        xyzs.append(c)
    try:
        ni_index2 = query_2nd.query.select('ion or name NI MN ZN CO CU MG FE')[0].getIndex()
        all_near2 = query_2nd.query.select('protein and within 2.84 of index ' + str(ni_index2))
        inds2 = all_near2.select('nitrogen or oxygen or sulfur').getResindices()
        query_contact_sc2 = query_2nd.query.select('protein and heavy and resindex ' + ' '.join([str(ind) for ind in inds2]))
    except:
        print('clashing fail: ' + query_2nd.getTitle())
        return True

    for c in query_contact_sc2.getCoords():
        xyzs.append(c)

    xyzs = np.vstack(xyzs)  
    dists = cdist(xyzs, xyzs)

    np.fill_diagonal(dists, 5)
    extracts = np.argwhere(dists <= clash_dist)

    first_len = len(query_contact_sc)
    extracts = [(ex[0], ex[1] - first_len) for ex in extracts if ex[0] < first_len and ex[1] > first_len]
    if len(extracts) > 0:
        return True
    return False
    '''

def query_target_clash(query, win, target, clash_dist = 2.0):
    xyzs = []

    #extract contact side chains. 
    try:
        ni_index = query.select('ion or name NI MN ZN CO CU MG FE')[0].getIndex()
        all_near = query.select('protein and within 2.83 of index ' + str(ni_index))
        inds = all_near.select('nitrogen or oxygen or sulfur').getResindices()
        query_contact_sc = query.select('sc and heavy and resindex ' + ' '.join([str(ind) for ind in inds]))
    except:
        print('clashing fail: ' + query.getTitle())
        return True


    for c in query_contact_sc.getCoords():
        xyzs.append(c)

    for c in target.select('bb and within 5 of resindex ' + ' '.join([str(ind) for ind in win]) + ' and not resindex ' + ' '.join([str(ind) for ind in win])).getCoords():
        xyzs.append(c)

    xyzs = np.vstack(xyzs)  
    dists = cdist(xyzs, xyzs)

    np.fill_diagonal(dists, 5)
    extracts = np.argwhere(dists <= clash_dist)

    first_len = len(query_contact_sc.select('sc and heavy'))
    extracts = [(ex[0], ex[1] - first_len) for ex in extracts if ex[0] < first_len and ex[1] > first_len]
    if len(extracts) > 0:
        return True
    return False

def geometry_filter(pdbs, contact_querys):
    '''
    There are only certain geometries for metal contact atoms.
    Only work for 3 aa + metal so far.
    '''
    contact_pdb = core.get_contact(pdbs)
    min_rmsd = 0.5
    min_query = None
    for query in contact_querys:
        if len(contact_pdb) != len(query.query):
            continue
        pr.calcTransformation(contact_pdb, query.query).apply(contact_pdb)
        rmsd = pr.calcRMSD(contact_pdb, query.query)

        if rmsd < min_rmsd:
            min_query = query
            min_rmsd = rmsd
    return min_query, min_rmsd
    

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

def check_hull_satisfy(target, wx, wy, query_x, query_y):
    '''
    For all vdms of the HIS GLU APS CYS, get the hull_ag which contains the metal points. 
    Then check the overlap of each pair on position wx, wy.
    query_x: The query contains hull_ag 
    query_y: The query contains hull_ag
    '''
    query_align_sel = 'name N CA C'


    inds = np.unique(query_x.query.getResindices())
    ind = inds[int(inds.shape[0]/2)-1]
    transform_x = pr.calcTransformation(query_x.query.select(query_align_sel + ' resindex ' + str(ind)), target.select(query_align_sel+ ' resindex ' + str(wx)))
    transform_x.apply(query_x.query)
    transform_x.apply(query_x.hull_ag)

    transform_y = pr.calcTransformation(query_y.query.select(query_align_sel+ ' resindex ' + str(ind)), target.select(query_align_sel+ ' resindex ' + str(wy)))
    transform_y.apply(query_y.query)
    transform_y.apply(query_y.hull_ag)
    
    hull_info = hull.calc_pairwise_hull(query_x.hull_ag.getCoords(), query_y.hull_ag.getCoords())
    if hull_info[0]:
        return True
    return False

def supperimpose_target_bb_bivalence(query, target, dist_array, id_array, win_filter = None, tolerance = 0.5, rmsd_cut = 0.5):
    '''
    Filter by contact map first with contact tolerance. 
    Then perform transformation and calculate RMSD. 
    return the Query array. 
    The win_filter is currently implemented in get_contact_map()
    '''
    cas = query.query.select('protein and name CA')
    dist = pr.calcDistance(cas[0], cas[1])
    left = np.searchsorted(dist_array, dist - tolerance)
    right = np.searchsorted(dist_array, dist + tolerance, side = 'right')

    new_querys = []
    if right >= left:
        for i in range(left, right+1):
            idi, idj = id_array[i]
            if win_filter and len(win_filter) > 0:
                if idi not in win_filter or idj not in win_filter:
                    continue
            target_sel = target.select('resindex ' + str(idi) + ' ' + str(idj))
            pr.calcTransformation(query.query.select('name N CA C O'), target_sel.select('name N CA C O')).apply(query.query)
            rmsd = pr.calcRMSD(target_sel.select('name N CA C O'), query.query.select('name N CA C O'))
            win = [idi, idj]
            if rmsd < rmsd_cut and query_target_clash(query, win, target):             
                new_query = Query(query.query.copy(), query.score, query.clu_num, query.clu_total_num, True, win, query.path)
                new_querys.append(new_query)

    return new_querys

def extract_win_filter_by_bivalence(query, target, dist_array, id_array, tolerance = 0.5, rmsd_cut = 0.5, validateOriginStruct = True):
    win_filter = set()
    cas = query.query.select('protein and name CA')
    dist = pr.calcDistance(cas[0], cas[1])
    left = np.searchsorted(dist_array, dist - tolerance)
    right = np.searchsorted(dist_array, dist + tolerance, side = 'right')

    if right >= left:
        for i in range(left, right+1):
            idi, idj = id_array[i]
            # if win_filter:
            #     if idi not in win_filter or idj not in win_filter:
            #         continue          
            target_sel = target.select('resindex ' + str(idi) + ' ' + str(idj))
            if validateOriginStruct:
                try:
                    if not target_sel.select('name CA').getResnames()[0] == cas.getResnames()[0]:
                        continue
                    if not target_sel.select('name CA').getResnames()[1] == cas.getResnames()[1]:
                        continue
                except:
                    print(query.query.getTitle())
                    continue
            if len(query.query.select('name N CA C O')) != len(target_sel.select('name N CA C O')):
                continue
            pr.calcTransformation(query.query.select('name N CA C O'), target_sel.select('name N CA C O')).apply(query.query)
            rmsd = pr.calcRMSD(target_sel.select('name N CA C O'), query.query.select('name N CA C O'))

            if rmsd < rmsd_cut:       
                #win_filter.add(idi)
                #win_filter.add(idj)
                win_filter.add((idi, idj))

    return win_filter


def extract_all_win_filter_by_bivalence(querys, target, tolerance = 0.5, rmsd_cut = 0.5, validateOriginStruct = True):
    win_filters = set()
    dist_array, id_array, dists = get_contact_map(target)
    for query in querys:
        win_filter = extract_win_filter_by_bivalence(query, target, dist_array, id_array, tolerance, rmsd_cut, validateOriginStruct)
        for w in win_filter:
            win_filters.add(w)
    print('win_filters {}'.format(win_filters))
    win_filter_inds = set([w for win in win_filters for w in win])
    #TO DO: The filter below is Not working yet.
    #win_filter_inds = filter_bivalence_win_by_comb(win_filters)
    return win_filter_inds

def pair_win_overlap_one(win_x, win_y):
    if win_x[0] in win_y and win_x[1] not in win_y:
        return True
    if win_x[1] in win_y and win_x[0] not in win_y:
        return True
    return False

def filter_bivalence_win_by_comb(win_filters_set):
    '''
    The idea is simple. A binding core must contain at least 3 contact amino acids if this is the case. 
    For example the extracted wins could be {[1, 2], [2, 3], [3, 1]} There mush have such an combination. 
    '''
    win_filter_inds = set() 
    win_filters = list(win_filters_set)
    for x in range(len(win_filters)):
        for y in range(x + 1, len(win_filters)):
            if not pair_win_overlap_one(win_filters[x], win_filters[y]):
                continue
            for z in range(y + 1, len(win_filters)):
                if pair_win_overlap_one(win_filters[x], win_filters[z]) and pair_win_overlap_one(win_filters[y], win_filters[z]):
                    win_filter_inds.add(win_filters[x][0])
                    win_filter_inds.add(win_filters[x][1])
                    win_filter_inds.add(win_filters[y][0])
                    win_filter_inds.add(win_filters[y][1])
                    win_filter_inds.add(win_filters[z][0])
                    win_filter_inds.add(win_filters[z][1])
    return win_filter_inds
                    

def generate_ind_combination_listoflist(_listoflist):
    _all_list = []
    for i in range(len(_listoflist)):
        _all_list.append(list(range(len(_listoflist[i]))))     
    all_inds = list(itertools.product(*_all_list))

    return all_inds

def get_combs_from_pair_extract(xys_len, extracts):

    xys = list(itertools.combinations(range(xys_len), 2))

    all_inds = generate_ind_combination_listoflist(extracts)

    comb_inds = []
    for inds in all_inds:
        extract = [extracts[j][inds[j]] for j in range(len(inds))]
        ext_inds = [None]*xys_len
        conflict = False
        for i in range(len(xys)):               
            x, y = xys[i]
            if ext_inds[x]:
                if ext_inds[x] != extract[i][0]:
                    conflict = True
                    break
            else:
                ext_inds[x] = extract[i][0]

            if ext_inds[y]:
                if ext_inds[y] != extract[i][1]:
                    conflict = True
                    break
            else:
                ext_inds[y] = extract[i][1]
        if not conflict:
            comb_inds.append(ext_inds)

    return comb_inds

def convert_query_2ndshellVdm(query):
    '''
    A query for 2nd shell can have different atom orders as an prody.atomGroup.
    The function is supposed to get the atomGroup of the query with the right atom order.
    '''

    metal = query.query.select(metal_sel)[0]
    metal_resind = metal.getResindex()

    contact_aa = query.query.select('protein and not carbon and not hydrogen and within 2.83 of resindex ' + str(metal_resind))
    _1stshell = query.query.select('name N C CA O and resindex ' + ' '.join([str(x) for x in contact_aa.getResindices()]))
    _1stshell_inds = _1stshell.getResindices()
    all_resinds = query.query.select('protein').getResindices()
    _2ndshell_resinds = [x for x in all_resinds if x not in _1stshell_inds]
    if len(_2ndshell_resinds) == 0:
        return None
    _2ndshell = query.query.select('name N C CA O and resindex ' + ' '.join([str(x) for x in _2ndshell_resinds]))

    neary_aas_coords = []
    neary_aas_coords.extend([x for x in _2ndshell.getCoords()])
    neary_aas_coords.extend([x for x in _1stshell.getCoords()])
    neary_aas_coords.append(metal.getCoords())
    coords = np.array(neary_aas_coords)

    names = []
    names.extend(_2ndshell.getNames())
    names.extend(_1stshell.getNames())
    names.append(metal.getName())

    ag = pr.AtomGroup(query.query.getTitle())
    ag.setCoords(coords)
    ag.setNames(names)

    return ag


def constructy_pseudo_2ndshellVdm(target, query, contact_resind):
    '''
    Find the resind of the vdm on target. Then extract resinds of atoms within a distance. 
    Followed by extracting the vdm resind and the atoms resind pairs together with the metal. 
    '''
        
    nearby_aas = target.select('protein and not carbon and not hydrogen and within 10 of resindex ' + str(contact_resind))
    nearby_aa_resinds = np.unique(nearby_aas.getResindices())    


    ags = []
    count = 0
    for resind in nearby_aa_resinds:
        if query.win and resind in query.win:
            continue
        neary_aas_coords = []
        neary_aas_coords.extend(target.select('name N C CA O and resindex ' + str(resind)).getCoords())
        neary_aas_coords.extend(query.query.select('bb or ion or name NI MN ZN CO CU MG FE').getCoords())
        coords = np.array(neary_aas_coords)

        names = []
        names.extend(target.select('name N C CA O and resindex ' + str(resind)).getNames())
        names.extend(query.query.select('bb or ion or name NI MN ZN CO CU MG FE').getNames())
        
        atom_contact_pdb = pr.AtomGroup('nearby_bb' + str(count))
        atom_contact_pdb.setCoords(coords)
        atom_contact_pdb.setNames(names)
        ags.append(atom_contact_pdb)
        count +=1

    return ags
    
    
def supperimpose_2ndshell(ag, query_2nd, rmsd_cut):
    '''
    supperimpose query to ag. 
    '''
    #print('supperimpose_2ndshell ' + query_2nd.query.getTitle())
    transform = pr.calcTransformation(query_2nd.ag, ag)
    transform.apply(query_2nd.ag)
    transform.apply(query_2nd.query)
    rmsd = pr.calcRMSD(ag, query_2nd.ag)

    if rmsd <= rmsd_cut:
        candidate = Query(query_2nd.query.copy(),  query_2nd.score, query_2nd.clu_num, query_2nd.clu_total_num, query_2nd.is_bivalent, query_2nd.win, query_2nd.path)
        return candidate
    return None


class Search_struct:
    '''
    The function to search comb
    '''
    def __init__(self, target_pdb, workdir, queryss, rmsd_cuts, dist_cuts, num_iter, qt_clash_dist, 
    qq_clash_dist, use_sep_aas, tolerance, fine_dist_cut = 0.3, win_filter = None, 
    contact_querys = None, secondshell_querys = None, validateOriginStruct = False,
    query_all_metal = None):
        if workdir:
            _workdir = os.path.realpath(workdir)
            if not os.path.exists(_workdir):
                os.mkdir(_workdir)
        else:
            _workdir = os.getcwd() + '/output_' + datetime.now().strftime('%Y-%m-%d-%H-%M-%S')          
            os.mkdir(_workdir)

        self.workdir = _workdir + '/'
        self.target = pr.parsePDB(target_pdb)
        self.rmsd_cuts = rmsd_cuts
        self.dist_cuts = dist_cuts
        self.num_iter = num_iter
        self.qt_clash_dist = qt_clash_dist
        self.qq_clash_dist = qq_clash_dist

        #Distance map for 2aa_sep database.
        self.dist_array, self.id_array, self.dists = get_contact_map(self.target, win_filter)
        self.use_sep_aas = use_sep_aas
        self.tolerance = tolerance

        self.fine_dist_cut = fine_dist_cut
        self.win_filter = win_filter
        self.validateOriginStruct = validateOriginStruct

        if len(queryss) < num_iter: 
            print('--You can only run win based search.')
        else:
            print('--You can run comb or win based search.')
        self.queryss = queryss

        #---------------------------------
        self.combs = []
        self.member_combs = dict()
        self.cquerysss = []  #check function generate_cquerys()
      
        xys = itertools.combinations(range(len(self.queryss)), 2)
        self.pair_extracts = [[0]*self.num_iter for i in range(self.num_iter)]
        for x, y in xys:
            self.pair_extracts[x][y] = dict()
        
        #contact-----------------------
        self.contact_querys = contact_querys

        #secondshell-----------------------
        self.secondshell_querys = secondshell_querys

        #Win based searching strategy----------
        self.win_querys = queryss[0]
        self.win_query_dict = dict() # 93: [<metalprot.apps.search_struct.Query at 0x7f7f58daeb20>, None, None, None, <metalprot.apps.search_struct.Query at 0x7f7f59275c10>]
        self.win_extract_dict = dict() # (33, 37): [(0, 24), (4, 24)]
        self.win_comb_dict = dict() # (33, 37, 56): [[0, 24, 0], [0, 24, 19], [4, 24, 0], [4, 24, 19]]

        #hull based searching strategy---------
        self.query_all_metal_x = None # Please check check_hull_satisfy()
        self.query_all_metal_y = None # Please check check_hull_satisfy()
        if query_all_metal:
            self.query_all_metal_x = query_all_metal
            self.query_all_metal_y = query_all_metal.copy()

        self.hull_query_dict = dict() # 93: [<metalprot.apps.search_struct.Query at 0x7f7f58daeb20>, None, None, None, <metalprot.apps.search_struct.Query at 0x7f7f59275c10>]
        self.hull_extract_dict = dict() # (33, 37): {(0, 24):(overlap, xs, ys), (4, 24):(overlap, xs,ys)}
        self.hull_comb_dict = dict() # Please check hull_overlap()
        self.hull_score_dict = dict()
        self.hull_geometry_dict = dict()

        #end---------------------------- 

    #region common functions

    def write_combs(self, combs, outpath = '/combs/'):      
        outdir = self.workdir + outpath
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        rank = 1
        for c in combs:
            count = 1
            for query in c.querys:
                pdb_path = outdir + str(rank) + '_' + str(count) + '_' + str(round(c.total_score, 2)) + '_' + query.query.getTitle()
                pr.writePDB(pdb_path, query.query)
                count+=1
            rank += 1


    def write_comb_info(self, combs, filename = '/_summary.txt'):
        with open(self.workdir + filename, 'w') as f:
            f.write('total_score\ttotal_clu_number\tquerys\tscores\tclu_nums\twins\tcontact_query\tcontact_score\tcontact_rmsd\tpair_dists\tpair_angles\n')
            for v in combs:
                f.write(v.to_tab_string() + '\n')  

    #endregion


    #region hull based search
    def run_hull_based_search(self):
        '''
        All functions need to run the hull based search.
        '''
        self.hull_generate_query_dict()
        self.hull_generate_pairwise_win_dict()

        self.hull_iter_win()

        self.hull_construct()       
        self.hull_calc_comb_score()
        self.hull_calc_geometry()
        self.hull_write()
        self.hull_write_summary()
        return

    
    def hull_generate_query_dict(self):
        '''
        self.query_dict is a dictionary [win, [query1, query2, query3, None,  ...]], where win is the target position. win could be int.
        '''
        # generate wins
        print('hull_generate_query_dict')
        wins = []
        if self.win_filter:
            wins.extend([w for w in self.win_filter])
        else:
            t = self.target.select('name CA').getResindices()
            wins.extend(([w for w in t]))
        
        #print('wins: {}'.format(wins))
        # for each win, add the querys into the self.hull_query_dict.
        for w in wins:     
            cquerys = []
            for ind in range(len(self.win_querys)):
                query = self.win_querys[ind]  
                cquery = supperimpose_target_bb(query, self.target, [w], self.rmsd_cuts[0], query_align_sel='name N CA C and resindex ' + str(query.contact_resind), target_align_sel='name N CA C', validateOriginStruct = self.validateOriginStruct)
                cquerys.append(cquery[0] if len(cquery) > 0 else None)
            if len(cquerys) > 0:
                self.hull_query_dict[w] = cquerys
        return


    def hull_generate_pairwise_win_dict(self):
        '''
        To generate pair win dict, we will first filter out the impossible win positions. 
        1. Filter by distance 
            check_pair_distance_satisfy()
        2. Filter by hull_ag overlap

        '''
        print('hull_generate_pairwise_win_dict')
        wins = list(self.hull_query_dict.keys())

        for inx in range(len(wins)):
            for iny in range(inx + 1, len(wins)):
                wx = wins[inx]
                wy = wins[iny]  
                dist_ok, ws = check_pair_distance_satisfy(wx, wy, self.dists)
                if not dist_ok:
                    continue
                # if self.query_all_metal and not check_hull_satisfy(self.target, wx, wy, self.query_all_metal_x, self.query_all_metal_y):
                #     continue
                extracts = self.hull_cal(wx, wy)
                if any(extracts):                 
                    self.hull_extract_dict[(wx, wy)] = extracts
        return


    def hull_cal(self, wx, wy):
        '''
        For each pair of aa positions, calculate the hull of all possible combination of members.
        extracts: (0, 24):(overlap, xs, ys)
        '''
        extracts = {}
        for i in range(len(self.hull_query_dict[wx])):
            for j in range(len(self.hull_query_dict[wy])):
                if not self.hull_query_dict[wx][i] or not self.hull_query_dict[wy][j]:
                    continue
                p1s = self.hull_query_dict[wx][i].get_hull_points()
                p2s = self.hull_query_dict[wy][j].get_hull_points()

                hull_info = hull.calc_pairwise_hull(p1s, p2s)
                if hull_info[0]:
                    extracts[(i, j)] = hull_info
        return extracts


    def hull_win2indcomb(self, win_comb):
        print('hull_win2indcomb')
        #For each pair of win combination, check if they satisfy dist_cut and contain extract. 
        _the_win_comb = []
        for i, j in itertools.combinations(range(self.num_iter), 2):
            wx = win_comb[i]
            wy = win_comb[j]
            if not (wx, wy) in self.hull_extract_dict.keys():
                return
            _the_win_comb.append((wx, wy))

        #For all extracts for the win_comb, check if any of the extract members are OK with the num_iter.
        extss = [self.hull_extract_dict[wc].keys() for wc in _the_win_comb]
        for xt in itertools.product(*extss):
            #xt: ((0, 0), (0, 2), (0, 2))
            ind_comb = dict()
            tuples = list(itertools.combinations(range(self.num_iter), 2))
            
            ab_ac_bc = True
            ind_xt = 0
            while ind_xt < len(xt) and ab_ac_bc:              
                i, j = tuples[ind_xt]              
                if i not in ind_comb.keys():
                    ind_comb[i] = xt[ind_xt][0]
                elif ind_comb[i] != xt[ind_xt][0]:
                    ind_comb[i] = -1
                    ab_ac_bc = False
                if j not in ind_comb.keys():
                    ind_comb[j] = xt[ind_xt][1]
                elif ind_comb[j] != xt[ind_xt][1]:
                    ind_comb[j] = -1
                    ab_ac_bc = False
                ind_xt +=1

            if not ab_ac_bc:
                continue
        
            self.hull_overlap(win_comb, _the_win_comb, tuple(ind_comb.values()), xt)
        return


    def hull_overlap(self, win_comb, win_comb_pair, clu_comb, clu_comb_pair):
        '''
        win_comb_pair: [(a, b), (a, c), (b, c)]. protein position a, b, c
        clu_comb_pair: ((0, 0), (0, 2), (0, 2)). for cluster a:0, b:0, c:2 

        hull_comb_dict[((a, b, c), (0, 0, 2))] = {(a, 0) : [T, F, T, F, F], (b, 0) : [T, F, T, F, F], (c, 2) : [T, F, T, F, F]}
        '''
        count = len(clu_comb_pair)

        overlap_mems = dict()
        for id in range(count):
            w1, w2 = win_comb_pair[id]
            c1, c2 = clu_comb_pair[id]

            firsts = self.hull_extract_dict[(w1, w2)][(c1, c2)][1]

            seconds = self.hull_extract_dict[(w1, w2)][(c1, c2)][2]

            if (w1, c1) in overlap_mems.keys():
                overlap_mems[(w1, c1)] = [overlap_mems[(w1, c1)][i] if firsts[i] else False for i in range(len(firsts))]
            else:
                overlap_mems[(w1, c1)] = firsts

            if (w2, c2) in overlap_mems.keys():
                overlap_mems[(w2, c2)] = [overlap_mems[(w2, c2)][i] if seconds[i] else False for i in range(len(seconds))]
            else:
                overlap_mems[(w2, c2)] = seconds
        # print(win_comb)
        # print(win_comb_pair)
        # print(clu_comb)
        # print(clu_comb_pair)
        if all([any(v) for v in overlap_mems.values()]):
            self.hull_comb_dict[(win_comb, clu_comb)] = overlap_mems
        return 

        
    def hull_iter_win(self):
        '''
        The combinations of positions are extracted from all possible positions with pairs.
        '''
        wins = list(self.hull_query_dict.keys())
        #wins = np.unique([k for key in self.hull_extract_dict.keys()for k in key]).sort()
        print('All wins with overlap {}'.format(wins))
        win_combs = itertools.combinations(wins, self.num_iter)
        for win_comb in win_combs:
            print(win_comb)
            self.hull_win2indcomb(win_comb)
        return


    def hull_construct(self):
        '''
        After the point overlap calculation to get self.hull_comb_dict. 
        Extract the member comb. 
        '''
        print('hull_construct comb')
        for key in self.hull_comb_dict.keys():       

            val = self.hull_comb_dict[key]
            tag = 'win_' + '-'.join([str(k) for k in key[0]]) + '_clu_' + '-'.join(str(k) for k in key[1])    
            for i in range(len(key[0])):
                win = key[0][i]
                clu = key[1][i] 
                query = self.hull_query_dict[win][clu]
                cand_bool = val[(win, clu)]
                query.candidates_metal_points = query.get_hull_points()[cand_bool]         
                query.realign_by_CCAN_candidates([i for i in range(len(cand_bool)) if cand_bool[i]])             
        return 


    def hull_calc_comb_score(self):
        '''
        The summed vdM score could not reflect the designability.
        Here is a new score method with weight added.
        '''
        print('hull_calc_comb_score')
        for key in self.hull_comb_dict.keys():
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


    def hull_calc_geometry(self):
        '''
        The hull overlap has several members for each query. The geometry is a centroid contact atom of each query's candidates.
        '''
        print('hull_calc_comb_score')
        for key in self.hull_comb_dict.keys():
            metal_coords = []
            coords = []   
            for i in range(len(key[0])):
                win = key[0][i]
                clu = key[1][i] 
                query = self.hull_query_dict[win][clu]
                query.extract_contact_atom()
                coords.append(pr.calcCenter(query.contact_ag))
                metal_coords.extend(query.candidates_metal_points)
            coords.append(pr.calcCenter(hull.transfer2pdb(metal_coords)))
            self.hull_geometry_dict[key] = hull.transfer2pdb(coords, ['NI' if i == len(coords)-1 else 'N' for i in range(len(coords))])
               
        return


    def hull_write(self):  
        print('hull_write')
        for key in self.hull_comb_dict.keys():       
            outpath = 'win_' + '-'.join([str(k) for k in key[0]]) + '/'
            outdir = self.workdir + outpath
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            tag = 'win_' + '-'.join([str(k) for k in key[0]]) + '_clu_' + '-'.join(str(k) for k in key[1])    
            for i in range(len(key[0])):
                win = key[0][i]
                clu = key[1][i] 
                query = self.hull_query_dict[win][clu]

                hull.write2pymol(query.candidates_metal_points, outdir, tag + '_win_' + str(win) +'_points.pdb')
                hull.write2pymol(query.get_hull_points(), outdir, tag + '_win_' + str(win) +'_all_points.pdb')

                pdb_path = outdir + tag + '_' + query.query.getTitle()
                pr.writePDB(pdb_path, query.query)

                ### If there are too many candidates, the output will be ugly. So it is better to write this only when specified. Or write it in a smart way.
                for c in query.candidates:
                    pdb_path = outdir + tag + '_w_' + str(win) + '_' + c.query.getTitle()
                    pr.writePDB(pdb_path, c.query)

            # Write geometry       
            pr.writePDB(outdir + tag +'_geometry.pdb', self.hull_geometry_dict[key])     
        return


    def hull_write_summary(self):
        print('hull_write_summary')

        with open(self.workdir + '_summary.tsv', 'w') as f:
            f.write('Wins\tClusterIDs\tTotalScore\taa_aa_dists\tmetal_aa_dists\tPair_angles\toverlap#\toverlaps#\tvdm_scores\ttotal_clu#\tclu_nums\tCentroids\n')
            for key in self.hull_comb_dict.keys(): 
                
                centroids = []
                vdm_scores = []
                overlaps = []
                clu_nums = []
                for i in range(len(key[0])):
                    win = key[0][i]
                    clu = key[1][i] 
                    query = self.hull_query_dict[win][clu]
                    centroids.append(query.query.getTitle())
                    vdm_scores.append(np.log(query.clu_num/query.clu_total_num))
                    overlaps.append(len(query.candidates))
                    clu_nums.append(query.clu_num)
                
                f.write('_'.join([str(x) for x in key[0]]) + '\t')
                f.write('_'.join([str(x) for x in key[1]]) + '\t')
                f.write(str(round(self.hull_score_dict[key], 2)) + '\t')
                aa_aa_pair, metal_aa_pair, angle_pair  = pair_wist_geometry(self.hull_geometry_dict[key])
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


    #region Win based search

    def run_win_based_search(self):

        self.generate_win_query_dict()

        self.iter_all_wins(self.dist_cuts[0]) 

        self.combs.extend(self.build_win_combs())

        self.write_combs(self.combs)

        self.write_comb_info(self.combs)
       

    def generate_win_query_dict(self):
        '''
        self.query_dict is a dictionary [win, [query1, query2, query3, None,  ...]], where win is the target position. win could be int or tuple for bivalent vdM.
        '''
        # check if we should consider bivalent vdm
        contain_bivalent_vdm = False
        ind = 0
        while ind < len(self.win_querys) - 1 and not contain_bivalent_vdm:
            query = self.win_querys[ind]
            if query.is_bivalent:
                contain_bivalent_vdm = True
            ind += 1

        # generate wins
        wins = []
        if self.win_filter:
            wins.extend([w for w in self.win_filter])
            if contain_bivalent_vdm:
                wins.extend(itertools.combinations(self.win_filter, 2))
        else:
            t = self.target.select('name CA').getResindices()
            wins.extend(([w for w in t]))
            if contain_bivalent_vdm:
                wins.extend(itertools.combinations(self.win_filter, 2)) 
        
        #print('wins: {}'.format(wins))
        # for each win, add the querys into the self.win_query_dict.
        for w in wins:     
            cquerys = []
            for ind in range(len(self.win_querys)):
                query = self.win_querys[ind]  
                if isinstance(w, tuple) and query.is_bivalent:
                    cquery = supperimpose_target_bb_bivalence(query, self.target, self.dist_array, self.id_array, w, self.tolerance, self.rmsd_cuts[0])
                elif not isinstance(w, tuple) and not query.is_bivalent:                    
                    cquery = supperimpose_target_bb(query, self.target, [w], self.rmsd_cuts[0], self.validateOriginStruct)
                cquerys.append(cquery[0] if len(cquery) > 0 else None)
            if len(cquerys) > 0:
                self.win_query_dict[w] = cquerys
        return


    def iter_all_wins(self, dist_cut):
        wins = list(self.win_query_dict.keys())
        for inx in range(len(wins)):
            for iny in range(inx + 1, len(wins)):
                wx = wins[inx]
                wy = wins[iny]  
                dist_ok, ws = check_pair_distance_satisfy(wx, wy, self.dists)
                if not dist_ok:
                    continue
                extracts = self.win_dist_cal(wx, wy, dist_cut)
                if len(extracts) > 0:
                    self.win_extract_dict[(wx, wy)] = extracts
        
        #TO DO: Not working for bivalent yet.
        win_combs = itertools.combinations(wins, self.num_iter)
        for win_comb in win_combs:
            #print(win_comb)
            ind_combs = self.win_2_ind_comb(win_comb)
            
            if len(ind_combs) > 0:
                self.win_comb_dict[win_comb] = ind_combs
        return


    def win_2_ind_comb(self, win_comb):
        ind_combs = []

        #For each pair of win combination, check if they satisfy dist_cut and contain extract. 
        _the_win_comb = []
        for i, j in itertools.combinations(range(self.num_iter), 2):
            wx = win_comb[i]
            wy = win_comb[j]
            if not (wx, wy) in self.win_extract_dict.keys():
                return ind_combs
            _the_win_comb.append((wx, wy))

        #For all extracts for the win_comb, check if any of the extract members are OK with the num_iter.
        extss = [self.win_extract_dict[wc] for wc in _the_win_comb]
        for xt in itertools.product(*extss):
            #xt: ((0, 0), (0, 2), (0, 2))
            ind_comb = dict()
            tuples = list(itertools.combinations(range(self.num_iter), 2))
            
            ab_ac_bc = True
            ind_xt = 0
            while ind_xt < len(xt) and ab_ac_bc:              
                i, j = tuples[ind_xt]              
                if i not in ind_comb.keys():
                    ind_comb[i] = xt[ind_xt][0]
                elif ind_comb[i] != xt[ind_xt][0]:
                    ind_comb[i] = -1
                    ab_ac_bc = False
                if j not in ind_comb.keys():
                    ind_comb[j] = xt[ind_xt][1]
                elif ind_comb[j] != xt[ind_xt][1]:
                    ind_comb[j] = -1
                    ab_ac_bc = False
                ind_xt +=1

            if not ab_ac_bc:
                continue
            ind_combs.append(list(ind_comb.values()))

        return ind_combs


    def win_dist_cal(self, wx, wy, dist_cut):

        xs = self.win_query_dict[wx]
        ys = self.win_query_dict[wy]

        xyzs = []
        for cq in xs:
            if not cq:
                xyzs.append(np.array([1000, 1000, 1000]))
            else:
                xyzs.append(cq.query.select('ion or name NI MN ZN CO CU MG FE')[0].getCoords())
        for cq in ys:
            if not cq:
                xyzs.append(np.array([2000, 2000, 2000]))
            else:
                xyzs.append(cq.query.select('ion or name NI MN ZN CO CU MG FE')[0].getCoords())

        xyzs = np.vstack(xyzs)  
        dists = cdist(xyzs, xyzs)

        np.fill_diagonal(dists, 5)
        extracts = np.argwhere(dists <= dist_cut)

        extracts = [(ex[0], ex[1] - len(xs)) for ex in extracts if ex[0] < len(xs) and ex[1] >= len(ys)]
        
        return extracts


    def build_win_combs(self): 
        '''
        build win based combs
        '''     
        combs = []
        for key in self.win_comb_dict.keys():
            for vs in self.win_comb_dict[key]:
                comb = []
                for i in range(len(key)):
                    k = key[i]
                    v = vs[i]
                    query = self.win_query_dict[k][v]                  
                    comb.append(query)
                combs.append(Comb(comb))
        if len(combs) > 0:
            combs.sort(key = lambda x: x.total_score, reverse = True) 
        return combs

    #endregion


    #region win based search for members

    def run_win_search_structure_member(self):
        print('combs count {}'.format(len(self.combs)))
        for i in range(len(self.combs)):          
            self.win_query_dict.clear()
            self.win_comb_dict.clear()
            wins = tuple(sorted([w for q in self.combs[i].querys for w in q.win]))
            if wins in self.member_combs.keys() and len(self.member_combs[wins]) > 0:
                continue

            for q in self.combs[i].querys:
                cvdms = []
                vdms = get_vdm_mem(q)
                for _query in vdms:            
                    target_sel = self.target.select('resindex ' + ' '.join([str(w) for w in _query.win]))
                    pr.calcTransformation(_query.query.select('name N CA C O'), target_sel.select('name N CA C O')).apply(_query.query)
                    cvdms.append(_query)
                self.win_query_dict[tuple(q.win)] = cvdms
            
            print(len(self.win_query_dict))           
            self.iter_all_wins(self.fine_dist_cut)
            print(len(self.win_comb_dict))
            combs = self.build_win_combs()          
            self.member_combs[wins] = combs
            print(len(self.member_combs))

        if len(self.member_combs) > 0:
            self.member_combs.sort(key = lambda x: x.total_score, reverse = True) 
        
        for c in self.member_combs:
            c.calc_pair_geometry()    

        self.write_combs(self.member_combs, outpath= '/mem_combs/')

        self.write_comb_info(self.member_combs, filename= '/_summary_mem.txt')
            
        return

    #endregion


    #region query based search

    def run_iter_search_structure(self):
        '''
        The searching step is follow '1st --> 2nd --> 3rd' 
        The searching speed is similar to run_search_struct. 
        '''     
        self.generate_cquerys(self.win_filter)
        comb_inds = self.get_iter_pair()

        self.combs.extend(self.build_combs(comb_inds))

        self.write_combs(self.combs, outpath= '/combs/')

        self.write_comb_info(self.combs)

    
    def get_iter_pair(self): 
        '''
        Get pair follow '1st-2nd --> 2nd-3rd --> 3rd-4th' 
        '''
        comb_inds = []
        for i in range(1, self.num_iter):
            all_inds = generate_ind_combination_listoflist(self.queryss[0:i+1])
            if len(comb_inds) > 0:
                all_inds = self.filter_all_inds(all_inds, comb_inds)
            comb_inds.clear()

            for inds in all_inds:
                extracts = self.get_pair_extracts(inds, self.dist_cuts)
                if extracts and len(extracts)>0:
                    combs = get_combs_from_pair_extract(i+1, extracts)
                    for comb in combs:
                        comb_inds.append((inds, comb))

        return comb_inds


    def filter_all_inds(self, all_inds, comb_inds):
        '''
        Used for get_iter_pair(). 
        For example, in each generation the generate_ind_combination_listoflist() will generate all possible combinations [1sts, 2nds, 3rds].
        This function will remove those that are not in [1sts, 2nds]. 
        If only [A, B1] is in comb_inds from last iteration, [A, B1, C] should be returned, but [A, B2, C] shouldn't. 
        '''
        inds_set = set([inds for (inds, comb) in comb_inds])

        filtered_all_inds = []
        for inds in all_inds:
            if tuple(inds[0:-1]) in inds_set:
                filtered_all_inds.append(inds)

        return filtered_all_inds


    def run_search_struct(self):
        self.generate_cquerys(self.win_filter)

        all_inds = generate_ind_combination_listoflist(self.queryss)
        print(len(all_inds))

        comb_inds = []
        for inds in all_inds:
            extracts = self.get_pair_extracts(inds, self.dist_cuts)
            if extracts and len(extracts)>0:
                combs = get_combs_from_pair_extract(self.num_iter, extracts)
                for comb in combs:
                    comb_inds.append((inds, comb))

        self.combs.extend(self.build_combs(comb_inds))

        self.write_combs(self.combs, outpath= '/combs/')

        self.write_comb_info(self.combs)


    def generate_cquerys(self, win_filter = None):
        '''
        self.cquerysss example
        [0, 1, 2]   num_iter
        [26, 26, 26]    candidate for each iter
        [88, 88, 88]    candidate possible superimpose match.

        '''
        for ind in range(len(self.queryss)):
            cqueryss = []
            for query in self.queryss[ind]:   
                if self.use_sep_aas[ind]:
                    cquerys = supperimpose_target_bb_bivalence(query, self.target, self.dist_array, self.id_array, win_filter, self.tolerance, self.rmsd_cuts[ind])
                else:
                    cquerys = supperimpose_target_bb(query, self.target, win_filter, self.rmsd_cuts[ind], self.validateOriginStruct)
                cqueryss.append(cquerys)
            self.cquerysss.append(cqueryss)


    def get_pair_extracts(self, inds, dist_cuts):

        xys = itertools.combinations(range(len(inds)), 2)

        extractss = []

        for x, y in xys:     
            #if already calculated, store
            if (inds[x], inds[y]) in self.pair_extracts[x][y].keys():
                #print('See it in the dictionary.')
                extracts_filtered = self.pair_extracts[x][y][(inds[x], inds[y])]
                if extracts_filtered:
                    extractss.append(extracts_filtered)
                    continue          
                else:
                    return
            
            cquerys_a = self.cquerysss[x][inds[x]]
            cquerys_b = self.cquerysss[y][inds[y]]

            if len(cquerys_a) <= 0 or len(cquerys_b) <=0: 
                #store
                self.pair_extracts[x][y][(inds[x], inds[y])] = None
                return

            extracts = self._metal_distance_extract(cquerys_a, cquerys_b, dist_cuts[y])
            extracts_filtered = self.distance_extracts_filter(cquerys_a, cquerys_b,extracts)
            #store
            self.pair_extracts[x][y][(inds[x], inds[y])] = extracts_filtered
            if extracts_filtered:
                extractss.append(extracts_filtered)
            else:
                return
        return extractss


    def _metal_distance_extract(self, cquerys_0, cquerys_ind, dist_cut):
        xyzs = []
        for cq in cquerys_0:
            xyzs.append(cq.query.select('ion or name NI MN ZN CO CU MG FE')[0].getCoords())
        for cq in cquerys_ind:
            xyzs.append(cq.query.select('ion or name NI MN ZN CO CU MG FE')[0].getCoords())

        xyzs = np.vstack(xyzs)  
        dists = cdist(xyzs, xyzs)

        np.fill_diagonal(dists, 5)
        extracts = np.argwhere(dists <= dist_cut)

        extracts = [(ex[0], ex[1] - len(cquerys_0)) for ex in extracts if ex[0] < len(cquerys_0) and ex[1] >= len(cquerys_0)]
        
        return extracts
    

    def distance_extracts_filter(self, cquerys_0, cquerys_ind, extracts):
        extracts_filtered = []
        if len(extracts) <= 0: 
            return
    
        for i, j in extracts:
            if simple_clash(cquerys_0[i], cquerys_ind[j], self.qq_clash_dist): 
                #print('two query clash.')
                continue
            if query_target_clash(cquerys_0[i].query, cquerys_0[i].win, self.target, self.qt_clash_dist) or query_target_clash(cquerys_ind[j].query, cquerys_ind[j].win, self.target, self.qt_clash_dist) :
                #print('query target clash.')
                continue
             
            extracts_filtered.append((i, j))
        if len(extracts_filtered) <= 0: 
            return
        return extracts_filtered

 
    def build_combs(self, comb_inds):
        '''
        comb_inds example: 
        [((18, 20, 22), [31, 59, 55]),
        ((18, 21, 22), [31, 60, 55])]
        
        '''
        combs = []
        check_dup = set()
        for inds, extracts in comb_inds:
            vdms = []
            for i in range(len(inds)):
                vdms.append(self.cquerysss[i][inds[i]][extracts[i]].copy())
            # TO DO: The remove duplicate here have bugs. The query.getTitle() may have same name but different position.
            if tuple([v.query.getTitle() for v in vdms]) in check_dup:
                continue
            for pm in itertools.permutations(range(len(vdms)), len(vdms)):
                check_dup.add(tuple([vdms[p].query.getTitle() for p in pm]))
            if self.contact_querys:
                min_query, min_rmsd = geometry_filter([v.query for v in vdms], self.contact_querys)           
                combs.append(Comb(vdms, min_query, min_rmsd))
            else:
                combs.append(Comb(vdms))
        if len(combs) > 0:
            combs.sort(key = lambda x: x.total_score, reverse = True)
        return combs 

    #endregion


    #region query based search for members

    def run_search_structure_member(self):
        '''
        Fine search step. To search into each cluster members. 
        '''
        #initialize self.cquerysss and self.pair_extracts
        self.cquerysss.clear()
        self.generate_cvdms(self.target)

        xys = itertools.combinations(range(len(self.queryss)), 2)
        self.pair_extracts = [[0]*self.num_iter for i in range(self.num_iter)]
        for x, y in xys:
            self.pair_extracts[x][y] = dict()

        all_inds = [[i]*self.num_iter for i in range(len(self.cquerysss[0]))]
        print(all_inds)

        comb_inds = []
        for inds in all_inds:
            extracts = self.get_pair_extracts(inds, dist_cuts = [self.fine_dist_cut]*self.num_iter)
            if extracts and len(extracts)>0:
                combs = get_combs_from_pair_extract(self.num_iter, extracts)
                for comb in combs:
                    comb_inds.append((inds, comb))

        self.combs.clear()
        
        self.combs.extend(self.build_combs(comb_inds))

        for c in self.combs:
            c.calc_pair_geometry()    

        self.write_combs(self.combs, outpath= '/mem_combs/')

        self.write_comb_info(self.combs, filename= '/_summary_mem.txt')

        return 

    def generate_cvdms(self, target):
        '''
        self.cvdmsss is in same structrue with self.cquerysss.
        cvdms is the members of each centroid candidates.
        '''
        for i in range(self.num_iter):
            cvdmss = []         
            for ind in range(len(self.combs)):
                query = self.combs[ind].querys[i] 
                cvdms = []             
                vdms = get_vdm_mem(query)
                for query in vdms:            
                    target_sel = target.select('resindex ' + ' '.join([str(w) for w in query.win]))
                    pr.calcTransformation(query.query.select('name N CA C O'), target_sel.select('name N CA C O')).apply(query.query)
                    cvdms.append(query)
                cvdmss.append(cvdms)
            self.cquerysss.append(cvdmss)


    #endregion


    #region bivalence based search

    def run_bivalence_search_structure(self):
        '''
        All querys in self.queryss are separate bivalence vdms.
        '''
        # self.generate_cquerys(self.win_filter)
                
        # comb_inds = self.get_bivalence_pair()

        # self.build_combs(comb_inds)

        # self.write_combs(outpath= '/combs/')

        # self.write_comb_info()

        self.generate_cquerys(self.win_filter)

        all_inds = generate_ind_combination_listoflist(self.queryss)
        print(len(all_inds))

        comb_inds = []
        for inds in all_inds:
            extracts = self.get_bivalence_extracts(inds)
            if extracts and len(extracts)>0:
                combs = get_combs_from_pair_extract(self.num_iter, extracts)
                for comb in combs:               
                    if self.overlap_all(inds, comb, self.num_iter):
                        #print(win)
                        comb_inds.append((inds, comb))

        self.combs.extend(self.build_combs(comb_inds))

        self.write_combs(self.combs, outpath= '/combs/')

        self.write_comb_info(self.combs)


    def get_bivalence_pair(self):
        '''
        Get ind follow '1st --> 2nd --> 3rd' 
        '''
        comb_inds = []
        for i in range(1, self.num_iter):
            all_inds = generate_ind_combination_listoflist(self.queryss[0:i+1])
            if len(comb_inds) > 0:
                all_inds = self.filter_all_inds(all_inds, comb_inds)
            comb_inds.clear()

            for inds in all_inds:
                extracts = self.get_bivalence_extracts(inds)      
                if extracts and len(extracts)>0:
                    #print(extracts)
                    combs = get_combs_from_pair_extract(i+1, extracts)

                    for comb in combs:
                        if len(comb) == self.num_iter: 
                            if self.overlap_all(inds, comb, self.num_iter):
                                #print(win)
                                comb_inds.append((inds, comb))
                        else:
                            comb_inds.append((inds, comb))
        #print(comb_inds)
        return comb_inds

    
    def get_bivalence_extracts(self, inds):
        '''
        For the bivalent vdMs, 
        '''
        xys = itertools.combinations(range(len(inds)), 2)
        extractss = []

        for x, y in xys:
            #if already calculated, store
            if (inds[x], inds[y]) in self.pair_extracts[x][y].keys():
                #print('See it in the dictionary.')
                extracts_filtered = self.pair_extracts[x][y][(inds[x], inds[y])]
                if extracts_filtered:                   
                    extractss.append(extracts_filtered)
                    continue          
                else:
                    return        

            cquerys_a = self.cquerysss[x][inds[x]]
            cquerys_b = self.cquerysss[y][inds[y]]

            if len(cquerys_a) <= 0 or len(cquerys_b) <=0: 
                #store
                self.pair_extracts[x][y][(inds[x], inds[y])] = None
                return

            #extracts = self._metal_distance_extract(cquerys_a, cquerys_b, dist_cuts[y])
            #extracts_filtered = self.distance_extracts_filter(cquerys_a, cquerys_b, extracts)
            extracts_filtered = self.bivalence_extract_filter(x, y, cquerys_a, cquerys_b, self.num_iter)

            #store
            self.pair_extracts[x][y][(inds[x], inds[y])] = extracts_filtered
            if extracts_filtered:
                extractss.append(extracts_filtered)
            else:
                return
        return extractss


    def bivalence_extract_filter(self, x, y, cquerys_a, cquerys_b, num_iter):
        '''
        Each of the bivalence vdM pair overlap with one of the amino acid, but not the other one.
        win_seen: dict(). {0:[1,2], 1:[2, 3], 2:[3,1]} is a valid one when self.num_iter==3.
        '''
        extracts = []
        win_seen = dict()
        for i in range(len(cquerys_a)):
            for j in range(len(cquerys_b)):
                win_seen.clear()
                win_seen[x] = cquerys_a[i].win
                win_seen[y] = cquerys_b[j].win
                if self.overlap(x, y, win_seen, self.num_iter):
                    extracts.append((i, j))
        return extracts


    def overlap(self, x, y, win_seen, num_iter):
        overlap = [value for value in win_seen[x] if value in win_seen[y]]    
        if len(overlap)>=2:
            return False

        if x + 1 == y:              
            if len(overlap)!= 1:
                return False
        elif x == 0 and y == num_iter -1:
            if len(overlap)!= 1:
                return False
        else:
            if len(overlap)>0:
                return False
        #print(x, y)
        #print(win_seen)
        return True


    def overlap_all(self, inds, comb, num_iter):
        win_seen = dict()                   
        xys = itertools.combinations(range(self.num_iter), 2)

        for x in range(self.num_iter):                              
            win_seen[x]=self.cquerysss[x][inds[x]][comb[x]].win
        #print(win_seen)

        if len(set([z for w in win_seen.values() for z in w])) > num_iter:
            return False

        for x, y in xys:
            if not self.overlap(x, y, win_seen, self.num_iter):
                return False
        return True

    #endregion bivalence based search


    #region 2nd shell search

    def run_search_2ndshells(self, outpath = '/mem_combs/', rmsd = 0.5):
        '''
        After find the self.combs, for each vdM in each combs, try to select the nearby aa bb to construct a pseudo 2ndshell vdM.
        Then compare with the 2ndshell vdM library. Keep the one within the rmsd limitation.
        '''
        for rank in range(len(self.combs)):
            self.search_2ndshell(rank, rmsd)

        for rank in range(len(self.combs)):
            self.write_2ndshell(rank, outpath)
    

    def search_2ndshell(self, rank, rmsd = 0.5):
        '''
        For the queries in each comb, search the 2nd shell vdms. Then store them in the query._2nd_shell. 
        '''
        comb = self.combs[rank]
        for query in comb.querys:
            contact_resind = query.win[int((len(query.win) - 1)/2)]
            ags = constructy_pseudo_2ndshellVdm(self.target, query, contact_resind)

            for ag in ags:                    
                candidates = self.search_2ndshellvmds(ag, rmsd)
                if len(candidates) > 0:
                    query._2nd_shells.extend(candidates)
        return


    def search_2ndshellvmds(self, ag, rmsd):
        candidates = []
        for query in self.secondshell_querys:
            if not query.ag:
                query_ag = convert_query_2ndshellVdm(query)
                if query_ag:
                    query.ag = query_ag
                else:
                    print('This query do not have 2nd shell: ' + query.query.getTitle())
                    continue
            candidate = supperimpose_2ndshell(ag, query, rmsd)
            if candidate:
                candidates.append(candidate)
        return candidates


    def write_2ndshell(self, rank, outpath = '/mem_combs/'):
        '''
        #Could be combined in write_comb_info.
        '''
        outdir = self.workdir + outpath
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        comb = self.combs[rank]
        count = 1
        for query in comb.querys:
            count2 = 1
            for _2ndshell in query._2nd_shells:
                pdb_path = outdir + str(rank + 1) + '_' + str(count) + '_2ndshell_' + str(count2) + '_' + str(round(_2ndshell.score, 2)) + '_' + _2ndshell.query.getTitle()
                pr.writePDB(pdb_path, _2ndshell.query)
                count2+=1
            count+=1

        return 

    #endregion 2nd shell search
