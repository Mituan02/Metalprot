import os
from typing import Dict
import numpy as np
import itertools
from numpy.lib.function_base import extract
import prody as pr
from scipy.spatial.distance import cdist, dice
import datetime
from .ligand_database import clu_info


class Query:
    def __init__(self, query, score, clu_num, clu_total_num):
        self.query = query
        self.score = score
        self.clu_num = clu_num
        self.clu_total_num = clu_total_num
        self.win = []

    def to_tab_string(self):
        query_info = self.query.getTitle() + '\t' + str(round(self.score, 2)) + '\t' + str(self.clu_num)  + '\t'+ str(self.clu_total_num)
        return query_info

class Comb:
    def __init__(self, querys):
        self.querys = querys       
        self.total_score = sum([q.score for q in querys])
        self.total_clu_number= sum([q.clu_num for q in querys])        
        self.scores = [q.score for q in querys]
        self.clu_nums= [q.clu_num for q in querys]

    def to_tab_string(self):
        vdm_info = str(round(self.total_score, 2)) + '\t' + str(self.total_clu_number) + '\t' + '||'.join([q.query.getTitle() for q in self.querys]) + '\t' + '||'.join([str(round(n, 2)) for n in self.scores]) + '\t' + '||'.join([str(s) for s in self.clu_nums]) 
        return vdm_info

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
                break
        if exist:
            filter_window_inds.append(inds)
    return np.array(filter_window_inds)


def supperimpose_target_bb(query, target, rmsd_cut = 0.5):
    '''
    Two possible way:1. Master search; 2. prody calcrmsd.
    query: prody pdb
    target: prody pdb
    '''
    try:
        query_len = len(query.query.select('protein and name CA'))
        target_len = len(target.select('protein and name CA'))
    except:
        print(query.query.getTitle())
        return []
    ind = np.arange(query_len)
    window_inds = np.array([ind + i for i in range(target_len- query_len + 1)])

    window_inds = connectivity_filter(target.select('protein and name CA'), window_inds)

    new_querys = []
    for win in window_inds:
        target_sel = target.select('resindex ' + ' '.join([str(w) for w in win]))
        if len(query.query.select('name N CA C O')) != len(target_sel.select('name N CA C O')):
            continue
        #TO DO: The calcTransformation here will change the position of pdb. 
        #This will make the output pdb not align well. Current solved by re align.
        pr.calcTransformation(query.query.select('name N CA C O'), target_sel.select('name N CA C O')).apply(query.query)
        rmsd = pr.calcRMSD(target_sel.select('name N CA C O'), query.query.select('name N CA C O'))

        if rmsd < rmsd_cut:
            new_query = Query(query.query.copy(), query.score, query.clu_num, query.clu_total_num)
            new_query.win.extend(win)          
            new_querys.append(new_query)

    return new_querys

def simple_clash(query, query_2nd, clash_dist = 2.0):
    '''
    If the two query has CA within 2, then it is a crash.
    '''
    xyzs = []
    try:
        ni_index = query.select('ion or name NI MN ZN CO CU MG FE')[0].getIndex()
        all_near = query.select('protein and within 2.84 of index ' + str(ni_index))
        inds = all_near.select('nitrogen or oxygen or sulfur').getResindices()
        query_contact_sc = query.select('protein and heavy and resindex ' + ' '.join([str(ind) for ind in inds]))
    except:
        print('clashing fail: ' + query.getTitle())
        return True

    for c in query_contact_sc.getCoords():
        xyzs.append(c)
    try:
        ni_index2 = query_2nd.select('ion or name NI MN ZN CO CU MG FE')[0].getIndex()
        all_near2 = query_2nd.select('protein and within 2.84 of index ' + str(ni_index2))
        inds2 = all_near2.select('nitrogen or oxygen or sulfur').getResindices()
        query_contact_sc2 = query_2nd.select('protein and heavy and resindex ' + ' '.join([str(ind) for ind in inds2]))
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

    for c in target.select('bb and not resindex ' + ' '.join([str(ind) for ind in win])).getCoords():
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

def geometry_filter():
    '''
    There are only certain geometries for metal contact atoms.
    '''

def get_contact_map(target):
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
    for i in range(len(xyzs[0])):
        for j in range(i+1, len(xyzs[0])):
            dist_array.append(dists[i, j])  
            id_array.append((i, j))

    dist_array, id_array = zip(*sorted(zip(dist_array, id_array)))
    return dist_array, id_array
    

def supperimpose_target_bb_bydist(query, target, dist_array, id_array, tolerance = 0.5, rmsd_cut = 0.5):
    '''
    Filter by contact map first with contact tolerance. 
    Then perform transformation and calculate RMSD. 
    return the Query array. 
    '''
    cas = query.query.select('protein and name CA')
    dist = pr.calcDistance(cas[0], cas[1])
    left = np.searchsorted(dist_array, dist - tolerance)
    right = np.searchsorted(dist_array, dist + tolerance, side = 'right')

    new_querys = []
    if right >= left:
        for i in range(left, right+1):
            idi, idj = id_array[i]
            target_sel = target.select('Resindex ' + str(idi) + ' ' + str(idj))
            pr.calcTransformation(query.query.select('name N CA C O'), target_sel.select('name N CA C O')).apply(query.query)
            rmsd = pr.calcRMSD(target_sel.select('name N CA C O'), query.query.select('name N CA C O'))

            if rmsd < rmsd_cut:
                new_query = Query(query.query.copy(), query.score, query.clu_num, query.clu_total_num)
                new_query.win.extend(win)          
                new_querys.append(new_query)

    return new_querys

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

class Search_struct:
    '''
    The function to search comb
    '''
    def __init__(self, target_pdb, workdir, queryss, rmsd_cuts, dist_cuts, num_iter, qt_clash_dist, qq_clash_dist, use_sep_aas, tolerance):
        if workdir:
            _workdir = os.path.realpath(workdir)
            if not os.path.exists(_workdir):
                os.mkdir(_workdir)
        else:
            _workdir = os.getcwd() + '/output_' + datetime.now().strftime('%Y-%m-%d-%H-%M-%S')          
            os.mkdir(_workdir)

        self.workdir = _workdir
        self.target = pr.parsePDB(target_pdb)
        self.rmsd_cuts = rmsd_cuts
        self.dist_cuts = dist_cuts
        self.num_iter = num_iter
        self.qt_clash_dist = qt_clash_dist
        self.qq_clash_dist = qq_clash_dist

        #Distance map for 2aa_sep database.
        dist_array, id_array = get_contact_map(self.target)
        self.dist_array = dist_array
        self.id_array = id_array
        self.use_sep_aas = use_sep_aas
        self.tolerance = tolerance

        if len(queryss) < num_iter:
            print('--Please includes the correct number of query list.')
        self.queryss = queryss

        #---------------------------------
        self.combs = []
        self.cquerysss = []

        xys = itertools.combinations(range(len(self.queryss)), 2)
        self.pair_extracts = [[None]* len(self.queryss)]* len(self.queryss)
        for x, y in xys:
            self.pair_extracts[x][y] = dict()

        #depre 
        #list of list of list(None), the len of list(None) equal to the len of queryss[0].
        self.extractsss = []
        for i in range(len(self.queryss)):
            extractss = []
            for j in range(len(self.queryss[i])):
                extractss.append([None for i in range(len(self.queryss[0]))])
            self.extractsss.append(extractss)
        self.queues = []
        

    def run_search_struct(self):
        self.generate_cquerys()

        all_inds = generate_ind_combination_listoflist(self.queryss)
        print(len(all_inds))

        comb_inds = []
        for inds in all_inds:
            extracts = self.get_pair_extracts(inds, dist_cut = 1)
            if extracts and len(extracts)>0:
                combs = get_combs_from_pair_extract(self.num_iter, extracts)
                for comb in combs:
                    comb_inds.append((inds, comb))

        self.build_combs(comb_inds)

        self.write_combs()

        self.write_comb_info()


    def generate_cquerys(self):
        for ind in range(len(self.queryss)):
            cqueryss = []
            for query in self.queryss[ind]:   
                if self.use_sep_aas[ind]:
                    cquerys = supperimpose_target_bb_bydist(query, self.target, self.dist_array, self.id_array, self.tolerance, self.rmsd_cuts[ind])
                else:
                    cquerys = supperimpose_target_bb(query, self.target, self.rmsd_cuts[ind])
                cqueryss.append(cquerys)
            self.cquerysss.append(cqueryss)


    def get_pair_extracts(self, inds, dist_cut):

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

            extracts = self._metal_distance_extract(cquerys_a, cquerys_b, dist_cut)
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
            if simple_clash(cquerys_0[i].query, cquerys_ind[j].query, self.qq_clash_dist): 
                #print('two query clash.')
                continue
            if query_target_clash(cquerys_0[i].query, cquerys_0[i].win, self.target, self.qt_clash_dist) or query_target_clash(cquerys_ind[j].query, cquerys_ind[j].win, self.target, self.qt_clash_dist) :
                #print('query target clash.')
                continue
             
            extracts_filtered.append((i, j))
        if len(extracts_filtered) <= 0: 
            return
        return extracts_filtered

 
    def build_combs(self, all_ind_extracts):
        
        for inds, extracts in all_ind_extracts:
            vdms = []
            for i in range(len(inds)):
                vdms.append(self.cquerysss[i][inds[i]][extracts[i]])
            self.combs.append(Comb(vdms))
        if len(self.combs) > 0:
            self.combs.sort(key = lambda x: x.total_score, reverse = True) 


    def write_combs(self):      
        outdir = self.workdir + '/combs/'
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        rank = 1
        for c in self.combs:
            count = 1
            for query in c.querys:
                pdb_path = outdir + str(rank) + '_' + str(count) + '_' + str(round(c.total_score, 2)) + '_' + query.query.getTitle()
                pr.writePDB(pdb_path, query.query)
                count+=1
            rank += 1


    def write_comb_info(self):
        with open(self.workdir + '/_summary.txt', 'w') as f:
            f.write('total_score\ttotal_clu_number\tquerys\tscores\tclu_nums\n')
            for v in self.combs:
                f.write(v.to_tab_string() + '\n')  
    

    def run_search_struct_depre(self):
        self.generate_cquerys()

        ind_exts = self.generate_combs_depre()

        self.build_combs(ind_exts)

        self.write_combs()

        self.write_comb_info()

    def generate_combs_depre(self):   
        _all_list = []
        for i in range(len(self.queryss)):
            _all_list.append(list(range(len(self.queryss[i]))))
        
        all_inds = list(itertools.product(*_all_list))

        all_ind_extracts = []
        for inds in all_inds:
            extractss = self.metal_distance_extract_depre(self.target, inds)
            if not extractss or len(extractss)<=0: continue
            for exts in extractss:
                all_ind_extracts.append((inds, exts))
        
        return all_ind_extracts

    def metal_distance_extract_depre(self, target, inds):
        '''
        Extract ids that satisfy distance limitation and clash filters.
        '''
        extract_binarys = []

        for index in range(1, len(inds)):

            extracts_filtered = self.extractsss[index][inds[index]][inds[0]]
            
            if extracts_filtered:
                extract_binarys.append(extracts_filtered)
                continue
            
            cquerys_0 = self.cquerysss[0][inds[0]]
            cquerys_ind = self.cquerysss[index][inds[index]]

            if len(cquerys_0) <= 0 or len(cquerys_ind) <=0: return

            extracts_filtered = self._metal_distance_extract(cquerys_0, cquerys_ind, self.dist_cuts[index])

            if not extracts_filtered: return 

            extract_binarys.append(extracts_filtered)
            self.extractsss[index][inds[index]][inds[0]] = extracts_filtered

        if len(inds) == 2:
            print(extract_binarys[0])
            return extract_binarys[0]

        #print(extract_binarys)
        extract_all_filtered = []
        for exb in list(itertools.product(*extract_binarys)):
            ex_1st = [x[0] for x in exb]
            if all(x == ex_1st[0] for x in ex_1st):  
                the_extract = [ex_1st[0]] + [x[1] for x in exb] 
                for i in range(1, len(inds)):
                    for j in range(i + 1, len(inds)):
                        if inds[i] == inds[j] and the_extract[i] == the_extract[j]:
                            return
                        if simple_clash(self.cquerysss[i][inds[i]][the_extract[i]].query, self.cquerysss[j][inds[j]][the_extract[j]].query, self.qq_clash_dist): 
                            print('two extra query clash.')
                            return

                print(the_extract)
                extract_all_filtered.append(the_extract)

        return extract_all_filtered