import os
import numpy as np
import itertools
import prody as pr
from scipy.spatial.distance import cdist
from .ligand_database import clu_info

class Query:
    def __init__(self, query, score, clu_num, clu_total_num):
        self.query = query
        self.score = score
        self.clu_num = clu_num
        self.clu_total_num = clu_total_num
        self.win = []

class Vdm:
    def __init__(self, querys):
        self.querys = querys       
        self.total_score = sum([q.score for q in querys])
        self.total_clu_number= sum([q.clu_num for q in querys])        
        self.scores = [q.score for q in querys]
        self.clu_nums= [q.clu_num for q in querys]

    def to_tab_string(self):
        vdm_info = str(round(self.total_score, 2)) + '\t' + str(self.total_clu_number) + '\t' + '||'.join([q.query.getTitle() for q in self.querys]) + '\t' + '||'.join([str(round(n, 2)) for n in self.scores]) + '\t' + '||'.join([str(s) for s in self.clu_nums]) 
        return vdm_info

def read_cluster_info(file_path):
    clu_infos = []
    if not os.path.exists(file_path):
        return clu_infos
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        for line in lines[1:]:
            ts = line[:-1].split('\t')
            if len(ts)!=7: continue
            clu_infos.append(clu_info(Metal=ts[0], clu_type = ts[1], clu_rmsd= ts[2], total_num=float(ts[3]), clu_rank= int(ts[4]), clu_num= int(ts[5]), score=float(ts[6])))
    return clu_infos

def filter_cluter_info(clu_infos, score_cut = 0, clu_num_cut = 10):
    filtered_infos = []
    for info in clu_infos:
        if info.score >= score_cut or info.clu_num >= clu_num_cut:
            filtered_infos.append(info)
    return filtered_infos

def extract_query(workdir, file_path = '_summary.txt', score_cut = 0, clu_num_cut = 10):
    clu_infos = read_cluster_info(workdir + file_path)
    filtered_infos = filter_cluter_info(clu_infos, score_cut, clu_num_cut)

    querys = []
    if len(filtered_infos) == 0:
        return querys

    for info in filtered_infos:
        centroid = [file for file in os.listdir(workdir  + str(info.clu_rank)) if 'centroid' in file][0]
        pbd = pr.parsePDB(workdir + str(info.clu_rank) + '/' + centroid)
        query = Query(pbd, info.score, info.clu_num, info.total_num)
        querys.append(query)

    return querys

def extract_centroid_pdb_in_clu(workdir, file_path = '_summary.txt', score_cut = 0, clu_num_cut = 10):  
    clu_infos = read_cluster_info(workdir + file_path)
    filtered_infos = filter_cluter_info(clu_infos, score_cut, clu_num_cut)
    pdb_paths = []
    if len(filtered_infos) == 0:
        return pdb_paths

    for info in filtered_infos:
        centroid = [file for file in os.listdir(workdir  + str(info.clu_rank)) if 'centroid' in file][0]
        pdb_paths.append(workdir + str(info.clu_rank) + '/' + centroid)
    return pdb_paths


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
    segids[segids == ''] = 'a' 
    arr[:, 2] = np.array(list(map(ord, segids))) 
    return window_inds[[_connectivity_filter(arr, inds) for inds in window_inds]]


def supperimpose_target_bb(query, target, rmsd_cut = 1):
    '''
    Two possible way:1. Master search; 2. prody calrmsd.
    query: prody pdb
    target: prody pdb
    '''
    query_len = len(query.query.select('protein and name CA'))
    target_len = len(target.select('protein and name CA'))

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
            new_querys.append(new_query)

    return new_querys

def metal_distance_extract(target, comb, query, dis_cut = 1):
    
    #only calculate the distance with the first vdm superimposed to the target.
    dist = pr.calcDistance(comb[0].query.select('ion'), query.query.select('ion'))

    if dist > dis_cut: 
        return 

    if query_target_clash(query.query, target):
            print('query target clash.')
            return
    for c in comb: 
        if simple_clash(c.query, query.query): 
            print('two query clash.')
            return
        
    comb.append(query)
    return comb 

def simple_clash(query, query_2nd):
    '''
    If the two query has CA within 2, then it is a crash.
    '''
    xyzs = []
    for c in query.select('protein and heavy').getCoords():
        xyzs.append(c)
    for c in query_2nd.select('protein and heavy').getCoords():
        xyzs.append(c)

    xyzs = np.vstack(xyzs)  
    dists = cdist(xyzs, xyzs)

    np.fill_diagonal(dists, 5)
    extracts = np.argwhere(dists <= 1.2)

    first_len = len(query.select('protein and heavy'))
    extracts = [(ex[0], ex[1] - first_len) for ex in extracts if ex[0] < first_len and ex[1] > first_len]
    if len(extracts) > 0:
        return True
    return False

def query_target_clash(query, target):
    xyzs = []
    for c in query.select('sc and heavy').getCoords():
        xyzs.append(c)
    for c in target.select('bb').getCoords():
        xyzs.append(c)

    xyzs = np.vstack(xyzs)  
    dists = cdist(xyzs, xyzs)

    np.fill_diagonal(dists, 5)
    extracts = np.argwhere(dists <= 1.2)

    first_len = len(query.select('sc and heavy'))
    extracts = [(ex[0], ex[1] - first_len) for ex in extracts if ex[0] < first_len and ex[1] > first_len]
    if len(extracts) > 0:
        return True
    return False

def geometry_filter():
    '''
    There are only certain geometries for metal contact atoms.
    '''

    
class Search_struct:
    '''
    The function to search comb
    '''
    def __init__(self, target_pdb, workdir, queryss, rmsd_cuts, dist_cuts, num_iter):
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

        if len(queryss) < num_iter:
            print('--Please includes the correct number of query list.')
        self.queryss = queryss

        #---------------------------------
        self.queues = []
        self.vdms = []

        self.cquerysss = []

        #list of list of list(None), the len of list(None) equal to the len of queryss[0].
        self.extractsss = []
        for i in range(len(self.queryss)):
            extractss = []
            for j in range(len(self.queryss[i])):
                extractss.append([None for i in range(len(self.queryss[0]))])
            self.extractsss.append(extractss)


    def run_search_struct(self):
        self.generate_cquerys()

        ind_exts = self.generate_combs()

        self.build_combs(ind_exts)

        self.write_combs()

        self.write_comb_info()

    def generate_cquerys(self):
        for ind in range(len(self.queryss)):
            cqueryss = []
            for query in self.queryss[ind]:    
                cquerys = supperimpose_target_bb(query, self.target, self.rmsd_cuts[ind])
                cqueryss.append(cquerys)
            self.cquerysss.append(cqueryss)


    def generate_combs(self):   
        _all_list = []
        for i in range(len(self.queryss)):
            _all_list.append(list(range(len(self.queryss[i]))))
        
        all_inds = list(itertools.product(*_all_list))

        all_ind_extracts = []
        for inds in all_inds:
            extractss = self.metal_distance_extract(self.target, inds)
            if not extractss or len(extractss)<=0: continue
            for exts in extractss:
                all_ind_extracts.append((inds, exts))
        
        return all_ind_extracts

    def _metal_distance_extract(self, cquerys_0, cquerys_ind, dist_cut):
        xyzs = []
        for cq in cquerys_0:
            xyzs.append(cq.query.select('ion').getCoords())
        for cq in cquerys_ind:
            xyzs.append(cq.query.select('ion').getCoords())

        xyzs = np.vstack(xyzs)  
        dists = cdist(xyzs, xyzs)

        np.fill_diagonal(dists, 5)
        extracts = np.argwhere(dists <= 1)

        extracts = [(ex[0], ex[1] - len(cquerys_0)) for ex in extracts if ex[0] < len(cquerys_0) and ex[1] > len(cquerys_0)]
        
        if len(extracts) <= 0: 
            return 
    
        extracts_filtered = []
        for i, j in extracts:
            if simple_clash(cquerys_0[i].query, cquerys_ind[j].query): 
                print('two query clash.')
                continue
            if query_target_clash(cquerys_0[i].query, self.target) or query_target_clash(cquerys_ind[j].query, self.target) :
                print('query target clash.')
                continue
            extracts_filtered.append([i, j])
        if len(extracts_filtered) <= 0: 
            return 
        return extracts_filtered


    def metal_distance_extract(self, target, inds):
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
            return extract_binarys[0]

        print(extract_binarys)
        extract_all_filtered = []
        for exb in list(itertools.product(*extract_binarys)):
            print(exb)
            ex_1st = [x[0] for x in exb]
            if all(x == ex_1st[0] for x in ex_1st):  
                the_extract = [ex_1st[0]] + [x[1] for x in exb] 
                print(the_extract)
                for i in range(1, len(inds)):
                    for j in range(i + 1, len(inds)):
                        if simple_clash(self.cquerysss[i][inds[i]][the_extract[i]].query, self.cquerysss[j][inds[j]][the_extract[j]].query): 
                            print('two extra query clash.')
                            return

                extract_all_filtered.append(the_extract)

        return extract_all_filtered


    def build_combs(self, all_ind_extracts):
        
        for inds, extracts in all_ind_extracts:
            vdms = []
            for i in range(len(inds)):
                vdms.append(self.cquerysss[i][inds[i]][extracts[i]])
            self.vdms.append(Vdm(vdms))
        if len(self.vdms) > 0:
            self.vdms.sort(key = lambda x: x.total_score, reverse = True) 


    def write_combs(self):      
        outdir = self.workdir + '/combs/'
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        rank = 1
        for c in self.vdms:
            count = 1
            for query in c.querys:
                pdb_path = outdir + str(rank) + '_' + str(count) + '_' + str(round(c.total_score, 2)) + '_' + query.query.getTitle()
                pr.writePDB(pdb_path, query.query)
                count+=1
            rank += 1

    def write_comb_info(self):
        with open(self.workdir + '/_summary.txt', 'w') as f:
            f.write('total_score\ttotal_clu_number\tquerys\tscores\tclu_nums\n')
            for v in self.vdms:
                f.write(v.to_tab_string() + '\n')  
    

    def generate_combs_depre(self):   
        for query in self.queryss[0]:
            new_querys = supperimpose_target_bb(query, self.target, self.rmsd_cuts[0])
            for q in new_querys:
                self.queues.append(([q], self.num_iter - 1))

        while len(self.queues) > 0:
            print('--In the queues')
            queue = self.queues.pop(0)
            self.construct_comb(queue[0], queue[1])
        
        self.write_combs()


    def construct_comb_depre(self, comb, num_iter):
        ind = self.num_iter - num_iter

        if ind == self.num_iter:
            print('Generate final combination.')
            self.vdms.append(Vdm(comb))
            return

        for query in self.queryss[ind]:
            _querys = supperimpose_target_bb(query, self.target, self.rmsd_cuts[ind])
            for q in _querys:
                _comb = metal_distance_extract(self.target, comb, q, self.dist_cuts[ind])
                if not _comb: continue
                self.queues.append((_comb, num_iter-1))


#The following functions are deprecated. 
def write_query_pdbs_depre(outdir, win_extract):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    count = 1
    for wx in win_extract:
        pr.writePDB(outdir + 'ex_' + str(count) + '_' + wx[0].getTitle(), wx[0])
        count += 1

def metal_distance_extract_depre_(target, win_extract, win_extract_2nd, distance_cut = 1):
    xyzs = []
    for win in win_extract:
        xyzs.append(win[0].select('ion').getCoords())
    for win in win_extract_2nd:
        xyzs.append(win[0].select('ion').getCoords())

    xyzs = np.vstack(xyzs)  
    dists = cdist(xyzs, xyzs)

    np.fill_diagonal(dists, 5)
    extracts = np.argwhere(dists <= 1)

    extracts = [(ex[0], ex[1] - len(win_extract)) for ex in extracts if ex[0] < len(win_extract) and ex[1] > len(win_extract)]
    
    extracts_filtered = []
    for i, j in extracts:
        if simple_clash(win_extract[i][0], win_extract_2nd[j][0]): 
            print('two query clash.')
            continue
        if query_target_clash(win_extract[i][0], target) or query_target_clash(win_extract_2nd[j][0], target) :
            print('query target clash.')
            continue
        extracts_filtered.append((i, j))

    return extracts_filtered

def write_cores_depre(outdir, win_extract, win_extract_2nd, extracts):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    count = 1

    for i, j in extracts:
        pr.writePDB(outdir + 'core_' + str(count) + '_' + win_extract[i][0].getTitle(), win_extract[i][0])
        pr.writePDB(outdir + 'core_' + str(count) + '_' + win_extract_2nd[j][0].getTitle(), win_extract_2nd[j][0])
        count += 1
    
def build_core_depre(outdir, target, query, query_2nd, rmsd_cut, rmsd_cut_2nd, distance_cut):
    
    print(query.getTitle())

    print(query_2nd.getTitle())
    
    win_extract = supperimpose_target_bb(query, target, rmsd_cut)

    win_extract_2nd = supperimpose_target_bb(query_2nd, target, rmsd_cut_2nd)

    print(len(win_extract))

    print(len(win_extract_2nd))

    if len(win_extract) == 0 or len(win_extract_2nd) ==0: 
        return len(win_extract), len(win_extract_2nd), None

    extracts = metal_distance_extract(target, win_extract, win_extract_2nd, distance_cut)
    if len(extracts) > 0:
        print(outdir)
        write_cores(outdir, win_extract, win_extract_2nd, extracts)

    return len(win_extract), len(win_extract_2nd), extracts



