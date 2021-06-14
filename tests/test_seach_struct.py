import pytest
import os
import shutil
import prody as pr
import itertools
import numpy as np

from metalprot.apps.core import Core
from metalprot.apps.search_struct import Query, constructy_pseudo_2ndshellVdm, convert_query_2ndshellVdm
from metalprot import search_struct, extract_vdm, ligand_database

def test_generate_ind_combination_listoflist():
    _listoflist = [[1, 2], [3, 4], [5, 6]]
    all_inds = search_struct.generate_ind_combination_listoflist(_listoflist)
    assert len(all_inds) == 8
    assert all_inds[1] == (0, 0, 1)

def test_get_combs_from_pair_extract():
    inds_len = 3
    extracts = [[(1, 2), (2, 4)], [(1, 3), (1, 5)], [(2, 3), (2, 5), (4, 6)]]
    #Test development   
    all_inds = search_struct.generate_ind_combination_listoflist(extracts)
    inds = all_inds[0]
    extract = [extracts[ind] for ind in inds]
    assert extract == [(1, 2), (1, 3), (2, 3)]

    combs = search_struct.get_combs_from_pair_extract(inds_len, extracts)
    assert combs == [[1, 2, 3], [1, 2, 5]]

def test_2ndshellVdm():

    #workdir = os.path.dirname(os.path.realpath(__file__)) + '/test_data/'
    workdir = '/mnt/e/GitHub_Design/Metalprot/tests/test_data/'
    outdir = workdir + 'output/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    target = pr.parsePDB(workdir + '5od1_zn.pdb')

    query_pdb = pr.parsePDB(workdir + '5_3_3.37_m1-1_cluster_22_mem_31_5od1_ZN_1_AAMetal_HIS_mem1.pdb')
    query = Query(query_pdb, 0, 0, 0)

    #The protein has a second shell vdm '57(GLN) -> 60(HIS) -> Zn'
    ags = constructy_pseudo_2ndshellVdm(target, query, 60)

    #The query_2nd_pdb is used as an centroid from a library.
    query_2nd_pdb = pr.parsePDB(workdir + '5od1_ZN_1_2nd__2ndShell_mem0.pdb')

    query_2nd = Query(query_2nd_pdb, 0, 0, 0)

    query_2nd_ag = convert_query_2ndshellVdm(query_2nd)

    candidates = []
    for ag in ags:
        candidate = search_struct.supperimpose_2ndshell(ag, query_2nd_ag, query_2nd, rmsd_cut = 0.5)
        if candidate:
            candidates.append(candidate)

    assert len(candidates) == 1

    '''
    #Develop supperimpose_2ndshell

    transform = pr.calcTransformation(query_2nd_ag, ag)
    transform.apply(query_2nd_ag)
    transform.apply(query_2nd.query)
    rmsd = pr.calcRMSD(ag, query_2nd_ag)
    if rmsd < 0.5:
        candidate = Query(query_2nd.query.copy(),  query_2nd.score, query_2nd.clu_num, query_2nd.clu_total_num, query_2nd.win, query_2nd.path)

    pr.writePDB(outdir + '_2nd_shell_search_' + candidate.query.getTitle() , candidates[0].query)
    '''        


def test_get_pair_extracts():
    #workdir = '/mnt/e/GitHub_Design/Metalprot/tests/test_data/'
    workdir = os.path.dirname(os.path.realpath(__file__)) + '/test_data/'

    queryss = []
    query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb/'

    #Get query pdbs 
    querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M1-1_AAMetalSc_HIS_cluster02'], score_cut = 1, clu_num_cut = 100)
    #assert len(querys) == 14
    queryss.append([querys[18], querys[20], querys[22]])
    queryss.append([querys[18], querys[20], querys[22]])
    queryss.append([querys[18], querys[20], querys[22]])
    # queryss.append(querys)
    # queryss.append(querys)
    # queryss.append(querys)

    contact_querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M8_AtomContact4_clusters'], score_cut = 0, clu_num_cut = 2)

    outdir = workdir + 'output_search_3vdm/'

    target_path = workdir + '5od1_zn.pdb'

    rmsd_cuts = [0.5, 0.5, 0.5]

    dist_cuts = [1, 1, 1]

    num_iter = 3

    clash_query_query = 2.3

    clash_query_target = 2.3

    use_sep_aas = [False, False, False]

    tolerance = 0.5

    ss = search_struct.Search_struct(target_path, outdir, queryss, rmsd_cuts, dist_cuts, num_iter, clash_query_query, clash_query_target, use_sep_aas, tolerance, contact_querys = contact_querys)

    #ss.run_search_struct()

    ss.generate_cquerys()

    ss.run_iter_search_structure()



    '''
    #Develop get_pair_extracts()
    all_inds = search_struct.generate_ind_combination_listoflist(ss.queryss)
    assert len(all_inds) == 2744

    inds = all_inds[0]
    xys = list(itertools.combinations(range(len(inds)), 2))

    x, y = xys[2]
    cquerys_a = ss.cquerysss[x][inds[x]]
    cquerys_b = ss.cquerysss[y][inds[y]]
    #assert len(cquerys_a) == len(cquerys_b)

    extracts = ss._metal_distance_extract(cquerys_a, cquerys_b, dist_cut = 1)
    extracts_filtered = ss.distance_extracts_filter(cquerys_a, cquerys_b,extracts)

    extracts = ss.get_pair_extracts(inds, dist_cut = 1)

    assert ss.pair_extracts[x][y]
    '''

    
    '''
    #Develop run_search_struct()  
    comb_inds = []
    for inds in all_inds:
        extracts = ss.get_pair_extracts(inds, dist_cut = 1)
        if extracts and len(extracts)>0:
            combs = search_struct.get_combs_from_pair_extract(ss.num_iter, extracts)
            for comb in combs:
                comb_inds.append((inds, comb))

    ss.build_combs(comb_inds)
    ss.write_combs()
    ss.write_comb_info()
    '''

    ss.cquerysss = []
    ss.generate_cvdms(ss.target)

    # if not os.path.exists(workdir + 'test/'):
    #     os.mkdir(workdir + 'test/')
    # for t in ss.cquerysss[0][0]:
    #     pr.writePDB(workdir + 'test/' + t.query.getTitle(), t.query)
    import itertools
    xys = list(itertools.combinations(range(ss.num_iter), 2))
    ss.pair_extracts = [[0]* ss.num_iter]* ss.num_iter
    ss.pair_extracts = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    for x, y in xys:
        ss.pair_extracts[x][y] = dict()

    all_inds = [[i]*ss.num_iter for i in range(len(ss.cquerysss[0]))]
    print(len(all_inds))

    comb_inds = []
    for inds in all_inds:
        extracts = ss.get_pair_extracts(inds, dist_cuts = [0.3]*ss.num_iter)
        if extracts and len(extracts)>0:
            combs = search_struct.get_combs_from_pair_extract(ss.num_iter, extracts)
            for comb in combs:
                comb_inds.append((inds, comb))


    ss.combs = []
    ss.build_combs(comb_inds)

    ss.write_combs( outpath= '/mem_combs/')

    ss.write_comb_info( filename= '/_summary_mem.txt')