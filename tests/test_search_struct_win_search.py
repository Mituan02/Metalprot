import pytest
import os
import shutil
import prody as pr
import itertools
import numpy as np

from metalprot.apps.core import Core
from metalprot.apps.search_struct import Query
from metalprot import search_struct, extract_vdm, ligand_database


def test_generate_win_query_dict():
    inds_len = 3
    extracts = [[(1, 2), (2, 4)], [(1, 3), (1, 5)], [(2, 3), (2, 5), (4, 6)]]
    #Test development   
    all_inds = search_struct.generate_ind_combination_listoflist(extracts)
    inds = all_inds[0]
    extract = [extracts[ind] for ind in inds]
    assert extract == [(1, 2), (1, 3), (2, 3)]

    combs = search_struct.get_combs_from_pair_extract(inds_len, extracts)
    assert combs == [[1, 2, 3], [1, 2, 5]]


def test_win_based_search():
    workdir = '/mnt/e/GitHub_Design/Metalprot/tests/test_data/'
        #workdir = os.path.dirname(os.path.realpath(__file__)) + '/test_data/'

    queryss = []
    query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb/'

    #Get query pdbs 
    querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M1-1_AAMetalSc_HIS_cluster02'], score_cut = 1, clu_num_cut = 100)
    #assert len(querys) == 14
    queryss.append(querys)

    contact_querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M8_AtomContact4_clusters'], score_cut = 1, clu_num_cut = 2)

    outdir = workdir + 'output_search_3vdm/'

    target_path = workdir + '5od1_zn.pdb'

    rmsd_cuts = [0.5, 0.5, 0.5]

    dist_cuts = [1.5, 1.5, 1.5]

    num_iter = 3

    clash_query_query = 2.3

    clash_query_target = 2.3

    use_sep_aas = [False, False, False]

    tolerance = 0.5

    ss = search_struct.Search_struct(target_path, outdir, queryss, rmsd_cuts, dist_cuts, num_iter, clash_query_query, clash_query_target, use_sep_aas, tolerance, contact_querys = contact_querys)

    ss.generate_win_query_dict()

    ss.iter_all_wins() 

    ss.build_win_combs()

    ss.write_combs()

    ss.write_comb_info()