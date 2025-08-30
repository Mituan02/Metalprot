import pytest
import os
import shutil
import prody as pr
import itertools
import numpy as np

from metalprot.apps.core import Core
from metalprot.apps.search_struct import Query, Search_struct
from metalprot import search_struct, extract_vdm, ligand_database

def test_overlap():
    win_seen = dict()

    win_seen[0] = [1, 2]
    win_seen[1] = [2, 3]
    win_seen[2] = [3, 1]

    test1 = Search_struct.overlap(0, 2, win_seen, 3)
    assert test1

    test2 = Search_struct.bivalence_extract_filter(0, 2, win_seen, 4)
    assert test2 == False
    

def test_supperimpose_target_bb_bivalence():
    #workdir = '/mnt/e/GitHub_Design/Metalprot/tests/test_data/'
    workdir = os.path.dirname(os.path.realpath(__file__)) + '/test_data/'

    outdir = workdir + 'output_search_3vdm/'

    target_path = workdir + '5od1_zn.pdb'

    target = pr.parsePDB(target_path)

    query_path = workdir + 'm5_bb_cluster_0_mem_324_centroid_5keb_ZN_1_AAdAAMetal_HIS_HIS_mem0.pdb'

    query = Query(pr.parsePDB(query_path), 0, 0, 0)

    dist_array, id_array = search_struct.get_contact_map(target)

    assert len(dist_array) == 4371

    win_filter = [30, 31, 32, 33, 34, 35, 36, 37, 38, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 94]

    new_querys = search_struct.supperimpose_target_bb_bivalence(query, target, dist_array, id_array, win_filter, tolerance = 0.5, rmsd_cut = 0.5)

    assert len(new_querys) == 3

