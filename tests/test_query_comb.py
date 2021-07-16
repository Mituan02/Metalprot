from metalprot.apps.core import Core
import pytest
import os
import shutil
import prody as pr

from metalprot import ligand_database
from metalprot import core
from metalprot.apps.search_struct import Query, Comb

def test_contact():
    workdir = os.path.dirname(os.path.realpath(__file__)) + '/test_data/'
    workdir = '/mnt/e/GitHub_Design/Metalprot/tests/test_data/'

    pdb1 = pr.parsePDB(workdir + '5_1_3.37_m1-1_cluster_18_mem_20_5od1_ZN_1_AAMetal_HIS_mem0.pdb')
    pdb2 = pr.parsePDB(workdir + '5_2_3.37_m1-1_cluster_20_mem_20_4k1r_ZN_5_AAMetal_HIS_mem0.pdb')
    pdb3 = pr.parsePDB(workdir + '5_3_3.37_m1-1_cluster_22_mem_31_5od1_ZN_1_AAMetal_HIS_mem1.pdb')

    querys = []
    querys.append(Query(pdb1, 0, 0, 0, win = [34]))
    querys.append(Query(pdb2, 0, 0, 0, win = [64]))
    querys.append(Query(pdb3, 0, 0, 0, win = [60]))

    comb = Comb(querys)
    comb.calc_pair_geometry()

    print(comb.pair_dists)
    print(comb.pair_angles)


