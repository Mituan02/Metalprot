from metalprot.apps.core import Core
import pytest
import os
import shutil
import prody as pr

from metalprot import ligand_database as ldb
from metalprot import core

def test_core_kMetal():

    workdir = os.path.dirname(os.path.realpath(__file__)) + '/test_data/'

    pdb_prody = pr.parsePDB(workdir + '5od1_zn.pdb')

    core = Core(pdb_prody)

    core.generate_AA_Metal(AA = 'HIS')

    core.generate_AA_kMetal(AA = 'HIS', k = 5, key = 'AA_Metal', key_out = 'AA_kMetal')

    core.write_vdM(workdir + 'output/', key = 'AA_kMetal')

    pdb_kMetal = pr.parsePDB(workdir + 'output/' + '5od1_zn_AA_kMetal_mem0.pdb')

    assert len(pdb_kMetal.select('name ZN')) == 5

    #shutil.rmtree(workdir + 'output/')