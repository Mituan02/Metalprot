import pytest
import os

from metalprot import ligand_database as ldb


def test_2ndshell():
    #workdir = '/mnt/e/GitHub_Design/Metalprot/tests/data/'
    workdir = os.path.dirname(os.path.realpath(__file__)) + '/data/'
    print(workdir)
    pdbs = ldb.get_all_pbd_prody(workdir)
    assert '5od1' in pdbs[0].getTitle()

    pdb_prody = pdbs[0]
    metal_sel = 'name ZN'

    #step test into 'get_metal_core_seq'
    nis = pdb_prody.select(metal_sel)
    ni = nis[0]
    ni_index = ni.getIndex()

    all_near = pdb_prody.select('protein and within 2.83 of index ' + str(ni_index))
    assert len(all_near.select('nitrogen or oxygen or sulfur')) >= 3

    inds = all_near.select('nitrogen or oxygen or sulfur').getResindices()
    assert len(inds) == 3

    #test 'extend_res_indices'
    ext_inds = ldb.extend_res_indices(inds, pdb_prody, extend=4)
    assert len(ext_inds) == 22

    #test 'get_2ndshell_indices'
    _2ndshell_resindices = ldb.get_2ndshell_indices(inds, pdb_prody, ni_index)
    assert _2ndshell_resindices[0] == 57

    #test 'get_metal_core_seq_2ndshell'
    core_2ndshell = ldb.get_metal_core_seq_2ndshell(pdb_prody, metal_sel)
    #ldb.writepdb(core_2ndshell, workdir + 'output/')


 