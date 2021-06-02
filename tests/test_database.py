import pytest
import os
import shutil
import prody as pr

from metalprot import ligand_database as ldb


def test_2ndshell():
    #workdir = '/mnt/e/GitHub_Design/Metalprot/tests/data/'
    workdir = os.path.dirname(os.path.realpath(__file__)) + '/test_data/'

    pdb_prody = pr.parsePDB(workdir + '5od1_zn.pdb')
    assert '5od1' in pdb_prody.getTitle()
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
    ldb.writepdb(core_2ndshell, workdir + 'output/')

    pdb_core = pr.parsePDB(workdir + 'output/' + '5od1_zn_ZN_1.pdb')
    #test 'extract_all_core_aa' with extract2ndshell=True
    aa_cores = ldb.extract_all_core_aa([pdb_core], metal_sel, aa = 'resname HIS', extract2ndshell=True)
    ldb.writepdb(aa_cores, workdir + 'output/his/')
    assert len(aa_cores[1].select('bb')) == 8 #There are two amino acid

    #shutil.rmtree(workdir + 'output/')

 