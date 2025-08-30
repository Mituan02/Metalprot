import pytest
import os
import shutil
import prody as pr
import itertools
import numpy as np

from metalprot.apps.core import Core
from metalprot.apps.search_struct import Query, constructy_pseudo_2ndshellVdm, convert_query_2ndshellVdm
from metalprot import search_struct, extract_vdm, ligand_database


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

