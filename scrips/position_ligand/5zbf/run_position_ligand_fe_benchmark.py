'''
Benchmark metal + ligand bindin protein.
'''

from metalprot.search import search_selfcenter
from metalprot.basic import filter, prody_ext
import metalprot.basic.constant as constant
from metalprot.combs import __search_ligand, position_ligand
import pickle
import time
import prody as pr
import os

import sys
sys.path.append(r'/wynton/home/degradolab/lonelu/GitHub_Design/Combs2')
import combs2

from metalprot.combs import __search_cg_vdms as metalprot_scv
import pandas as pd
import shutil

def load_vdms(input_dir, pdb_gly, pdb_ala, resnums, use_exist_resfile = False, path_to_resfile = '', path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/vdMs/'):
    '''
    Load vdMs on the current protein target.
    '''
    os.makedirs(input_dir, exist_ok = True)

    template = combs2.design.template.Template(pdb_gly)
    template.set_alpha_hull(pdb_ala, alpha=9)

    if not use_exist_resfile:
        getResnums = pdb_ala.ca.getResnums() #grabs all residue numbers contained in the allALA backbone.
        if len(resnums) <= 0:
            resnums = getResnums.tolist() #converts all residue numbers to a list format that can be read
        chains = ['A'] * len(resnums)
        segs = [''] * len(resnums)
        segs_chains_resnums = zip(segs, chains, resnums)

        # for x in segs_chains_resnums:
        #     print(x)
        cgs = ['coo', 'bb_cco', 'phenol', 'ph']

        cg_max_dists = dict(ph=0.9)
        combs2.design.functions.write_resfile(template, CGs=cgs,
                                                outpath=input_dir,
                                                filename='resfile', tag='',
                                                resindices=None, segs_chains_resnums=segs_chains_resnums,
                                                pikaa_dict=None, bb_dep=1,
                                                use_enriched_vdMs=False, CA_burial_distance=None, exclude_exposed=False,
                                                exclude_intermed=False,
                                                exclude_buried=False, top_exposed=None, top_intermed=None, top_buried=None,
                                                alpha_hull_radius=10,
                                                use_propensities=True,
                                                propensity_threshold=0.9, use_abple=True, use_dssp=False,
                                                path_to_pdb_for_dssp=None,
                                                allowed_exposed='LYQEKH', allowed_intermed='LYQEKH',
                                                allowed_buried='LYQEKH',
                                                hb_only_residues='', all_contact_residues='')

        path_to_resfile= input_dir + 'resfile.txt'
    print(path_to_resfile)
    sc = combs2.design._sample.Sample(**dict(path_to_resfile=path_to_resfile,
                                            path_to_database=path_to_database))
    sc.read_resfile()

    sc.load_vdms(template, filter_by_phi_psi=False, run_parallel=True)
    return sc

def run_combs(outdir, ligands, sc, input_dict, rmsd, clash_radius, benchmark_filters):
    '''
    Search vdms which are close to the ligands cgs.
    '''
    os.makedirs(outdir, exist_ok= True)

    labels_cgs = {}
    df_cgs = {}
    dist_ind_cgs = {}

    cg_ids = input_dict.keys()
    for cg_id in cg_ids:
        metalprot_scv.search_vdm(sc.cg_dict, ligands, cg_id, input_dict, labels_cgs, df_cgs, dist_ind_cgs, rmsd = rmsd)

    CgCombInfoDict = metalprot_scv.construct_vdm_write(outdir, ligands, labels_cgs, input_dict, df_cgs, dist_ind_cgs, clash_radius = clash_radius, benchmark_filters = benchmark_filters)

    metalprot_scv.write_summary(outdir, CgCombInfoDict, name = '_summary.tsv')
    return



workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/5zbf_fe/'

outdir = workdir + 'benchmark1-/'
os.makedirs(outdir, exist_ok = True)

target_path = workdir + '5zbf.pdb'
target = pr.parsePDB(target_path)

lig_path = workdir + 'eno_fe2.pdb'
ligand = pr.parsePDB(lig_path)

ideal_geo_o_path = workdir + 'fe_geo_o.pdb'



pdb_gly = prody_ext.target_to_all_gly_ala(target, 'all_gly.pdb', [], aa = 'GLY')
pdb_ala = prody_ext.target_to_all_gly_ala(target, 'all_ala.pdb', [], aa = 'ALA')

pr.writePDB(outdir + 'all_gly.pdb', pdb_gly)
pr.writePDB(outdir + 'all_ala.pdb', pdb_ala)

pdb_gly = pr.parsePDB(outdir + 'all_gly.pdb')
pdb_ala = pr.parsePDB(outdir + 'all_ala.pdb')


input_dict = {
    (0, 0):{
        'cg' : 'coo',
        'lgd_sel' : ['C1', 'O1', 'O2'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : [ 'CG', 'OD1', 'OD2']
    },
    (0, 1):{
        'cg' : 'coo',
        'lgd_sel' : [ 'C1', 'O1', 'O2'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : [ 'CG', 'OD2', 'OD1']
    },    
    (0, 2):{
        'cg' : 'coo',
        'lgd_sel' : [ 'C1', 'O1', 'O2'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : [ 'CG', 'OE1', 'OE2']
    },
    (0, 3):{
        'cg' : 'coo',
        'lgd_sel' : ['C1', 'O1', 'O2'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CG', 'OE2', 'OE1']
    },
    (1, 0):{
        'cg' : 'phenol',
        'lgd_sel' : ['C6', 'C7', 'O3'],
        'represent_name' : 'OH',
        'correspond_resname' : 'TYR',
        'correspond_names' : ['CE1', 'CZ', 'OH']
    },
    (1, 1):{
        'cg' : 'phenol',
        'lgd_sel' : ['C6', 'C7', 'O3'],
        'represent_name' : 'OH',
        'correspond_resname' : 'TYR',
        'correspond_names' : ['CE2', 'CZ', 'OH']
    },
    (2, 0):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C3', 'C2', 'O4'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (2, 1):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C3', 'C2', 'O4'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (2, 2):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C3', 'C2', 'O4'],
        'represent_name' : 'O',
        'correspond_resname' : 'PRO',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (3, 0):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C2', 'C1', 'O1'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (3, 1):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C2', 'C1', 'O1'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (3, 2):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C2', 'C1', 'O1'],
        'represent_name' : 'O',
        'correspond_resname' : 'PRO',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (4, 0):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C2', 'C1', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (4, 1):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C2', 'C1', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (4, 2):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C2', 'C1', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'PRO',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (5, 0):{
        'cg' : 'ph',
        'lgd_sel' : ['C4', 'C6', 'C8'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CG', 'CE1', 'CE2']
    },
    (5, 1):{
        'cg' : 'ph',
        'lgd_sel' : ['C4', 'C6', 'C8'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CG', 'CE2', 'CE1']
    },
    (5, 2):{
        'cg' : 'ph',
        'lgd_sel' : ['C4', 'C6', 'C8'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE1', 'CG', 'CE2']
    },
    (5, 3):{
        'cg' : 'ph',
        'lgd_sel' : ['C4', 'C6', 'C8'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE1', 'CE2', 'CG']
    },
    (5, 4):{
        'cg' : 'ph',
        'lgd_sel' : ['C4', 'C6', 'C8'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE2', 'CG', 'CE1']
    },
    (5, 5):{
        'cg' : 'ph',
        'lgd_sel' : ['C4', 'C6', 'C8'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE2', 'CE1', 'CG']
    }
}


resnums = [30, 40, 41, 43, 100]
path_to_resfile = workdir + 'resfile_benchmark_bbindep.txt'
print(path_to_resfile)
print(pdb_gly)
print(pdb_ala)
sc = load_vdms(outdir, pdb_gly, pdb_ala, resnums, use_exist_resfile= True, path_to_resfile = path_to_resfile)

rmsd = 1.2
clash_radius = 2.7

benchmark_filters = {
    ('', 'A', 30): 'GLN',
    ('', 'A', 40): 'LEU',
    ('', 'A', 41): 'GLU',
    ('', 'A', 43): 'HIS',
    ('', 'A', 100): 'LYS'
}

run_combs(outdir + 'combs_out/', [ligand], sc, input_dict, rmsd, clash_radius, benchmark_filters)

