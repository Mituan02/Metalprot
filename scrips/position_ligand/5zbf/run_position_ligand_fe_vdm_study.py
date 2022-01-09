'''
Vdm Study metal + ligand bindin protein.
'''

from metalprot.search import search_selfcenter
from metalprot.basic import filter, prody_ext
import metalprot.basic.constant as constant
from metalprot.combs import position_ligand, search_ligand
import pickle
import time
import prody as pr
import os

import sys
sys.path.append(r'/wynton/home/degradolab/lonelu/GitHub_Design/Combs2')
import combs2

from metalprot.combs import search_cg_vdms as metalprot_scv
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
    
    sc = combs2.design._sample.Sample(**dict(path_to_resfile=path_to_resfile,
                                            path_to_database=path_to_database))
    sc.read_resfile()

    sc.load_vdms(template, filter_by_phi_psi=False, run_parallel=True)
    return sc


def write_all_vdms(outdir, ligand, sc, clash_radius, benchmark_filters):
    for seg_chain_resnum in benchmark_filters.keys():
        aa = benchmark_filters[seg_chain_resnum][0]
        for cgkey in benchmark_filters[seg_chain_resnum][1]:
            if cgkey not in sc.cg_dict.keys():
                print('The cg vdm is not loaded.')
            dfa = sc.cg_dict[cgkey]
            print('_'.join([str(z) for z in seg_chain_resnum]) + '_' + aa)
            print(cgkey)
            print(dfa.shape)
            vdM_output_dir_loada = outdir + '_'.join([str(z) for z in seg_chain_resnum]) + '_' + cgkey + '/'
            extract_vdm_by_aa(dfa, seg_chain_resnum, aa, ligand, vdM_output_dir_loada, clash_radius)
    return


def extract_vdm_by_aa(dfa, seg_chain_resnum, aa, ligand, vdM_output_dir_loada, clash_radius):
    labels_loada = dfa[(dfa['chain'] == 'X') & (dfa['resname'] == aa) & (dfa['name'] == 'CB') & (dfa['seg_chain_resnum'] == seg_chain_resnum)][['CG', 'rota', 'probe_name', 'seg_chain_resnum']]
    print(labels_loada.shape)
    for ind in range(labels_loada.shape[0]):
        x = labels_loada.iloc[ind]
        v = dfa[(dfa['CG'] == x['CG']) & (dfa['rota'] == x['rota']) & (dfa['probe_name'] == x['probe_name']) & (dfa['seg_chain_resnum'] == x['seg_chain_resnum'])]
        if metalprot_scv.vdm_ligand_clash(v, ligand, clash_radius):
            #print('clash ' + str(ind))
            continue
        combs2.design.functions.print_dataframe(v, outpath=vdM_output_dir_loada, tag = '_' + str(ind))


workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/5zbf_fe/'

outdir = workdir + 'vdm_study4/'
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

resnums = [30, 40, 41, 43, 100]

path_to_resfile = workdir + 'resfile_benchmark_bbindep.txt'
sc = load_vdms(outdir, pdb_gly, pdb_ala, resnums, use_exist_resfile= True, path_to_resfile = path_to_resfile)

clash_radius = 2.7

# benchmark_filters = {
#     ('', 'A', 30): 'GLN',
#     ('', 'A', 40): 'LEU',
#     ('', 'A', 41): 'GLU',
#     ('', 'A', 43): 'HIS',
#     ('', 'A', 100): 'LYS'
# }


benchmark_filters = {
    ('', 'A', 30): ('GLN', ['bb_cco', 'coo']),
    ('', 'A', 40): ('LEU', ['phenol']),
    ('', 'A', 41): ('GLU', ['phenol']),
    ('', 'A', 43): ('HIS', ['phenol', 'ph']),
    ('', 'A', 100): ('LYS', ['bb_cco', 'coo'])
}

# benchmark_filters = {
#     ('', 'A', 43): ('HIS', ['phenol', 'ph'])
# }


write_all_vdms(outdir, ligand, sc, clash_radius, benchmark_filters)