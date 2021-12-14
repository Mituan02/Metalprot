import sys
sys.path.append(r'/wynton/home/degradolab/lonelu/GitHub_Design/Combs2')
import combs2

from metalprot.combs import search_cg_vdms as metalprot_scv
import prody as pr
import pandas as pd
import os
import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import lil_matrix

pd.set_option("display.max_columns", None)

input_dir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2/'

gly_pdb_name = '1cqs_allg.pdb'
ala_pdb_name = '1cqs_alla.pdb'

pdb_gly_path = input_dir + gly_pdb_name
pdb_ala_path = input_dir + ala_pdb_name

pdb_gly = pr.parsePDB(pdb_gly_path)
pdb_ala = pr.parsePDB(pdb_ala_path)

template = combs2.design.template.Template(pdb_gly)
template.set_alpha_hull(pdb_ala, alpha=9)

getResnums = pdb_ala.ca.getResnums() #grabs all residue numbers contained in the allALA backbone.
resnums = getResnums.tolist() #converts all residue numbers to a list format that can be read
resnums = [20, 28, 31, 32, 36, 38, 40, 47, 53, 56, 57, 60, 61, 66, 68, 86, 88, 99, 101, 103, 105, 113, 116, 118, 120]
chains = ['A'] * len(resnums)
segs = [''] * len(resnums)
segs_chains_resnums = zip(segs, chains, resnums)

# for x in segs_chains_resnums:
#     print(x)

cgs = ['coo', 'bb_cco', 'phenol']
cg_max_dists = dict(ph=0.9)
outpath = input_dir

combs2.design.functions.write_resfile(template, CGs=cgs,
                                        outpath=outpath,
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
                                        allowed_exposed='ACDEFGHIKLMNPQRSTVWY', allowed_intermed='ACDEFGHIKLMNPQRSTVWY',
                                        allowed_buried='ACDEFGHIKLMNPQRSTVWY',
                                        hb_only_residues='', all_contact_residues='')


path_to_resfile= input_dir + 'resfile.txt'
path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/vdMs/'
s = combs2.design._sample.Sample(**dict(path_to_resfile=outpath + 'resfile.txt',
                                        path_to_database=path_to_database))
s.read_resfile()


s.load_vdms(template, filter_by_phi_psi=False, run_parallel=True)

######################################################################

### Testing vdMs.
'''
dfa = s.cg_dict['coo']

x = dfa[['CG', 'rota', 'probe_name', 'seg_chain_resnum']].iloc[0]

v = dfa[(dfa['CG'] == x['CG']) & (dfa['rota'] == x['rota']) & (dfa['probe_name'] == x['probe_name']) & (dfa['seg_chain_resnum'] == x['seg_chain_resnum'])]
vdM_output_dir = input_dir + 'vdM/'
combs2.design.functions.print_dataframe(v, outpath=vdM_output_dir, tag='_ph')


vdM_output_dir_loada = input_dir + 'vdM_loaded_a/'
labels_loada = dfa[(dfa['chain'] == 'X') & (dfa['resname'] == 'GLN') & (dfa['name'] == 'CB') & (dfa['seg_chain_resnum'] == dfa['seg_chain_resnum'][59191])][['CG', 'rota', 'probe_name', 'seg_chain_resnum']]
for ind in range(labels_loada.shape[0]):
    x = labels_loada.iloc[ind]
    v = dfa[(dfa['CG'] == x['CG']) & (dfa['rota'] == x['rota']) & (dfa['probe_name'] == x['probe_name']) & (dfa['seg_chain_resnum'] == x['seg_chain_resnum'])]
    combs2.design.functions.print_dataframe(v, outpath=vdM_output_dir_loada, tag = '_' + str(ind))

vdM_output_dir_loada = input_dir + 'vdM_loaded_a200/'
labels_loada = dfa[(dfa['chain'] == 'X') & (dfa['resname'] == 'THR') & (dfa['name'] == 'CB') & (dfa['seg_chain_resnum'] == dfa['seg_chain_resnum'][1000000])][['CG', 'rota', 'probe_name', 'seg_chain_resnum']]
for ind in range(labels_loada.shape[0]):
    x = labels_loada.iloc[ind]
    v = dfa[(dfa['CG'] == x['CG']) & (dfa['rota'] == x['rota']) & (dfa['probe_name'] == x['probe_name']) & (dfa['seg_chain_resnum'] == x['seg_chain_resnum'])]
    combs2.design.functions.print_dataframe(v, outpath=vdM_output_dir_loada, tag = '_' + str(ind))

'''

######################################################################
ligands = [pr.parsePDB(input_dir + 'filtered_ligs/' + x) for x in os.listdir(workdir + 'filtered_ligs/')]

dfa = s.cg_dict['coo']

lgd_sel = ['C9', 'O3', 'O4']  # correpond to coo: ASP ['CG', 'OD1', 'OD2'] or GLU: ['CD', 'OE1', 'OE2']

ligand_coords = scv.get_ligand_coords(input_dir + 'filtered_ligs/', lgd_sel)

represent_name = 'OD2'

correspond_resname = 'ASP'

correspond_names = ['CG', 'OD1', 'OD2']

labels, vdm_coords = metalprot_scv.get_vdm_labels_coords(dfa, represent_name, correspond_resname, correspond_names)

'''

represent_name = 'OD2'

correspond_resname = 'ASP'

correspond_names = ['CG', 'OD1', 'OD2']

labels2, vdm_coords2 = get_vdm_labels_coords(dfa, represent_name, correspond_resname, correspond_names)

labels= pd.concat([labels, labelsa2])

vdm_coords = np.concatenate((vdm_coords, vdm_coords_a2))

'''

inds = metalprot_scv.get_nearest_vdms(vdm_coords, ligand_coords, radius = 2)

outdir = input_dir + 'vdM_candidate_a/'
prefix = 'a_'

write_vdms(outdir, inds, labels, dfa, prefix)

#################################################################


dfa = s.cg_dict['phenol']

lgd_sel = ['C3', 'O2'] # correpond to coo: ASP ['CG', 'OD1', 'OD2'] or GLU: ['CD', 'OE1', 'OE2']

ligand_coords = get_ligand_coords(input_dir + 'filtered_ligs/', lgd_sel)

represent_name = 'OH'

correspond_resname = 'TYR'

correspond_names = ['CZ', 'OH']

labels, vdm_coords = get_vdm_labels_coords(dfa, represent_name, correspond_resname, correspond_names)

inds = get_nearest_vdms(vdm_coords, ligand_coords, radius = 1)

outdir = input_dir + 'vdM_candidate_c/'
prefix = 'c_'

write_vdms(outdir, inds, labels, dfa, prefix)


##################################################################################################

input_dict = {
    (0, 0):{
        'cg' : 'coo',
        'lgd_sel' : ['C9', 'O3', 'O4'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CG', 'OD1', 'OD2']
    },
    (0, 1):{
        'cg' : 'coo',
        'lgd_sel' : ['C9', 'O3', 'O4'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CG', 'OD2', 'OD1']
    },    
    (0, 2):{
        'cg' : 'coo',
        'lgd_sel' : ['C9', 'O3', 'O4'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CG', 'OE1', 'OE2']
    },
    (0, 3):{
        'cg' : 'coo',
        'lgd_sel' : ['C9', 'O3', 'O4'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CG', 'OE2', 'OE1']
    },
    (1, 0):{
        'cg' : 'phenol',
        'lgd_sel' : ['C2', 'C3', 'O2'],
        'represent_name' : 'OH',
        'correspond_resname' : 'TYR',
        'correspond_names' : ['CE1', 'CZ', 'OH']
    },
    (1, 1):{
        'cg' : 'phenol',
        'lgd_sel' : ['C2', 'C3', 'O2'],
        'represent_name' : 'OH',
        'correspond_resname' : 'TYR',
        'correspond_names' : ['CE2', 'CZ', 'OH']
    },
    (2, 0):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C7', 'C8', 'O1'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'C', 'O']
    },
    (2, 1):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C7', 'C8', 'O1'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'C', 'O']
    }
}

ligands = [pr.parsePDB(input_dir + 'filtered_ligs/' + x) for x in os.listdir(input_dir + 'filtered_ligs/')]

labels_cgs = {}
df_cgs = {}
dist_ind_cgs = {}

cg_ids = input_dict.keys()
for cg_id in cg_ids:
    search_vdm(s.cg_dict, ligands, cg_id, input_dict, labels_cgs, df_cgs, dist_ind_cgs, rmsd = 0.7)

outdir = input_dir + 'output_07_all/'
os.makedirs(outdir, exist_ok= True)
CgCombInfoDict = construct_vdm_write(outdir, ligands, labels_cgs, df_cgs, dist_ind_cgs, clash_radius = 2.7)

write_summary(outdir, CgCombInfoDict, name = '_summary.tsv')


###################################################################################

ligands = [pr.parsePDB(input_dir + 'filtered_ligs/' + x) for x in os.listdir(input_dir + 'filtered_ligs/')]

labels_cgs = {}
df_cgs = {}
dist_ind_cgs = {}

cg_ids = input_dict.keys()
for cg_id in cg_ids:
    metalprot_scv.search_vdm(s.cg_dict, ligands, cg_id, input_dict, labels_cgs, df_cgs, dist_ind_cgs, rmsd = 0.7)

outdir = input_dir + 'output_07_all/'
os.makedirs(outdir, exist_ok= True)
CgCombInfoDict = metalprot_scv.construct_vdm_write(outdir, ligands, labels_cgs, df_cgs, dist_ind_cgs, clash_radius = 2.7)

metalprot_scv.write_summary(outdir, CgCombInfoDict, name = '_summary.tsv')