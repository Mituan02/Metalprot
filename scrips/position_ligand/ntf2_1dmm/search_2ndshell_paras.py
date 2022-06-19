import os
import pickle
import numpy as np

def get_file_path(on_wynton):
    if on_wynton:

        path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/'

        #workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_sel/'
        workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_2rd_sel/'

        lig_path = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe/tts_fe_rdkit.pdb'

    else:

        path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

        #workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/'
        #workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/_rosetta_2ndRound/output_F55D_sel/'
        workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta_16-20-28/_rosetta_tts_r2_876/output_F55D_sel/'
        lig_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/tts_fe_rdkit.pdb'

    return workdir, path_to_database, lig_path


task_type = 'search_2ndshell'

### Target strcture 
#predefined_win_filters = [19]
predefined_win_filters = [('A', 15), ('A', 19), ('A', 27)]
lig_cgs = [['hid', 'hie'], ['hid', 'hie'], ['coo']]
#lig_cgs = [['hid', 'hie']]

### Database para
use_enriched = True
use_abple=True


### Search ligands paramters.
predefined_resnums = [('A', 11), ('A', 12), ('A', 16), ('A', 20), ('A', 24), ('A', 28), 
        ('A', 30), ('A', 31), ('A', 35),('A', 37), ('A', 39), ('A', 46), ('A', 52), ('A', 55), 
        ('A', 56), ('A', 60), ('A', 65), ('A', 67), ('A', 69), ('A', 83), ('A', 85), ('A', 87), 
        ('A', 98), ('A', 100), ('A', 102), ('A', 104), ('A', 106), ('A', 112), ('A', 115), ('A', 117), ('A', 119)]
#predefined_resnums = [12, 16, 20, 60, 65, 67, 69, 83, 85, 87, 102]


load_cg_aa_vdm_dict = {
    'coo': [('AGKNQRSTY', True, False)], # (aas, filter_hb, filter_cc)
    'hid': [('AGDEKNQRST', True, False)],
    'hie': [('AGDEKNQRST', True, False)],
}

rmsd = 1.0

#TO DO: rename the atoms
vdm_cg_aa_atommap_dict = {
    ('coo_0'):{
        'cg' : 'coo',
        'lgd_sel' : ['CB', 'CG', 'OD1', 'OD2'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CB', 'CG', 'OD1', 'OD2']
    },
    ('coo_1'):{
        'cg' : 'coo',
        'lgd_sel' : ['CB', 'CG', 'OD1', 'OD2'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CG', 'CD', 'OE1', 'OE2']
    },
    ('hid_0'):{
        'cg' : 'hid',
        'lgd_sel' : ['ND1', 'CE1', 'CG', 'CD2'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['ND1', 'CE1', 'CG', 'CD2']
    },
    ('hie_0'):{
        'cg' : 'hie',
        'lgd_sel' : ['NE2', 'CE1', 'CD2', 'CG'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['NE2', 'CE1', 'CD2', 'CG']
    }
}