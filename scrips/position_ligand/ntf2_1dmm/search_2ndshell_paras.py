import os
import pickle
import numpy as np

def get_file_path(on_wynton):
    if on_wynton:
        with open('/wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/metalprot/constants/ideal_alanine_bb_only.pkl', 'rb') as f:
            ideal_alanine_bb_only = pickle.load(f)
        ideal_ala_coords = np.array(ideal_alanine_bb_only[['c_x', 'c_y', 'c_z']])

        path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/'

        #workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_sel/'
        workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_2rd_sel/'

        lig_path = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe/tts_fe_rdkit.pdb'

    else:
        with open('/mnt/e/GitHub_Design/Metalprot/metalprot/constants/ideal_alanine_bb_only.pkl', 'rb') as f:
            ideal_alanine_bb_only = pickle.load(f)
        ideal_ala_coords = np.array(ideal_alanine_bb_only[['c_x', 'c_y', 'c_z']])

        path_to_database='/mnt/e/DesignData/Combs/Combs2_database/'

        #workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/'
        workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/_rosetta_2ndRound/output_F55D_sel/'

        lig_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/tts_fe_rdkit.pdb'

    return workdir, path_to_database, ideal_ala_coords, lig_path


task_type = 'search_2ndshell'

### Target strcture 
#predefined_win_filters = [19]
predefined_win_filters = [15, 19, 27]
lig_cgs = [['hid', 'hie'], ['hid', 'hie'], ['coo']]
#lig_cgs = [['hid', 'hie']]

### Database para
use_enriched = True
use_abple=True


### Search ligands paramters.
predefined_resnums = [11, 12, 16, 20, 24, 28, 30, 31, 35, 37, 39, 46, 52, 55, 56, 60, 65, 67, 69, 83, 85, 87, 98, 100, 102, 104, 106, 112, 115, 117, 119]
#predefined_resnums = [12, 16, 20, 60, 65, 67, 69, 83, 85, 87, 102]


load_cg_aa_vdm_dict = {
    'coo': [('AGKNQRSTY', True, False)], # (aas, filter_hb, filter_cc)
    'hid': [('AGDEKNQRST', True, False)],
    'hie': [('AGDEKNQRST', True, False)],
}

rmsd = 1.0

#TO DO: rename the atoms
vdm_cg_aa_cc_dict = {
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