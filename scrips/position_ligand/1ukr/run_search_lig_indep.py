'''
The script is to design ligand-metal binding enzyme using the vdM database. 
The ligand is superimposed on the metal with all possible conformers. 
And the ligs (different conformers) are then searched against the vdM database.
'''

import os
import sys
from metalprot.combs import search_lig_indep, search_lig_indep_wrap 

import prody as pr
import datetime

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/1ukr/run_search_lig_indep.py local 1

'''

#>>> Parameters
class ParaLig:
    #>>> Ligand generation
    ro1 = ['C8', 'C7']
    rest1 = ['H7', 'C9', 'O1', 'O3', 'O2', 'FE1']
    ro2 = ['C7', 'C6']
    rest2 = ['H5', 'H6', 'H7', 'C8', 'C9', 'O1', 'O3', 'O2', 'FE1']
    rot_degree = [8, 8]
    interMolClashSets = [(['O1'], ['C1', 'C5']), (['O2'], ['C1', 'C5', 'C6'])]

    #>>> Ligand superimpose
    lig_connects = [['O1','O3', 'FE1'], ['O3','O1', 'FE1']]
    #geo_sel = 'chid X and name O1 O3 FE1'
    geo_sel = 'chid X and name O2 O3 FE1'
    clash_dist = 3.0
    lig_name = 'TTS'


class Para:

    ### TaskType
    task_type = 'search_unknow' #>>> 'search_unknow', 'search_eval', 'search_exists'

    #>>> iurk
    predefined_win_filters = [('A',68), ('A',77), ('A',79)]
    predefined_resnums = [4, 6, 8, 9, 24, 34, 35, 37, 62, 64, 66, 68, 70, 77, 79, 81, 87, 109, 111, 113, 116, 119, 122, 125, 127, 129, 161, 162, 166, 168, 170]
    
    predefined_chidres = [('A', x) for x in predefined_resnums]

    #>>> Database para
    use_enriched = True
    use_abple=True
    rmsd = 0.75

    vdm_cg_aa_atommap_dict = {
        ('coo_0'):{
            'cg' : 'coo',
            'lgd_sel' : ['C8', 'C9', 'O3', 'O2'],
            'represent_name' : 'OD2',
            'correspond_resname' : 'ASP',
            'correspond_names' : ['CB', 'CG', 'OD1', 'OD2'],
            'aas' : 'GAWSTYNQDEKRH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('coo_1'):{
            'cg' : 'coo',
            'lgd_sel' : ['C8', 'C9', 'O3', 'O2'],
            'represent_name' : 'OD2',
            'correspond_resname' : 'ASP',
            'correspond_names' : ['CB', 'CG', 'OD2', 'OD1'],
            'aas' : 'GAWSTYNQDEKRH',
            'filter_hb' : True,
            'filter_cc' : False
        },    
        ('coo_2'):{
            'cg' : 'coo',
            'lgd_sel' : ['C8', 'C9', 'O3', 'O2'],
            'represent_name' : 'OE2',
            'correspond_resname' : 'GLU',
            'correspond_names' : ['CG', 'CD', 'OE1', 'OE2'],
            'aas' : 'GAWSTYNQDEKRH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('coo_3'):{
            'cg' : 'coo',
            'lgd_sel' : ['C8', 'C9', 'O3', 'O2'],
            'represent_name' : 'OE2',
            'correspond_resname' : 'GLU',
            'correspond_names' : ['CG', 'CD', 'OE2', 'OE1'],
            'aas' : 'GAWSTYNQDEKRH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('phenol_0'):{
            'cg' : 'phenol',
            'lgd_sel' : ['C2', 'C3', 'C4', 'O4'],
            'represent_name' : 'OH',
            'correspond_resname' : 'TYR',
            'correspond_names' : ['CE1', 'CZ', 'CE2', 'OH'],
            'aas' : 'GAWSTYNQDEKRH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('phenol_1'):{
            'cg' : 'phenol',
            'lgd_sel' : ['C2', 'C3', 'C4', 'O4'],
            'represent_name' : 'OH',
            'correspond_resname' : 'TYR',
            'correspond_names' : ['CE2', 'CZ', 'CE1', 'OH'],
            'aas' : 'GAWSTYNQDEKRH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_0'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C8', 'C9', 'O2'],
            'represent_name' : 'O',
            'correspond_resname' : 'GLY',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'GAWSTYNQDEKRH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_1'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C8', 'C9', 'O2'],
            'represent_name' : 'O',
            'correspond_resname' : 'ALA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'GAWSTYNQDEKRH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_2'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C8', 'C9', 'O2'],
            'represent_name' : 'O',
            'correspond_resname' : 'PRO',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'GAWSTYNQDEKRH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('ph_0'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CG', 'CD1', 'CD2', 'CZ'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_1'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CG', 'CD2', 'CD1', 'CZ'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_2'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CD1', 'CG', 'CE1', 'CE2'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_3'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CD1', 'CE1', 'CG', 'CE2'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_4'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CE1', 'CD1', 'CZ', 'CD2'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_5'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CE1', 'CZ', 'CD1', 'CD2'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_6'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CZ', 'CE1', 'CE2', 'CG'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_7'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CZ', 'CE2', 'CE1', 'CG'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_8'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CE2', 'CZ', 'CD2', 'CD1'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_9'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CE2', 'CD2', 'CZ', 'CD1'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_10'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CD2', 'CE2', 'CG', 'CE1'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_11'):{
            'cg' : 'ph',
            'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
            'represent_name' : 'CG',
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CD2', 'CG', 'CE2', 'CE1'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        }
    }


#>>> Run Search ligands.

def run_all(file, workdir, path_to_database, lig_path, para, para_lig):
    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 

    target_path = workdir + file
    outdir = workdir + para.task_type  + '_'+ file.split('.')[0] + '_' + time_tag + '/'
    os.makedirs(outdir, exist_ok=True)

    target = pr.parsePDB(target_path)

    ligs = search_lig_indep_wrap.prepare_ligs(outdir, target, task_type = para.task_type, lig_path = lig_path, para_lig = para_lig)
    
    target, chidres2ind = search_lig_indep.prepare_rosetta_target(outdir, target_path, para.predefined_win_filters)
    
    print('Filtered Ligs: ' + str(len(ligs)))
    if len(ligs) <= 0:
        return

    outdir_uni = outdir + 'vdms_output_uni/'
    os.makedirs(outdir_uni, exist_ok=True)

    outdir_all = outdir + 'vdms_output_all/'
    os.makedirs(outdir_all, exist_ok=True)
    print('number of ligs: {}'.format(len(ligs)))

    lig_vdm_dict = search_lig_indep_wrap.run_search(target, ligs, path_to_database, para, para.predefined_chidres, lig_cg_2ndshell = None)
    
    print('lig_vdm_dict size {}'.format(len(list(lig_vdm_dict.keys()))))
    search_lig_indep_wrap.write_vdm(outdir, outdir_all, outdir_uni, target, lig_vdm_dict,  para)

    return


def run_local():
    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'
    lig_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/tts_fe_rdkit.pdb'

    #>>>
    workdir = '/mnt/e/DesignData/Metalloenzyme/1ukr/68-77-79_DHH/'
    

    para = Para()
    para_lig = ParaLig()
    print('Task: ' + para.task_type)

    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    ind = int(sys.argv[2]) -1
    if ind > len(pdb_files) -1:
        return
    print(pdb_files[ind])
    
    run_all(pdb_files[ind], workdir,  path_to_database, lig_path, para, para_lig)
    return


def run_wynton():
    path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/'
    workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_sel/'
    #workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_2rd_sel/'
    lig_path = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe/tts_fe_rdkit.pdb'
    
    para = Para()
    para_lig = ParaLig()
    print('Task: ' + para.task_type)

    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    ind = int(sys.argv[2]) -1
    if ind > len(pdb_files) -1:
        return
    print(pdb_files[ind])
    
    run_all(pdb_files[ind], workdir,  path_to_database, lig_path, para, para_lig)
    return

if __name__=='__main__':
    if sys.argv[1] == 'wynton': 
        run_wynton()
    elif sys.argv[1] == 'local':
        run_local()

