import os
import sys
from metalprot.combs import search_lig_indep_wrap 
import datetime

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/run_search_2ndshell.py 1
''' 

class Para:
    #>>> On wynton path
    #path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/'
    #workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_2rd_sel/'

    #>>> On local path
    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'
    workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta_16-20-28/_rosetta_tts_r2_876/output_F55D_sel/'
    lig_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/tts_fe_rdkit.pdb'

    #>>> There are three type of tasks: 'search_unknow', 'search_eval', or 'search_2ndshell'
    task_type = 'search_2ndshell'

    ### Target strcture 
    predefined_win_filters = [('A', 15), ('A', 19), ('A', 27)]
    lig_cg_2ndshell = [['hid', 'hie'], ['hid', 'hie'], ['coo']]

    #>>> vdM database filter.
    use_enriched = True
    use_abple=True

    #>>> Search ligands only with vdMs on the predefined chidres. If task_type is 'search_2ndshell', set it to None.
    predefined_chidres = None

    rmsd = 0.75

    #TO DO: rename the atoms
    vdm_cg_aa_atommap_dict = {
        ('coo_0'):{
            'cg' : 'coo',
            'lgd_sel' : ['CB', 'CG', 'OD1', 'OD2'],
            'represent_name' : 'OD2',
            'correspond_resname' : 'ASP',
            'correspond_names' : ['CB', 'CG', 'OD1', 'OD2'],
            'aas' : 'AGKNQRSTY',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('coo_1'):{
            'cg' : 'coo',
            'lgd_sel' : ['CB', 'CG', 'OD1', 'OD2'],
            'represent_name' : 'OE2',
            'correspond_resname' : 'GLU',
            'correspond_names' : ['CG', 'CD', 'OE1', 'OE2'],
            'aas' : 'AGKNQRSTY',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('hid_0'):{
            'cg' : 'hid',
            'lgd_sel' : ['ND1', 'CE1', 'CG', 'CD2'],
            'represent_name' : 'ND1',
            'correspond_resname' : 'HIS',
            'correspond_names' : ['ND1', 'CE1', 'CG', 'CD2'],
            'aas' : 'AGDEKNQRST',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('hie_0'):{
            'cg' : 'hie',
            'lgd_sel' : ['NE2', 'CE1', 'CD2', 'CG'],
            'represent_name' : 'NE2',
            'correspond_resname' : 'HIS',
            'correspond_names' : ['NE2', 'CE1', 'CD2', 'CG'],
            'aas' : 'AGDEKNQRST',
            'filter_hb' : True,
            'filter_cc' : False
        }
    }


def run_local():
    para = Para()
    print('Task: ' + para.task_type)
    workdir, path_to_database, lig_path = para.workdir, para.path_to_database, None

    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    ind = int(sys.argv[1]) -1
    if ind > len(pdb_files) -1:
        return

    print(pdb_files[ind])

    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 
    outdir = workdir + para.task_type  + '_result/output_' + '_'+ pdb_files[ind] + '_' + time_tag + '/'
    os.makedirs(outdir, exist_ok=True)
    
    search_lig_indep_wrap.run_search_2ndshell(pdb_files[ind], workdir, outdir, path_to_database, lig_path, para)
    return


if __name__=='__main__':
    run_local()