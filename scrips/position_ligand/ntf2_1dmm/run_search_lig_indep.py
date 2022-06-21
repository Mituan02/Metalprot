'''
The script is to design ligand-metal binding enzyme using the vdM idea. 
Specially, here is to search the combs vdM library without using Combs2.
'''

import os
import sys
import pdb
import prody as pr
import pickle
import pandas as pd
import numpy as np
from metalprot.basic import constant, utils
from metalprot.combs import search_lig_indep, search_lig_indep_wrap, search_lig_indep_2ndshell 
from metalprot.combs import position_ligand
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import lil_matrix
import gc
import shutil
import datetime
import importlib.machinery

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/run_search_lig_indep.py 0 /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/search_lig_paras_rdkit.py 7

python run_search_lig_indep.py 0 search_lig_paras_eval.py 7
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/run_search_lig_indep.py 0 /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/search_2ndshell_paras.py 5
'''


### Run Search ligands.

def run_all(file, workdir, path_to_database, lig_path, para):
    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 

    target_path = workdir + file
    outdir = workdir + para.task_type  + '_result/output_' + '_'+ file + '_' + time_tag + '/'
    os.makedirs(outdir, exist_ok=True)

    target, chidres2ind = search_lig_indep.prepare_rosetta_target(outdir, target_path, para.predefined_win_filters)
    
    if para.task_type == 'search_2ndshell':
        ligs = search_lig_indep_2ndshell.extract_2ndshell_cgs(outdir, target, para.predefined_win_filters)
    else:
        ligs = search_lig_indep.prepare_ligs(outdir, target, lig_path, para)

    print('Filtered Ligs: ' + str(len(ligs)))
    if len(ligs) <= 0:
        return

    outdir_uni = outdir + 'vdms_output_uni/'
    os.makedirs(outdir_uni, exist_ok=True)

    outdir_all = outdir + 'vdms_output_all/'
    os.makedirs(outdir_all, exist_ok=True)
    print('number of ligs: {}'.format(len(ligs)))
    for i in range(len(ligs)):
        lig = ligs[i]
        lig_vdm_dict = search_lig_indep_wrap.run_search(target, lig, path_to_database, constant.ideal_ala_coords, para, para.lig_cgs[i])
        print('lig_vdm_dict size {}'.format(len(list(lig_vdm_dict.keys()))))
        search_lig_indep.write_vdm(outdir, outdir_all, outdir_uni, target, lig_vdm_dict, i, para)

    return


def main():
    #path = '/mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/search_lig_paras.py'
    on_wynton = bool(int(sys.argv[1]))
    path = sys.argv[2]
    para = importlib.machinery.SourceFileLoader('para', path).load_module()
    print(path)
    print('Task: ' + para.task_type)
    workdir, path_to_database, lig_path = para.get_file_path(on_wynton)
    print('on_wynton: ' + str(on_wynton))

    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    ind = int(sys.argv[3]) -1
    if ind > len(pdb_files) -1:
        return
    print(pdb_files[ind])
    
    run_all(pdb_files[ind], workdir,  path_to_database, lig_path, para)
    return


if __name__=='__main__':
    main()

