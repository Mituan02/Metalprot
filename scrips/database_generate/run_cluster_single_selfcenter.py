import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract as ldb
from metalprot.database import database_cluster as ldb_clu
from metalprot.database import database_vdMAtomOrder

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/database_generate/run_cluster_single_selfcenter.py
'''

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

workdir = "/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20211008/20211009_selfcenter/"


for folder in os.listdir(workdir):

    AA = folder.split('_')[1]
    len_sel = 9

    #Change len_sel according to aa sidechain number. His:6, Glu: 5, Asp: 4, Cys: 2
    if AA == 'HIS':
        len_sel_sc = len_sel + 6 
    elif AA == 'GLU':
        len_sel_sc = len_sel + 5
    elif AA == 'ASP':
        len_sel_sc = len_sel + 4 
    elif AA == 'CYS':
        len_sel_sc = len_sel + 2 

    print(workdir + folder)
    _pdbs = ldb.get_all_pbd_prody(workdir + folder + '/')

    ldb_clu.run_cluster(_pdbs, workdir, 'AAMetalPhiPsi_' + AA + '_cluster05/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel_sc, align_sel = 'heavy', min_cluster_size = 0, tag = 'AAMetalPhiPsi_' + AA + '_')

