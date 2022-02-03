'''
For the generation of single vdMers. 
This is the third step. To run this, run the 'database_extract_vdm_from_core' first. And then run '.../search_selfcenter/database_selfcenter_pickle.py'
'''

import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
from metalprot.database import database_cluster as ldb_clu


'''
python /mnt/e/GitHub_Design/Metalprot/scrips/database_generate/run_cluster_single.py
'''

def cluster_single(workdir, outdir, is_self_center, metal_sel, align_sel, len_sel):

    #aas = ['HIS', 'GLU', 'ASP', 'CYS']
    aas = ['HIS', 'GLU', 'ASP']
    for AA in aas:
        #Change len_sel according to aa sidechain number. His:6, Glu: 5, Asp: 4, Cys: 2
        if AA == 'HIS':
            len_sel_sc = len_sel + 6 
        elif AA == 'GLU':
            len_sel_sc = len_sel + 5
        elif AA == 'ASP':
            len_sel_sc = len_sel + 4 
        elif AA == 'CYS':
            len_sel_sc = len_sel + 2 

        ### Align core with phipsi---------------- 

        if AA == 'ASP' or AA == 'GLU':
            _pdbs = database_extract.get_all_pbd_prody(workdir + 'AAMetalPhiPsi_' + AA + '_reps_shift/')
        else:
            _pdbs = database_extract.get_all_pbd_prody(workdir + 'AAMetalPhiPsi_' + AA + '_reps/')
        print(len(_pdbs))
        
        ldb_clu.run_cluster(_pdbs, outdir, 'AAMetalPhiPsi_' + AA + '_cluster05/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel_sc, align_sel = align_sel, min_cluster_size = 0, tag = 'AAMetalPhiPsi_' + AA + '_', is_self_center=is_self_center)


### set up parameters

workdir = "/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20211013/20211013_vdm_reps/"
workdir = '/mnt/e/DesignData/ligands/all/20220116_FE_MN_CO/20220116_vdm_reps/'
workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/20220128_1stshell/'
is_self_center = True


metal_sel = 'name NI MN ZN CO CU FE' 
#align_sel = 'heavy'
#len_sel = 9
align_sel = 'heavy and not name NI MN ZN CO CU FE'
len_sel = 8

outdir = '/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20211013/20211209_selfcenter_nometal/'
outdir = '/mnt/e/DesignData/ligands/all/20220116_FE_MN_CO/20220116_selfcenter/'
outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/20220128_1stshell/20220128_selfcenter/'
os.makedirs(outdir, exist_ok = True)

cluster_single(workdir, outdir, is_self_center, metal_sel, align_sel, len_sel)
        
### Align aa core w/o sc------------------------------------------------------------------------------------
'''
# AAMetal_HIS
extract_single_vdm(cores, workdir + 'M1_AAMetal_' + AA + '_reps/', AA = AA, key = 'AAMetal_HIS', basic = True)

_pdbs = ldb.get_all_pbd_prody(workdir + 'M1_AAMetal_' + AA + '_reps/')

ldb.run_cluster(_pdbs, workdir, 'M1_AAMetal_' + AA + '_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = len_sel, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm1_')

ldb.run_cluster(_pdbs, workdir, 'M1-1_AAMetalSc_' + AA + '_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = len_sel_sc, align_sel = 'heavy', min_cluster_size = 2, tag = 'm1-1_') 

#ldb.run_cluster(_pdbs, workdir, 'M1_AAMetalTest_HIS_cluster05/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 5, align_sel = 'name C CA N CB ZN', min_cluster_size = 2, tag = 'm1-1_') 
'''

'''
# AA5Metal_HIS
extract_single_vdm(cores, workdir + "M1_AA5Metal_HIS_reps/", AA = AA, key = 'AAMetal_HIS', basic = False, n = 5, key_out='AA5Metal_HIS')

_pdbs2 = ldb.get_all_pbd_prody(workdir + "M1_AA5Metal_HIS_reps/")

ldb.run_cluster(_pdbs2, workdir, 'M1_AA5Metal_HIS_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 9, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm1_')

ldb.run_cluster(_pdbs2, workdir, 'M1_AA5MetalSc_HIS_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 15, align_sel = 'heavy', min_cluster_size = 2, tag = 'm1-1_') 
'''

### Align his core +- 1 AA------------------------------------------------------------------------------------ 
'''
# AAExtMetal_HIS_ext1
extract_single_vdm(cores, workdir + 'M3-1_AAextMetal_' + AA + '_reps/', AA = AA, key = 'AAextMetal_' + AA + '_ext3', basic = False, extention=1, n = None, key_out=None)

_pdbs = ldb.get_all_pbd_prody(workdir + 'M3-1_AAextMetal_' + AA + '_reps/')

ldb.run_cluster(_pdbs, workdir, 'M3-1_AAextMetal_' + AA + '_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 13, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm3-1_')
'''


'''
# AAExt5Metal_HIS_ext1
extract_single_vdm(cores, workdir + "M3-1_AAext5Metal_' + AA + '_reps/", AA = AA, key = 'AAextMetal_' + AA + '_ext3', extention=1, n = 5, key_out='AAext5Metal_' + AA + '_ext3')

_pdbs = ldb.get_all_pbd_prody(workdir + 'M3-1_AAext5Metal_' + AA + '_reps/')

ldb.run_cluster(_pdbs, workdir, 'M3-1_AAext5Metal_' + AA + '_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 17, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm3-1_')
'''

