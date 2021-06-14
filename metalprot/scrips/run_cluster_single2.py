import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import ligand_database as ldb


def extract_vdm(cores, outdir, AA, key, extention, n, key_out):
    for c in cores:   
        if extention:
            c.generate_AA_ext_Metal(AA, extention, key)
        else:
            c.generate_AA_Metal(AA, key)

        if n:
            c.generate_AA_kMetal(AA, n, key, key_out)
            c.write_vdM(outdir, key = key_out)
        else:
            c.write_vdM(outdir, key = key)


### set up parameters

workdir = "/mnt/e/DesignData/ligands/FE_MN_CO/20210610/"

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 
align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE or ion'

cores = ldb.load_cores(workdir + '_Seq_cores_reps/')


AA = 'CYS'

len_sel = 5

#Change len_sel according to aa sidechain number. His:6, Glu: 5, Asp: 4, Cys: 2
if AA == 'HIS':
    len_sel_sc = len_sel + 6 
elif AA == 'GLU':
    len_sel_sc = len_sel + 5
elif AA == 'ASP':
    len_sel_sc = len_sel + 4 
elif AA == 'CYS':
    len_sel_sc = len_sel + 2 

        
### Align aa core w/o sc------------------------------------------------------------------------------------

# AAMetal_HIS
extract_vdm(cores, workdir + 'M1_AAMetal_' + AA + '_reps/', AA = AA, key = 'AAMetal_HIS', extention = None, n = None, key_out=None)

_pdbs = ldb.get_all_pbd_prody(workdir + 'M1_AAMetal_' + AA + '_reps/')

ldb.run_cluster(_pdbs, workdir, 'M1_AAMetal_' + AA + '_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = len_sel, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm1_')

ldb.run_cluster(_pdbs, workdir, 'M1-1_AAMetalSc_' + AA + '_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = len_sel_sc, align_sel = 'heavy', min_cluster_size = 2, tag = 'm1-1_') 

#ldb.run_cluster(_pdbs, workdir, 'M1_AAMetalTest_HIS_cluster05/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 5, align_sel = 'name C CA N CB ZN', min_cluster_size = 2, tag = 'm1-1_') 

'''
# AA5Metal_HIS
extract_vdm(cores, workdir + "M1_AA5Metal_HIS_reps/", AA = AA, key = 'AAMetal_HIS', n = 5, key_out='AA5Metal_HIS')

_pdbs2 = ldb.get_all_pbd_prody(workdir + "M1_AA5Metal_HIS_reps/")

ldb.run_cluster(_pdbs2, workdir, 'M1_AA5Metal_HIS_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 9, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm1_')

ldb.run_cluster(_pdbs2, workdir, 'M1_AA5MetalSc_HIS_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 15, align_sel = 'heavy', min_cluster_size = 2, tag = 'm1-1_') 
'''

### Align his core +- 1 AA------------------------------------------------------------------------------------ 

# AAExtMetal_HIS_ext1
extract_vdm(cores, workdir + 'M3-1_AAextMetal_' + AA + '_reps/', AA = AA, key = 'AAextMetal_' + AA + '_ext3', extention=1, n = None, key_out=None)

_pdbs = ldb.get_all_pbd_prody(workdir + 'M3-1_AAextMetal_' + AA + '_reps/')

ldb.run_cluster(_pdbs, workdir, 'M3-1_AAextMetal_' + AA + '_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 13, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm3-1_')

'''
# AAExt5Metal_HIS_ext1
extract_vdm(cores, workdir + "M3-1_AAext5Metal_' + AA + '_reps/", AA = AA, key = 'AAextMetal_' + AA + '_ext3', extention=1, n = 5, key_out='AAext5Metal_' + AA + '_ext3')

_pdbs = ldb.get_all_pbd_prody(workdir + 'M3-1_AAext5Metal_' + AA + '_reps/')

ldb.run_cluster(_pdbs, workdir, 'M3-1_AAext5Metal_' + AA + '_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 17, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm3-1_')
'''
