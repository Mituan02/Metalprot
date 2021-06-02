import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import ligand_database as ldb


def extract_vdm(cores, outdir, AA, key, extention, k, key_out):
    for c in cores:   
        if extention:
            c.generate_AA_ext_Metal(AA, extention, key)
        else:
            c.generate_AA_Metal(AA, key)

        if k:
            c.generate_AA_kMetal(AA, k, key, key_out)
            c.write_vdM(outdir, key = key_out)
        else:
            c.write_vdM(outdir, key = key)


### set up parameters

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb/"

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 
align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE or ion'

cores = ldb.load_cores(workdir + '2_seq_cores_reps/')

        
### Align aa core w/o sc------------------------------------------------------------------------------------

# AAMetal_HIS
extract_vdm(cores, workdir + "M1_AAMetal_HIS_reps/", AA = 'HIS', key = 'AAMetal_HIS', k = None, key_out=None)

_pdbs = ldb.get_all_pbd_prody(workdir + "M1_AAMetal_HIS_reps/")

ldb.run_cluster(_pdbs, workdir, 'M1_AAMetal_HIS_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 5, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm1_')

ldb.run_cluster(_pdbs, workdir, 'M1_AAMetalSc_HIS_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 11, align_sel = 'heavy', min_cluster_size = 2, tag = 'm1-1_') 

#ldb.run_cluster(_pdbs, workdir, 'M1_AAMetalTest_HIS_cluster05/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 5, align_sel = 'name C CA N CB ZN', min_cluster_size = 2, tag = 'm1-1_') 


# AA5Metal_HIS
extract_vdm(cores, workdir + "M1_AA5Metal_HIS_reps/", AA = 'HIS', key = 'AAMetal_HIS', k = 5, key_out='AA5Metal_HIS')

_pdbs2 = ldb.get_all_pbd_prody(workdir + "M1_AA5Metal_HIS_reps/")

ldb.run_cluster(_pdbs2, workdir, 'M1_AA5Metal_HIS_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 9, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm1_')

ldb.run_cluster(_pdbs2, workdir, 'M1_AA5MetalSc_HIS_cluster02/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 15, align_sel = 'heavy', min_cluster_size = 2, tag = 'm1-1_') 


### Align his core +- 1 AA------------------------------------------------------------------------------------ 

# AAExtMetal_HIS_ext1
extract_vdm(cores, workdir + "M3-1_AAextMetal_HIS_reps/", AA = 'HIS', key = 'AAextMetal_HIS_ext3', extention=1, k = None, key_out=None)

_pdbs = ldb.get_all_pbd_prody(workdir + 'M3-1_AAextMetal_HIS_reps/')

ldb.run_cluster(_pdbs, workdir, 'M3-1_AAextMetal_HIS_cluster05/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 13, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm3-1_')

# AAExt5Metal_HIS_ext1
extract_vdm(cores, workdir + "M3-1_AAext5Metal_HIS_reps/", AA = 'HIS', key = 'AAextMetal_HIS_ext3', extention=1, k = 5, key_out='AAext5Metal_HIS_ext3')

_pdbs = ldb.get_all_pbd_prody(workdir + 'M3-1_AAext5Metal_HIS_reps/')

ldb.run_cluster(_pdbs, workdir, 'M3-1_AAext5Metal_HIS_cluster05/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 17, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm3-1_')


