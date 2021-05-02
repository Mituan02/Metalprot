import os
import sys
import prody as pr
sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import ligand_database as ldb

workdir = "/mnt/e/DesignData/ligands/CU_NI/"

# According to the prody atom flag (http://prody.csb.pitt.edu/manual/reference/atomic/flags.html#flags), NI MN ZN CO CU MG FE CA are not all flag as ion. 
# Note that calcium is read as CA, which is same as alpha carbon in prody selection. 
# So to select calcium, we need to add ion before it.

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE or ion'

aa = 'resname HIS' #Change len_sel according to aa in the following code.

aa_name = 'HIS' 

len_sel = 5

len_sel_sc = len_sel + 6 #Change len_sel according to aa sidechain number. His:6, Glu: 5, Asp: 4, Cyc: 2

len_sel_ps = len_sel + 4 

len_sel_ps_sc = len_sel_ps + 6 #Change len_sel according to aa sidechain number. His:6, Glu: 5, Asp: 4, Cyc: 2


'''
# Extract aa core (metal with contact aa) with aas from cores_seq for manual check.

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

cores = extract_all_core(pdbs, metal_sel)

writepdb(cores, workdir + '3_aa_cores_reps/')

'''

pdbs = ldb.get_all_pbd_prody(workdir + '2_seq_cores_reps/')

# Align aa core

aa_cores = ldb.extract_all_core_aa(pdbs, metal_sel, aa = aa)
ldb.writepdb(aa_cores, workdir + 'M1_' + aa_name + '_cores_reps/')

_pdbs = ldb.get_all_pbd_prody(workdir + 'M1_' + aa_name + '_cores_reps/')
ldb.run_cluster(_pdbs, workdir, 'M1_' + aa_name + '_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm1_')

# Align aa core + sidechain

ldb.run_cluster(_pdbs, workdir, 'M1-1_' + aa_name + '_sc_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel_sc, align_sel = 'heavy', min_cluster_size = 2, tag = 'm1-1_') 

# Align aa core + phipsi.

aa_ps_cores = ldb.extract_all_core_aa(pdbs, metal_sel, aa = aa, consider_phipsi = True)
ldb.writepdb(aa_ps_cores, workdir + 'M2_' + aa_name + '_ps_cores_reps/')

_pdbs = ldb.get_all_pbd_prody(workdir + 'M2_' + aa_name + '_ps_cores_reps/')
ldb.run_cluster(_pdbs, workdir, 'M2_' + aa_name + '_ps_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel_ps, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm2_')

# Align aa core + sidechain + phipsi

ldb.run_cluster(_pdbs, workdir, 'M2-2_' + aa_name + '_sc_ps_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel_ps_sc, align_sel = 'heavy', min_cluster_size = 2, tag = 'm2-1_') 


# Align his core + 1 AA 

aa_2aa_cores = ldb.extract_all_core_aa(pdbs, metal_sel, aa = aa, extention= 1)
ldb.writepdb(aa_2aa_cores, workdir + 'M3-1_' + aa_name + '_1aa_cores_reps/')

_pdbs = ldb.get_all_pbd_prody(workdir + 'M3-1_' + aa_name + '_1aa_cores_reps/')
ldb.run_cluster(_pdbs, workdir, 'M3-1_' + aa_name + '_1aa_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 13, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm3-1_')

# Align his core + 2 AA 

aa_2aa_cores = ldb.extract_all_core_aa(pdbs, metal_sel, aa = aa, extention= 2)
ldb.writepdb(aa_2aa_cores, workdir + 'M3-2_' + aa_name + '_2aa_cores_reps/')

_pdbs = ldb.get_all_pbd_prody(workdir + 'M3-2_' + aa_name + '_2aa_cores_reps/')
ldb.run_cluster(_pdbs, workdir, 'M3-2_' + aa_name + '_2aa_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 21, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm3-2_')

# Align his core + 3 AA 

aa_2aa_cores = ldb.extract_all_core_aa(pdbs, metal_sel, aa = aa, extention= 3)
ldb.writepdb(aa_2aa_cores, workdir + 'M3-3_' + aa_name + '_3aa_cores_reps/')

_pdbs = ldb.get_all_pbd_prody(workdir + 'M3-3_' + aa_name + '_3aa_cores_reps/')
ldb.run_cluster(_pdbs, workdir, 'M3-3_' + aa_name + '_3aa_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 29, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm3-3_')
