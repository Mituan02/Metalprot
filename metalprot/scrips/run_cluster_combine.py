import os
import sys
import prody as pr
sys.path.append(r'/mnt/e/GitHub_Design/MetalDesign')
from metalprot import ligand_database as ldb


workdir = "/mnt/e/DesignData/ligands/ZN_rcsb/"

# According to the prody atom flag (http://prody.csb.pitt.edu/manual/reference/atomic/flags.html#flags), NI MN ZN CO CU MG FE CA are not all flag as ion. 
# Note that calcium is read as CA, which is same as alpha carbon in prody selection. 
# So to select calcium, we need to add ion before it.

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE or ion'


pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')


# Align 2 his core 
def extract_2aa_core(workdir, pdbs, metal_sel, aa_name, extention_out, aas):
    aa_2aa_cores = ldb.extract_all_core_aa(pdbs, metal_sel, aa = None, extention= 3, extention_out = extention_out, extract2aa=True, aas = aas)
    ldb.writepdb(aa_2aa_cores, workdir + 'M4_' + aa_name + '_cores_reps/')

    _pdbs = ldb.get_all_pbd_prody(workdir + 'M4_' + aa_name + '_cores_reps/')
    for i in range(2, 10):
        subworkdir = 'M4_' + aa_name + '_cores_cluster_' + str(i) + '/'
        ldb.run_cluster(_pdbs, workdir, subworkdir, rmsd = 0.5, metal_sel = metal_sel, len_sel = i*4 + 1, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm1_')


# Align 2 his core 
all_aas = ['HIS', 'ASP', 'GLU', 'CYS']
extention_out = 0
for i in range(len(all_aas)):
    for j in range(i, len(all_aas)):
        aas = [all_aas[i], all_aas[j]]
        aa_name = '_'.join(aas)
        extract_2aa_core(workdir, pdbs, metal_sel, aa_name, extention_out, aas)


