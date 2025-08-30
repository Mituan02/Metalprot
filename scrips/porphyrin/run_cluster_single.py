from genericpath import exists
import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import ligand_database as ldb


def extract_single_vdm(cores, outdir, AA, key, basic = True, extention = None, n = None, key_out = None, phipsi = False):
    for c in cores:   
        if basic:
            c.generate_AA_Metal(AA, key)
        elif extention:
            c.generate_AA_ext_Metal(AA, extention, key)          
        elif n:
            c.generate_AA_kMetal(AA, n, key, key_out) 
        elif phipsi:
            c.generate_AA_phipsi_Metal(AA, key)
        c.write_vdM(outdir, key = key)
    return

### set up parameters

workdir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/'

metal_sel = 'name NI MN ZN CO CU MG FE' 


cores = ldb.load_cores(workdir + 'contact_core_reps_aligned/')


AA = 'MET'
align_sel= 'name NI MN ZN CO CU MG FE or resname MET and not name O'
len_sel = 4

#Change len_sel according to aa sidechain number. His:6, Glu: 5, Asp: 4, Cys: 2
if AA == 'HIS':
    len_sel_sc = len_sel + 6 
elif AA == 'GLU':
    len_sel_sc = len_sel + 5
elif AA == 'ASP':
    len_sel_sc = len_sel + 4 
elif AA == 'CYS':
    len_sel_sc = len_sel + 2 
elif AA == 'MET':
    len_sel_sc = len_sel + 4

        
### Align aa core w/o sc------------------------------------------------------------------------------------

# AAMetal_HIS
extract_single_vdm(cores, workdir + 'M1_AAMetal_' + AA + '_reps/', AA = AA, key = 'AAMetal'+ AA, basic = True)

_pdbs = ldb.get_all_pbd_prody(workdir + 'M1_AAMetal_' + AA + '_reps/')

ldb.run_cluster(_pdbs, workdir, 'M1-1_AAMetalSc_' + AA + '_cluster05/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel_sc, align_sel = align_sel, min_cluster_size = 0, tag = 'm1-1_') 

