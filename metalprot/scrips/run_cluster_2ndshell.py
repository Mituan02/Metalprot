import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import ligand_database as ldb

### set up parameters

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb/20210608/"

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 
align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE or ion'

cores = ldb.load_cores(workdir + '_Seq_cores_with_2ndshell_reps/')

for c in cores:   
    c.generate_AA_2ndShell_Metal(key = 'AA2sMetal-HIS', filter_AA =True, AA = 'HIS')
    c.write_vdM(workdir + 'M7_AA2sMetal-HIS_reps/',  key = 'AA2sMetal-HIS')


_pdbs = ldb.get_all_pbd_prody(workdir + 'M7_AA2sMetal-HIS_reps/')
ldb.run_cluster(_pdbs, workdir, 'M7_AA2sMetal-HIS_clusters/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 9, align_sel = align_sel_backbone, min_cluster_size = 2, tag = 'm7_')


