import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract as ldb
from metalprot.database import database_cluster as ldb_clu

### set up parameters

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211018/"
outdir = workdir + '20211026_2ndshell_selfcenter/'

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 
align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE or ion'

cores = ldb.load_cores(workdir + '_Seq_core_2ndshell_date_reps/')

for aa in ['HIS', 'GLU', 'ASP']:
    _key = 'AA2ndS-' + aa
    for c in cores:   
        c.generate_AA_2ndShell_Metal(key = _key, filter_AA =True, AA = aa)
        c.write_vdM(outdir + _key + '_reps/',  key = _key)


    _pdbs = ldb.get_all_pbd_prody(outdir + _key + '_reps/')
    ldb_clu.run_cluster(_pdbs, outdir, _key + '_clusters/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 9, align_sel = align_sel_backbone, min_cluster_size = 2, tag = _key + '_', is_self_center = True)


