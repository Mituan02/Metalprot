'''
2ndshell Hbond interaction is extracted from Probe program.
The script cluster the 2ndshell vdms. To run the script, run database_prep_2ndshell_probe.py first. Then run '.../search_2ndshell_pickle.py' to get the pickle form.
'''

import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
from metalprot.database import database_cluster as ldb_clu
from metalprot.basic import vdmer_2ndshell 

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/"

outdir = workdir + '_Seq_core_2ndshell_date_reps_probe2ndshell/'

metal_sel = 'name NI MN ZN CO CU MG FE' 
align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE'


for aa in ['HIS', 'GLU', 'ASP']:
#for aa in ['CYS']:
    _key = aa

    #Organise the vdm in right prody order.
    _pdbs = database_extract.get_all_pbd_prody(outdir + _key + '/')
    outdir_reps_x = outdir + _key + '_repsx/'
    os.makedirs(outdir_reps_x, exist_ok=True)
    for _pdb in _pdbs:
        _pdb_x = vdmer_2ndshell.organize_2ndshellVdm_probe(_pdb)
        # if not _pdb_x:
        #     continue
        pr.writePDB(outdir_reps_x + _pdb.getTitle(), _pdb_x)

    _pdbsx = database_extract.get_all_pbd_prody(outdir_reps_x)
    select_poses =[] 
    ldb_clu.run_cluster(_pdbsx, outdir, _key + '_clustersx/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 9, align_sel = align_sel_backbone, min_cluster_size = 2, tag = _key + '_', is_self_center = True, select_poses=select_poses)
