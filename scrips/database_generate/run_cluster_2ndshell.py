import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
from metalprot.database import database_cluster as ldb_clu
from metalprot.basic import vdmer_2ndshell 

### set up parameters
workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/"


metal_sel = 'name NI MN ZN CO CU MG FE' 
align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE'

cores = database_extract.load_cores(workdir + '_Seq_core_2ndshell_date_reps/')


### extract reps and cluster
outdir = workdir + '20220116_selfcenter_bb2ndshell_connect/'
os.makedirs(outdir, exist_ok=True)
only_bb_2ndshell = False
generate_sse = False #The generate_sse can only be True when only_bb_2ndshell = True.

for aa in ['HIS', 'GLU', 'ASP']:
    _key = 'AA2ndS-' + aa
    for c in cores:   
        if not generate_sse:
            c.generate_AA_2ndShell_Metal(key = _key, filter_AA =True, AA = aa, only_bb_2ndshell = only_bb_2ndshell)
        else:
            c.generate_AA_2ndShell_connect_Metal(key = _key, filter_AA =True, AA = aa, only_bb_2ndshell = only_bb_2ndshell, extend = 4 ,_2nd_extend = 4, generate_sse = generate_sse)
        c.write_vdM(outdir  + _key + '_reps/',  key = _key)

    #Organise the vdm in right prody order.
    _pdbs = database_extract.get_all_pbd_prody(outdir + _key + '_reps/')
    outdir_reps_x = outdir + _key + '_repsx/'
    os.makedirs(outdir_reps_x, exist_ok=True)
    for _pdb in _pdbs:
        _pdb_x = vdmer_2ndshell.organize_2ndshellVdm(_pdb)
        if not _pdb_x:
            continue
        pr.writePDB(outdir_reps_x + _pdb.getTitle(), _pdb_x)

    _pdbsx = database_extract.get_all_pbd_prody(outdir_reps_x)
    select_poses =[] 
    if generate_sse:
        for _pdb in _pdbsx:
            inds = _pdb.select('protein and name CA').getResindices()
            ind_metal = _pdb.select(metal_sel).getResindices()[0]
            select_poses.append((inds[0], inds[-1], ind_metal))
    ldb_clu.run_cluster(_pdbsx, outdir, _key + '_clustersx/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 9, align_sel = align_sel_backbone, min_cluster_size = 2, tag = _key + '_', is_self_center = True, select_poses=select_poses)


### Cluster for probe extracted 2ndshell.

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


#for aa in ['HIS', 'GLU', 'ASP']:
for aa in ['CYS']:
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
