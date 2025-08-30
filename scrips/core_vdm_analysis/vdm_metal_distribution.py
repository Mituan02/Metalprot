import os
import prody as pr

from metalprot.database import database_extract as ldb
from metalprot.database import database_evaluate as ev

### 

workdir = '/mnt/e/DesignData/ligands/TYR_HIS_vdM_sc_only/'
centroid_pdb_name = 'TYR_HIS_155_C=2.350_41_4.pdb'
outdir = workdir + 'test_out/'
if not os.path.exists(outdir):
    os.mkdir(outdir)

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/AAMetalPhiPsi_HIS_cluster05/0/'
centroid_pdb_name = 'AAMetalPhiPsi_HIS_cluster_0_mem_236_centroid_2018_5yt6_ZN_3_mem1.pdb'
outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/AAMetalPhiPsi_HIS_cluster05/'

### Start

pdbs = ldb.get_all_pbd_prody(workdir)

centroid_pdb = pr.parsePDB(workdir + centroid_pdb_name)

metal_coords, dists = ev.ev_atom_distribution(pdbs, centroid_pdb, atom_sel = 'name ND1 or name NE2', align_sel = 'heavy')

ev.plt_dist_his(dists, outdir + 'dists_origin.png')

ev.plt_3d_points(metal_coords, outdir + 'metal3d_origin.png')

with open(outdir + 'cluster_0_dists.txt', 'w') as f:
    for d in dists:
        f.write(str(round(d, 2)) + '\n')

for pdb in pdbs:

    pr.calcTransformation(pdb.select('index 1 2 3 4 5 6 12 15 16 17 18'), centroid_pdb.select('index 1 2 3 4 5 6 12 15 16 17 18')).apply(pdb)
    
    pr.writePDB(outdir + pdb.getTitle(), pdb)

metal_coords, dists = ev.ev_atom_distribution(pdbs, centroid_pdb, atom_sel = 'name ND1 or name NE2', align_sel = 'index 1 2 3 4 5 6 12')

ev.plt_dist_his(dists, outdir + 'dists.png')

ev.plt_3d_points(metal_coords, outdir + 'metal3d.png')

for pdb in pdbs:

    pr.calcTransformation(pdb.select('index 1 2 3 4 5 6 12'), centroid_pdb.select('index 1 2 3 4 5 6 12')).apply(pdb)
    
    pr.writePDB(outdir + '_' + pdb.getTitle(), pdb)
