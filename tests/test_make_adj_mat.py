import os
import numpy as np
from numpy.core.defchararray import center
import prody as pr
from metalprot.basic import cluster
from prody.proteins.pdbfile import parsePDB



workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/AAMetalPhiPsi_HIS_cluster05/0/'

pdbs = []
coords = []


ind = 0
for file in os.listdir(workdir):
    if '.pdb' not in file:
        continue
    if 'centroid' in file:
        centroid = parsePDB(workdir + file) 
        centroid_ind = ind
        
    ind += 1
    pdb = parsePDB(workdir + file)
    pdbs.append(pdb)
    coords.append(pdb.select('heavy').getCoords())

pdb_coords = np.array(coords, dtype = 'float32')
rmsd_mat = cluster._make_pairwise_rmsd_mat_no_superposition(pdb_coords)


rmsds = []
for pdb in pdbs:
    rmsd = pr.calcRMSD(centroid.select('heavy'), pdb.select('heavy'))
    rmsds.append(rmsd)




