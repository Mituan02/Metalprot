'''
SAHA is flexible ligand, if we can identify two cgs, 
one with metal binding, the other with vdm binding, 
we can search a possible conformation from MD generated ligs.

The procedure can be generalized for design of flexible ligand binding protein.
'''
import os
import numpy as np
import prody as pr

workdir = '/mnt/e/DesignData/Metalloenzyme/SAHA_Vorinostat/run_design_cgs3/SAHA_Rosetta20220730_2/'
outdir = workdir + 'MaestroLigs_1/'
os.makedirs(outdir, exist_ok=True)
lig_path = '/mnt/e/DesignData/Metalloenzyme/SAHA_Vorinostat/Maestro/'


lig_cg1 = pr.parsePDB(workdir + 'UNK.pdb')
lig_cg1_coords = np.array([lig_cg1.select('name ' + x).getCoords()[0] for x in ['N1', 'O1', 'C7', 'C9']])

lig_cg2 = pr.parsePDB(workdir + 'SAHA_cgs2_x2_y3_z0.pdb.gz')
lig_cg2_coords = np.array([lig_cg2.select('name ' + x).getCoords()[0] for x in ['N2', 'O2', 'O3', 'C8']])

lig_sel_coords = np.concatenate((lig_cg1_coords, lig_cg2_coords), axis = 0)

lig_0 = pr.parsePDB(lig_path + '5311.pdb')
lig_0_coords = np.array([lig_0.select('name ' + x).getCoords()[0] for x in ['N1', 'O1', 'C7', 'C9']])

tr = pr.calcTransformation(lig_0_coords, lig_cg1_coords)

ligs = []
for file in os.listdir(lig_path):
    if not '.pdb' in file:
        continue
    _lig = tr.apply(pr.parsePDB(lig_path + file))
    ligs.append(_lig)

min_rmsd = 3.0
for _lig in ligs:
    _lig_coords = np.array([_lig.select('name ' + x).getCoords()[0] for x in ['N1', 'O1', 'C7', 'C9', 'N2', 'O2', 'O3', 'C8']])
    tr = pr.calcTransformation(_lig_coords, lig_sel_coords).apply(_lig)
    _rmsd = pr.calcRMSD(_lig_coords, lig_sel_coords) 
    #print(_rmsd)
    if _rmsd < min_rmsd:
        min_rmsd = min_rmsd
        print(_lig.getTitle())
        pr.writePDB(outdir + _lig.getTitle(), _lig)

