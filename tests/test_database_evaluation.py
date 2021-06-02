import pytest
import os
import shutil
import prody as pr
import pickle

from metalprot import ligand_database as ldb
from metalprot import database_evaluate as ev

# How to generate test pickle file from pdbs.
'''
workdir = '/mnt/e/DesignData/ligands/ZN_rcsb/M5_2aa_sep_cores_sc_cluster21/0/'
outdir = workdir + 'test_out/'
if not os.path.exists(outdir):
    os.mkdir(outdir)

pdbs = ldb.get_all_pbd_prody(workdir)

with open(outdir + 'all.pkl', 'wb') as outfile:
    pickle.dump(file=outfile, obj=pdbs)

infile = open(outdir + 'all.pkl','rb')
pdbs = pickle.load(infile)
infile.close()
'''

def test_db_ev():
    workdir = os.path.dirname(os.path.realpath(__file__)) + '/test_data/'
    infile = open(workdir + 'M1_AA5MetalSc_HIS_cluster02.pkl','rb')
    pdbs = pickle.load(infile)
    infile.close()
    
    for pdb in pdbs:
        if 'centroid' in pdb.getTitle():
            centroid_pdb = pdb

    outdir = workdir + 'test_db_ev/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    metal_sel = 'ion or name NI MN ZN CO CU MG FE' 
    metal_coords, dists = ev.ev_atom_distribution(pdbs, centroid_pdb, metal_sel, align_sel = None)
    ev.plt_dist_his(dists, outdir + 'dists.png')
    ev.plt_3d_points(metal_coords, outdir + 'metal3d.png')

    metal_coords, dists = ev.ev_atom_distribution(pdbs, centroid_pdb, metal_sel, align_sel = 'bb')
    ev.plt_dist_his(dists, outdir + 'dists_origin.png')
    ev.plt_3d_points(metal_coords, outdir + 'metal3d_origin.png')

    shutil.rmtree(outdir)
