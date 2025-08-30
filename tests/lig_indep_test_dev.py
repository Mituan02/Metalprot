'''
# Testing clashing. loading files.
import prody as pr
from sklearn.neighbors import NearestNeighbors

workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/output_selfcenter_o1_1dmm_16-20-28_H-H-D_a_820__20220210-101859/represents_combs/'

target = pr.parsePDB(workdir + 'W_15-19-27_H-H-D_1000-404-467_allgly.pdb')

outdir = workdir + 'W_15-19-27_H-H-D_1000-404-467_/vdms_output/'

vdm = pr.parsePDB(outdir + 'phenol_PHE_102_2763.pdb')

outdir_ligs = workdir + 'W_15-19-27_H-H-D_1000-404-467_/ligs_inorder/'

lig = pr.parsePDB(outdir_ligs + 'Lig-150_Geo_5_tts_fe_C8-C7_255_C7-C6_355.pdb')

'''