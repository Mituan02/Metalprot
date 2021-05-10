import os
import sys
import prody as pr

#You can either add the python package path.
sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import generate_sse as gss


workdir = '/mnt/e/DesignData/ligands/C4/'

target_path = workdir + 'm3-3_cluster_1_mem_20_centroid_3fms_NI_1_HIS_2.pdb'

outdir = workdir + 'output/'

metal_sel = 'name NI CU'

gss.generate_c4s(outdir, target_path, metal_sel)