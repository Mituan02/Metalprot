'''
Try to generate c2 helix bundles with two different directions to build antiparallel 4 helix bundles.
'''

import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import generate_sse as gss
from metalprot import extract_vdm


# Generate queryss

query_dir = '/mnt/e/DesignData/ligands/CU_NI/'

print(query_dir)

querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['cluster', 'M3-3_HIS'], score_cut = 0, clu_num_cut = 10)



workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C2/'

outdir = workdir + 'output/'

metal_sel = 'name NI CU'

if not os.path.exists(outdir):
    os.mkdir(outdir)

for q in querys:
    phi_180, psi_180, seq = gss.cal_phipsi(q.query)
    dssps = gss.cal_dssp(phi_180, psi_180, seq)
    if len(dssps) > 1 and len([x for x in dssps if x =='A']) <= 2 :
        print(q.query.getTitle() + ' is not helix.')
        continue

    x_rotations = [0]
    x2_rotations = [[0, 0]]

    gss._generate_cn_symmetry_helix(outdir, q.query, name = 'C2_W_' + q.query.getTitle(), metal_sel = metal_sel, n = 2, x_rotations = [0], x2_rotations = x2_rotations )
    gss._generate_cn_symmetry_helix(outdir, q.query, name = 'C2_M_' + q.query.getTitle(), metal_sel = metal_sel, n = 2, x_rotations = [180], x2_rotations = x2_rotations, z_rotation = 90 )

    
