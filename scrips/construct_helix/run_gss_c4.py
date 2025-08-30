import os
import sys
import prody as pr

#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import generate_sse as gss
from metalprot import extract_vdm


def generate_c4s(outdir, target_path, metal_sel, pre_flag = 'target_c4_', x = 20, y = 10, delta = 3):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    target = pr.parsePDB(target_path)

    x_rotations = []
    y_rotations = []
    for i in range(-x, x+1):
        x_rotations.append(i*delta)

    for j in range(-y, y+1):
        y_rotations.append([j*delta, 0, j*delta,  0])

    gss._generate_cn_symmetry_helix(outdir, target, name = pre_flag, metal_sel = metal_sel, n = 4, x_rotations = x_rotations, y_rotations = y_rotations )

'''
### Generate for a special one
workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C4/'
outdir = workdir + 'output_delta1/'
target_path = workdir + 'm3-3_cluster_1_mem_20_centroid_3fms_NI_1_HIS_2.pdb'

name = 'C4_'+ str(count) + '_' + target_path.split('/')[-1].split('centroid_')[-1].split('.')[0]
generate_c4s(outdir, target_path, metal_sel, name, x = 60, y = 30, delta = 1)
'''


### Extract querys

query_dir = '/mnt/e/DesignData/ligands/CU_NI/'
workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C4_all_CYS/'

print(query_dir)

querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['cluster', 'M3-3_CSY'], score_cut = 0, clu_num_cut = 3)

for q in querys:
    phi_180, psi_180, seq = gss.cal_phipsi(q.query)
    dssps = gss.cal_dssp(phi_180, psi_180, seq)
    if len(dssps) > 1 and len([x for x in dssps if x =='A']) <= 2 :
        print(q.query.getTitle() + ' is not helix.')
        continue
    pr.writePDB(workdir + q.query.getTitle(), q.query)

### 

workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C4_all_CYS/'

outdir = workdir + 'output_x40_y20/'

metal_sel = 'name NI CU'

if not os.path.exists(outdir):
    os.mkdir(outdir)

### Generate c4s

pdb_paths = []
for pdb in os.listdir(workdir):
    if not pdb.endswith(".pdb"):
        continue
    pdb_paths.append(workdir + pdb)

count = 0
for target_path in pdb_paths:
    name = 'C4_'+ str(count) + '_' + target_path.split('/')[-1].split('centroid_')[-1].split('.')[0]
    generate_c4s(outdir, target_path, metal_sel, name)
    count += 1

### write XYZ for CCCP-fitting

pdbs = []
for pdb in os.listdir(outdir):
    if not pdb.endswith(".pdb"):
        continue
    pdbs.append(pr.parsePDB(outdir + pdb))

gss.write_XYZs(outdir, pdbs)