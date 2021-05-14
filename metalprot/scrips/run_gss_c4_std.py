'''
Generate a std CC helix witch contains 7 aa. 
extract the HIS-Metal M1 VdMs and try to match the bb, supperimpose them to the 3r aa of the std helix by bb, keep those with RMSD <= 0.15
Mutate the std CC at 3r to the VdM to generate queries.
Use the queries to generate C4.
'''

import os
import sys
import prody as pr
import itertools
from scipy.spatial.distance import cdist
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import generate_sse as gss
from metalprot import extract_vdm, search_struct

### Extract queryss

def generate_c4s(outdir, target_path, metal_sel, pre_flag = 'target_c4_', delta = 3):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    target = pr.parsePDB(target_path)

    x_rotations = []
    y_rotations = []
    for i in range(-20, 21):
        x_rotations.append(i*delta)

    for j in range(-10, 11):
        y_rotations.append([j*delta, 0, j*delta,  0])

    gss._generate_cn_symmetry_helix(outdir, target, name = pre_flag, metal_sel = metal_sel, n = 4, x_rotations = x_rotations, y_rotations = y_rotations )

### Extract querys

query_dir = '/mnt/e/DesignData/ligands/CU_NI/'

print(query_dir)

querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['cluster', 'M1-1_HIS_sc'], score_cut = 0, clu_num_cut = 20)

for q in querys:
    mobile = q.query.copy()
    pr.calcTransformation(mobile.select('bb'), target.select('resindex 3 and bb')).apply(mobile)
    if pr.calcRMSD(mobile.select('bb'), target.select('resindex 3 and bb')) <= 0.15:
        pr.writePDB(workdir + mobile.getTitle(), mobile)

### supperimpose 1 aa VdM to std helix with 7 aas.

workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C4_std/'

outdir = workdir + 'output_a/'

metal_sel = 'name NI CU'

if not os.path.exists(outdir):
    os.mkdir(outdir)

target = pr.parsePDB(workdir + 'std_a.pdb')

pdbs = []
for pdb in os.listdir(workdir):
    if not pdb.startswith('m1-1') or not pdb.endswith(".pdb"):
        continue
    pdbs.append(pr.parsePDB(workdir + pdb))

count = 0
for pdb in pdbs:
    _target = target.copy()
    t1 = _target.select('resindex 0 1 2').toAtomGroup()

    t2 = pdb.select('protein').toAtomGroup()
    t2.setResnums([3 for i in range(len(t2))])
    t2.setChids([t1.getChids()[0] for i in range(len(t2))])
    t2.setSegnames([t1.getSegnames()[0] for i in range(len(t2))])

    t3 = _target.select('resindex 4 5 6').toAtomGroup()  

    t4 = pdb.select(metal_sel).toAtomGroup()  
    t4.setResnums([1000])
    t4.setChids([t1.getChids()[0] ])
    t4.setSegnames([t1.getSegnames()[0]])

    t = t1 + t2 + t3 + t4
    pr.writePDB(workdir + 'candidate_' + str(count)  + '.pdb', t)
    count += 1

### Generate c4s

candidates = []
for pdb in os.listdir(workdir):
    if not pdb.startswith('candidate_') or not pdb.endswith(".pdb"):
        continue
    candidates.append(workdir + pdb)

count = 0
for c in candidates:
    generate_c4s(outdir, c, metal_sel, pre_flag = 'cdd' + str(count))
    count += 1

### write XYZ for CCCP-fitting

pdbs = []
for pdb in os.listdir(outdir):
    if not pdb.endswith(".pdb"):
        continue
    pdbs.append(pr.parsePDB(outdir + pdb))

write_XYZs(outdir, pdbs)