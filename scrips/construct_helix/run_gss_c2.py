'''
Try to generate c2 helix bundles with two different directions to build antiparallel 4 helix bundles.
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

query_dir = '/mnt/e/DesignData/ligands/CU_NI/'

print(query_dir)

querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['cluster', 'M3-3_HIS'], score_cut = 0, clu_num_cut = 10)

### 

workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C2/'

outdir = workdir + 'output_all4/'

metal_sel = 'name NI CU'

if not os.path.exists(outdir):
    os.mkdir(outdir)

#TO DO: write the ind to a txt file.
pairs = []
count = 1
for q in querys:
    #Filter the VdMs which are not helix
    phi_180, psi_180, seq = gss.cal_phipsi(q.query)
    dssps = gss.cal_dssp(phi_180, psi_180, seq)
    if len(dssps) > 1 and len([x for x in dssps if x =='A']) <= 2 :
        print(q.query.getTitle() + ' is not helix.')
        continue
    for i in range(-10, 11):
        x_rotation = i*3
        gss._generate_cn_symmetry_helix(outdir, q.query, name = 'C2_W_' + str(count), metal_sel = metal_sel, n = 2, x_rotations = [0 + x_rotation], y_rotations = [[0, 0]] )
        for j in range(-10, 11):
            y_rotations = [[-3*j, -3*j]]
            gss._generate_cn_symmetry_helix(outdir, q.query, name = 'C2_M_' + str(count) + '_'+ str(j+10), metal_sel = metal_sel, n = 2, x_rotations = [180 + x_rotation], y_rotations = y_rotations, z_rotation = 90 )
        pairs.append((count, q.query.getTitle()))
        count += 1

with open(outdir + '_index.txt', 'w') as f:
    for p in pairs:
        f.write(str(p[0]) + '\t' + p[1] + '\n')

#BUILD antiparallel 4 helix by c2 pair
def load_c2_pdbs(workdir):
    Ws = []
    Ms = []

    for pdb in os.listdir(workdir):
        if not pdb.endswith(".pdb"):
            continue
        if 'W' in pdb:
            Ws.append(pr.parsePDB(workdir + pdb))
        elif 'M' in pdb:
            Ms.append(pr.parsePDB(workdir + pdb))
    return Ws, Ms


def generate_antiparallel_4helix_by_c2pair(outdir, Ws, Ms):
    merge_pdbs = []

    for x in itertools.product(Ws, Ms):
        pdb1 = x[0].copy()
        pdb2 = x[1].copy()
        if gss.quick_clash(pdb1, pdb2, 5.0):
            print('clash')
            continue
        if not gss.filter_z_angle(pdb1) or not gss.filter_z_angle(pdb2):
            print('z angle too large')
            continue

        merge = gss.MergeAtomGroupAuto(outdir, [pdb1, pdb2], pdb1.getTitle().split('_x')[0] + '_' + pdb2.getTitle().split('x')[0], write_pdb = True)
        merge_pdbs.append(merge)
        print('merged.')

    return merge_pdbs

def temp_mv_exist_par(workdir, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    for par in os.listdir(workdir):
        if not par.endswith(".par"):
            continue
        os.rename(workdir + par, outdir + par)
        pdb = par.split('.')[0] + '.pdb'
        os.rename(workdir + pdb, outdir + pdb)
        txt = par.split('.')[0] + '.txt'
        os.rename(workdir + txt, outdir + txt)


'''
workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C2/output_all4/'
outdir = workdir + 'merges/'

if not os.path.exists(outdir):
    os.mkdir(outdir)

Ws, Ms = load_c2_pdbs(workdir)
generate_antiparallel_4helix_by_c2pair(outdir, Ws, Ms)

pdbs = []
for pdb in os.listdir(outdir):
    if not pdb.endswith(".pdb"):
        continue
    pdbs.append(pr.parsePDB(outdir + pdb))
write_XYZs(outdir, pdbs)
'''