
import os
import prody as pr
from dataclasses import dataclass
import numpy as np 
from sklearn.neighbors import NearestNeighbors

@dataclass
class vdM_info:
    lig_id: str
    score: float
    score_pos_diff:float
    count: int
    aa_diff: int
    pos_diff: int



workdir = '/Users/lonelu/DesignData/APEX/Design_210/'
target_path = workdir + 'KH210_af3_add5rf.pdb'
target = pr.parsePDB(target_path)
ligand = pr.parsePDB(workdir + 'hem_trimless_af3.pdb')
# outdir = 'output_20240502_KH210/'
# lig_name = 'v_3_ARG_Lig_6_k0-2_r_0.48_1.29__1_1_6M4B_biomol_1_A_A.pdb'
# lig = pr.parsePDB(workdir + outdir + lig_name)
# lig.select('chid Y or (bb and chid X and resnum 10)').getNames()


def vdm_ligand_clash(vdm, ligand, clash_radius = 1.2):
    '''
    The vdm Y sidechain may clash with the ligand.
    return True if clash.
    '''
    vdm_coords = vdm.select('heavy and sc and resnum 10').getCoords()

    ligand_coords = ligand.select('heavy').getCoords()

    if len(vdm_coords) == 0:
        print(vdm[['CG', 'rota', 'probe_name', 'chain', 'name']])
        return False
    nbr = NearestNeighbors(radius=clash_radius).fit(np.array(vdm_coords))
    adj_matrix = nbr.radius_neighbors_graph(ligand_coords).astype(bool)

    if adj_matrix.sum() > 0:
        return True
    return False


def vdm_target_clash(vdm, target, clash_radius = 2.6):
    '''
    The vdm X sidechain may clash with the target.
    return True if clash.
    '''
    if vdm.select('chid X and resnum 10').getResnames()[0] == 'GLY':
        return False
    vdm_sel = vdm.select('resnum 10 and sc and heavy')
    # if vdm_sel == None:
    #     print('----------------------------')
    #     print(vdm.getTitle())
    #     return False
    vdm_coords = vdm_sel.getCoords()
    resnum = vdm.getTitle().split('_')[1]
    target_coords = target.select('name N C CA O and not resnum ' + str(resnum)).getCoords()

    nbr = NearestNeighbors(radius=clash_radius).fit(np.array(vdm_coords))
    adj_matrix = nbr.radius_neighbors_graph(target_coords).astype(bool)

    if adj_matrix.sum() > 0:
        return True
    return False

clashed = []

for pdb_path in os.listdir(workdir + 'output_240823/'):
    if '.pdb' not in pdb_path:
        continue
    pdb = pr.parsePDB(workdir + 'output_240823/' + pdb_path)

    if True:
        #vdM target clash.
        if vdm_ligand_clash(pdb, ligand, 1.2):
            clashed.append(pdb)
            continue
len(clashed)
        
clashed = []
_pdbs = dict()
for pdb_path in os.listdir(workdir + 'output_240823/'):
    if '.pdb' not in pdb_path:
        continue
    #'vdm_3_ARG_Lig-1_key_0-0_rmsd_0.73_0.4__6_1_3AP1_biomol_1_A_A.pdb'
    #'v_66_ARG_K0-0_C_0.26__2_1_6R2I_biomol_1_A_A.pdb'
    key = '_'.join(pdb_path.split('_')[0:3])

    pos_aa = '_'.join(pdb_path.split('_')[0:3])
    pos = '_'.join(pdb_path.split('_')[0:2])
    if key not in _pdbs.keys():
        c = vdM_info(lig_id = key, score = 0, score_pos_diff = 0, count = 0, aa_diff = 0, pos_diff = 0)
        _pdbs[key] = [c, None, [], dict()]
    #try:
    pdb = pr.parsePDB(workdir + 'output_240823/' + pdb_path)

    if True:
        #vdM target clash.
        if vdm_ligand_clash(pdb, target, 1.2):
            clashed.append(pdb)
            continue

        score = float(pdb_path.split('_')[5])

        #update score_pos_diff

        #update count
        _pdbs[key][0].count = _pdbs[key][0].count + 1
 
    # except:
    #     print('error: {}'.format(pdb_path))

outdir = 'output_20240504_KH212_ext/'
os.makedirs(workdir + outdir, exist_ok= True)

with open(workdir + outdir + 'summary.tsv', 'w') as f:
    f.write('lig_id\tscore\tscore_pos_diff\tcount\taa_diff\tpos_diff\n')
    for key in _pdbs.keys():
        f.write(_pdbs[key][0].lig_id + '\t' 
                          + str(_pdbs[key][0].score) + '\t' 
                          + str(_pdbs[key][0].score_pos_diff) + '\t'
                          + str(_pdbs[key][0].count) + '\t'
                          + str(_pdbs[key][0].aa_diff) + '\t'
                          + str(_pdbs[key][0].pos_diff) + '\n')


import pandas as pd

table = pd.read_table(workdir + outdir + 'summary.tsv')

table.columns

def getLigTitle(x):
    xs = x.split('_')
    xss= []
    xss.extend(xs[3:5])
    xss.extend(xs[5:])
    y = '_'.join(xss)
    return y

def getVdmTitle(x):
    xs = x.split('_')
    xss= []
    xss.extend(xs[3:5])
    xss.extend(xs[0:3])
    xss.extend(xs[5:])
    y = '_'.join(xss)
    return y

for key in table['score_pos_diff'].nlargest(n=30).keys():
    lig_key = table.iloc[key]['lig_id']
    titles = _pdbs[lig_key][1].getTitle().split('mem')
    xtitle = titles[0].replace('-', '_') + 'mem' + titles[1]
    pr.writePDB(workdir + outdir + getLigTitle(xtitle), _pdbs[lig_key][1])

    for _x in _pdbs[lig_key][2]:
        pr.writePDB(workdir + outdir + getVdmTitle(_x.getTitle()), _x)

outdir = 'output_20240504_KH212_ext_pos/'
os.makedirs(workdir + outdir, exist_ok= True)

for key in table['pos_diff'].nlargest(n=30).keys():
    lig_key = table.iloc[key]['lig_id']
    titles = _pdbs[lig_key][1].getTitle().split('mem')
    xtitle = titles[0].replace('-', '_') + 'mem' + titles[1]
    pr.writePDB(workdir + outdir + getLigTitle(xtitle), _pdbs[lig_key][1])

    for _x in _pdbs[lig_key][2]:
        pr.writePDB(workdir + outdir + getVdmTitle(_x.getTitle()), _x)
