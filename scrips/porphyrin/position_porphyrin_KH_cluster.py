
import os
import prody as pr
from metalprot.database import database_cluster 


workdir = '/Users/lonelu/DesignData/labmate_runs/KH_heme_240501/'

outdir = 'output_20240502_KH210_clu/'

_pdbs = dict()
for pdb_path in os.listdir(workdir + 'output_20240502_KH210/'):
    if 'biomol' not in pdb_path or '.pdb' not in pdb_path:
        continue
    key = '_'.join(pdb_path.split('_')[0:3])
    try:
        pdb = pr.parsePDB(workdir + 'output_20240502_KH210/' + pdb_path)
        if key not in _pdbs.keys():
            _pdbs[key] = [pdb]
        else:
            _pdbs[key].append(pdb)
    except:
        print(pdb_path)



#Cluster.
for key in _pdbs.keys():
    key = 'v_63_TYR'
    _outdir = outdir + key + '/'
    os.makedirs(_outdir, exist_ok=True)

    rmsd = 0.5
    metal_sel = 'FE'
    len_sel = 8
    align_sel = 'chid Y or (bb and chid X and resnum 10)'
    database_cluster.run_cluster(_pdbs[key], workdir, _outdir, rmsd, metal_sel, len_sel, align_sel)





############
#Test case
workdir = '/Users/lonelu/DesignData/labmate_runs/KH_heme_240501/'

outdir = 'output_20240502_KH210/'

lig_name = 'v_3_ARG_Lig_6_k0-2_r_0.48_1.29__1_1_6M4B_biomol_1_A_A.pdb'

lig = pr.parsePDB(workdir + outdir + lig_name)

lig.select('chid Y or (bb and chid X and resnum 10)').getNames()


###############################################
#Information Extraction.


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



workdir = '/Users/lonelu/DesignData/labmate_runs/KH_heme_240501/'
target_path = workdir + 'KH212.pdb'
target = pr.parsePDB(target_path)

# outdir = 'output_20240502_KH210/'
# lig_name = 'v_3_ARG_Lig_6_k0-2_r_0.48_1.29__1_1_6M4B_biomol_1_A_A.pdb'
# lig = pr.parsePDB(workdir + outdir + lig_name)
# lig.select('chid Y or (bb and chid X and resnum 10)').getNames()



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
_pdbs = dict()
for pdb_path in os.listdir(workdir + 'output_20240504_KH212/'):
    if '.pdb' not in pdb_path:
        continue
    is_lig = False
    if 'biomol' not in pdb_path:
        is_lig = True
        key = '_'.join(pdb_path.replace('-','_').split('_')[3:5])
    else:
        key = '_'.join(pdb_path.split('_')[3:5])

    pos_aa = '_'.join(pdb_path.split('_')[0:3])
    pos = '_'.join(pdb_path.split('_')[0:2])
    if key not in _pdbs.keys():
        c = vdM_info(lig_id = key, score = 0, score_pos_diff = 0, count = 0, aa_diff = 0, pos_diff = 0)
        _pdbs[key] = [c, None, [], dict()]
    #try:
    pdb = pr.parsePDB(workdir + 'output_20240504_KH212/' + pdb_path)

    if is_lig:
        _pdbs[key][1] = pdb
        continue
    else:   
        #vdM target clash.
        if vdm_target_clash(pdb, target, 2.6):
            clashed.append(pdb)
            continue

        score = float(pdb_path.split('_')[8])

        #update score
        if score > 0:
            _pdbs[key][0].score = _pdbs[key][0].score + score
        
        #update score_pos_diff
        if pos not in _pdbs[key][3].keys():
            _pdbs[key][3][pos] = max(0, score)
        else:
            _pdbs[key][3][pos] = max(score, _pdbs[key][3][pos])

        _pdbs[key][0].score_pos_diff = sum([_pdbs[key][3][_x] for _x in _pdbs[key][3].keys()])

        #update count
        _pdbs[key][0].count = _pdbs[key][0].count + 1
 
        #update count_pos_diff, count_aa_diff
        see_pos = set()
        see_aa = set()
        for _xvdm_pdb in _pdbs[key][2]:
            see_pos.add('_'.join(_xvdm_pdb.getTitle().split('_')[0:2]))
            see_aa.add('_'.join(_xvdm_pdb.getTitle().split('_')[0:3]))
        if pos not in see_pos:
            _pdbs[key][0].pos_diff = _pdbs[key][0].pos_diff + 1
        if pos_aa not in see_aa:
            _pdbs[key][0].aa_diff = _pdbs[key][0].aa_diff + 1
        
        _pdbs[key][2].append((pdb))
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


############################
list(_pdbs.keys())[0:5]

_pdbs['Lig_4131']



