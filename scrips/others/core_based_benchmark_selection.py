'''
Based 'E:\DesignData\ligands\ZN_rcsb_datesplit\20210624\reason\core_aa_info.tsv'
In folder 'E:\DesignData\ligands\ZN_rcsb_datesplit\20210624\_Seq_core_date_reps\'
Extract the 3-contact binding (HIS, ASP, GLU) core pdbs from and copy it into a new folder.

Then calculate the pairwise infomation. 
'''

import os
import shutil

root_dir = '/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20210624/'

workdir = root_dir + '_Seq_core_date_reps/'

outdir = root_dir + '_Seq_core_date_3contact/'

os.makedirs(outdir, exist_ok=True)

core_info_file = root_dir + 'reason/core_aa_info.tsv'

pdb_dict = {}

with open(core_info_file, 'r') as f:
    for line in f.readlines():
        info = line.split('\t')
        if len(info[1].split('-')) != 3:
            continue
        cs = info[1].split('\n')[0]
        sati = True
        for c in cs.split('-'):
            if not c in ['HIS', 'ASP', 'GLU']:
                sati = False
        if sati:      
            pdb_dict[info[0]] = cs

for p in os.listdir(workdir):
    if p.split('.')[0] in pdb_dict.keys():
        shutil.copyfile(workdir + p, outdir + p)




import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
import itertools
import numpy as np
from numpy.core.fromnumeric import argmin

'''
Extract binary metal binding dists and angles. 
'''
metal_sel = 'ion or name NI MN ZN CO CU MG FE' 


def calc_pair_geometry(pdb):
    '''
    Calc paired query angle and distance. 
    '''
    #Get metal, binding atom for each binding atom
    cts = []

    metal = pdb.select(metal_sel)[0]
    _contact_aas = pdb.select('protein and not carbon and not hydrogen and within 2.83 of resindex ' + str(metal.getResindex()))
    #For each aa, only select one contact atom. 
    resindices = np.unique(_contact_aas.getResindices())
    for rid in resindices:
        #TO DO: Not the first one, but the closest one for the  such ASP contains two contact atoms.
        _ct = _contact_aas.select('resindex ' + str(rid))
        if len(_ct) > 1:
            dists = [None]*len(_ct)
            for i in range(len(_ct)):
                dist = pr.calcDistance(metal, _ct[i])
                dists[i] = dist
            ct = _ct[argmin(dists)]
        else:
            ct = _ct[0]
        cts.append(ct)

    ct_len = len(cts)

    pair_infos = []

    for i, j in itertools.combinations(range(ct_len), 2):  
        pair_info = [] 

        pair_info.append((cts[i].getResname(), cts[j].getResname()))
        pair_info.append(str(cts[i].getResindex()) + '-' + str(cts[j].getResindex()))
        dist = pr.calcDistance(cts[i], cts[j])
        pair_info.append(dist)
        angle = pr.calcAngle(cts[i], metal, cts[j])
        pair_info.append(angle)

        pair_infos.append(pair_info)
    return pair_infos

workdir = outdir

pdbs = database_extract.get_all_pbd_prody(workdir)

dist_angle = []

for pdb in pdbs:
    pair_infos = calc_pair_geometry(pdb)
    for info in pair_infos:
        dist_angle.append((info[0][0], info[0][1], info[1], pdb.getTitle(), info[2], info[3]))


with open(workdir + '_pairwise_angle_dist_info.tsv', 'w') as f:
    f.write('aa\tinds\tTitle\tDist\tAngle\n')
    for da in dist_angle:
        x = da[0] + '-' + da[1]
        if da[0] < da[1]:
            x = da[1] + '-' + da[0]
        f.write(x + '\t' + '\t'.join([str(x) for x in da[2:]]) + '\n')
