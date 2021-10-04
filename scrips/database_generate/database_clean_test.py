import os
import sys
import shutil
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
import itertools
import numpy as np
from numpy.core.fromnumeric import argmin


root_dir = '/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20210624/'

workdir = root_dir + '_Seq_core_date_3HIScontact/'

pdbs = database_extract.get_all_pbd_prody(workdir)

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 


def get_metal_contact_atoms(pdb):
    '''
    Calc paired query angle and distance. 
    '''
    #Get metal, binding atom for each binding atom
    cts = []

    metal = pdb.select(metal_sel)[0]
    cts.append(metal)
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

    return cts

'''
### Development
pdb = pdbs[0]

cts = get_metal_contact_atoms(pdb)

bs = []
for c in cts:
    bs.append(c.getBeta())
    
'''

filtered_pdbs = []

cut_off = 45

for pdb in pdbs:
    cts = get_metal_contact_atoms(pdb)

    bs = []
    for c in cts:
        bs.append(c.getBeta())

    if all([b <= cut_off for b in bs]):
        filtered_pdbs.append(pdb)


### copy the filter pdb into a new folder.
outdir = root_dir + '_Seq_core_date_3HIScontact_B45/'
os.makedirs(outdir, exist_ok=True)
filtered_pdbs_dict = set([pdb.getTitle() for pdb in filtered_pdbs])
for p in os.listdir(workdir):
    if '.pdb' not in p:
        continue
    if p.split('.')[0] in filtered_pdbs_dict:
        shutil.copyfile(workdir + p, outdir + p)


### calc pairwise angle dist info.

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

dist_angle = []

for pdb in filtered_pdbs:
    pair_infos = calc_pair_geometry(pdb)
    for info in pair_infos:
        dist_angle.append((info[0][0], info[0][1], info[1], pdb.getTitle(), info[2], info[3]))


with open(workdir + '_pairwise_angle_dist_info_filteredB45.tsv', 'w') as f:
    f.write('aa\tinds\tTitle\tDist\tAngle\n')
    for da in dist_angle:
        x = da[0] + '-' + da[1]
        if da[0] < da[1]:
            x = da[1] + '-' + da[0]
        f.write(x + '\t' + '\t'.join([str(x) for x in da[2:]]) + '\n')