'''
It would be interesting to  check the vdm distribution of vdmers regard to ABPLE and phi psi.
Still need to further calc the backbone ABPLE. 
'''

import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.search import search, search_eval
from metalprot.basic import filter
import pickle
from metalprot.database import database_extract as ldb
from metalprot.basic import utils, vdmer

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/search_eval/run_selfcenter_eval_search.py
'''

### without CYS
query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/pickle_all/'

with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)


outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/reason_CYS/'
with open(outdir + 'vdm_phipsi_abple_all.tsv', 'w') as f:
    f.write('title\ttype\tphi\tpsi\tabple\n')
    for v in all_querys:
        f.write(v.query.getTitle() + '\t' + v.aa_type + '\t' + str(round(v.phi, 3)) + '\t'+ str(round(v.psi, 3)) + '\t'+ v.abple + '\n')


metal_sel = 'ion or name NI MN ZN CO CU MG FE' 
'''
Pay attention to this function. The dehydral angle depend on the atoms we selected for each vdM.
'''
def calc_AB_angle(vdm):
    metal = vdm.query.select(metal_sel)[0]
    contact_atom = vdmer.get_contact_atom(vdm.query)
    third_atom = None
    forth_atom = None
    if vdm.aa_type == 'H':    
        third_atom = vdm.query.select('name CE1')[0]
        if contact_atom.getName() == 'ND1':
            forth_atom = vdm.query.select('name NE2')[0]

        elif contact_atom.getName() == 'NE2':
            forth_atom = vdm.query.select('name ND1')[0]
        else:
            forth_atom = vdm.query.select('name NE2')[0]
    elif vdm.aa_type == 'E':
        third_atom = vdm.query.select('name CD')[0]
        forth_atom = vdm.query.select('name OE1')[0]

    elif vdm.aa_type == 'D':
        third_atom = vdm.query.select('name CG')[0]
        forth_atom = vdm.query.select('name OD1')[0]

    elif vdm.aa_type == 'C':
        third_atom = vdm.query.select('name CB')[0]
        forth_atom = vdm.query.select('name CA')[0]
    try:
        A = pr.calcAngle(metal, contact_atom, third_atom)
        B = pr.calcDihedral(metal, contact_atom, third_atom, forth_atom)
        return A, B
    except:
        print(vdm.query.getTitle())
        print(third_atom)
        print(forth_atom)
        return None, None


def calc_AB_angle2(vdm):
    metal = vdm.query.select(metal_sel)[0]
    contact_atom = vdmer.get_contact_atom(vdm.query)
    third_atom = None
    forth_atom = None
    if vdm.aa_type == 'H':    
        #third_atom = vdm.query.select('name CE1')[0]
        if contact_atom.getName() == 'ND1':
            third_atom = vdm.query.select('name CG')[0]
            #forth_atom = vdm.query.select('name NE2')[0]
            forth_atom = vdm.query.select('name CB')[0]
        elif contact_atom.getName() == 'NE2':
            third_atom = vdm.query.select('name CD2')[0]
            #forth_atom = vdm.query.select('name ND1')[0]
            forth_atom = vdm.query.select('name CG')[0]
        else:
            third_atom = vdm.query.select('name CG')[0]
            forth_atom = vdm.query.select('name NE2')[0]
    elif vdm.aa_type == 'E':
        third_atom = vdm.query.select('name CD')[0]
        #forth_atom = vdm.query.select('name OE1')[0]
        forth_atom = vdm.query.select('name CG')[0]
    elif vdm.aa_type == 'D':
        third_atom = vdm.query.select('name CG')[0]
        #forth_atom = vdm.query.select('name OD1')[0]
        forth_atom = vdm.query.select('name CB')[0]
    elif vdm.aa_type == 'C':
        third_atom = vdm.query.select('name CB')[0]
        forth_atom = vdm.query.select('name CA')[0]
    try:
        A = pr.calcAngle(metal, contact_atom, third_atom)
        B = pr.calcDihedral(metal, contact_atom, third_atom, forth_atom)
        return A, B
    except:
        print(vdm.query.getTitle())
        print(third_atom)
        print(forth_atom)
        return None, None


# Extra angle
outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/reason_CYS/'
with open(outdir + 'vdm_contact_angles_all_dehedralCarbon.tsv', 'w') as f:
    f.write('title\ttype\tA\tB_dihidral\n')
    for v in all_querys:
        A, B = calc_AB_angle(v)
        if not A:
            continue
        f.write(v.query.getTitle() + '\t' + v.aa_type + '\t' + str(round(A, 3)) + '\t'+ str(round(B, 3)) + '\n')


'''
### vdm with +-3 aa (AAextMetal). To calc phi psi of 5 aa, and determine their abple.

workdir = "/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20211013/20211015_AAext3/"
pdbs = []
for dir in os.listdir(workdir):
    _pdbs = ldb.get_all_pbd_prody(workdir + dir + '/')
    pdbs.extend(_pdbs)

with open(outdir + 'vdmext_phipsi_abple.tsv', 'w') as f:
    f.write('title\ttype\tphi\tpsi\tabple\tphis\tpsis\tabples\tguess\n')
    for pdb in pdbs:
        try:
            abples, phipsi = utils.seq_get_ABPLE(pdb)
            f.write(pdb.getTitle() + '\t' + pdb.select('name CA').getResnames()[3] + '\t')
            f.write(str(round(phipsi[3][0],2)) + '\t')
            f.write(str(round(phipsi[3][1],2)) + '\t')
            f.write(abples[3] + '\t')
            f.write('||'.join([str(round(phipsi[i][0],2)) for i in range(7)]) + '\t')
            f.write('||'.join([str(round(phipsi[i][1],2)) for i in range(7)]) + '\t')
            f.write('||'.join([abples[i] for i in range(7)]) + '\t')
            guess = abples[3]
            if len([True for x in abples if x == guess]) < 3:
                guess = 'Lg'
            f.write(guess + '\n')
        except:
            print(pdb.getTitle())
'''