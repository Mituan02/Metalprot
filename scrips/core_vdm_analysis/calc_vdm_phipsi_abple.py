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
from metalprot.basic import utils

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/search_eval/run_selfcenter_eval_search.py
'''

### without CYS
query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/pickle_noCYS_alignBB/'

with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)


outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/reason/'
with open(outdir + 'vdm_phipsi_abple_noCYS.tsv', 'w') as f:
    f.write('title\ttype\tphi\tpsi\tabple\n')
    for v in all_querys:
        f.write(v.query.getTitle() + '\t' + v.aa_type + '\t' + str(round(v.phi, 3)) + '\t'+ str(round(v.psi, 3)) + '\t'+ v.abple + '\n')

### with CYS
query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/pickle_alignBB/'

with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)


outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/reason/'
with open(outdir + 'vdm_phipsi_abple_all.tsv', 'w') as f:
    f.write('title\ttype\tphi\tpsi\tabple\n')
    for v in all_querys:
        f.write(v.query.getTitle() + '\t' +  v.aa_type + '\t' + str(round(v.phi, 3)) + '\t'+ str(round(v.psi, 3)) + '\t'+ v.abple + '\n')


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
