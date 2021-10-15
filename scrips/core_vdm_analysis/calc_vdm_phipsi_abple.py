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