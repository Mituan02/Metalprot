'''
For the two vdMs that are connecting with each other, 
extract the connect aa type, 
calc the number of aa between them. 
It turned out there are some AAcAA bivdms are next to each other.
'''

import os
import sys

from metalprot.basic.quco import Query, Comb
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract



workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/"

cores = database_extract.load_cores(workdir + 'AAcAAMetal_reps/')

#pdbs = database_extract.get_all_pbd_prody(workdir)

infos = []

for c in cores:
    if len(c.contact_aa_resinds) < 2:
        print(c.full_pdb.getTitle())
        continue
    num = abs(c.contact_aa_resinds[1] - c.contact_aa_resinds[0]) - 1
    aa0 = c.full_pdb.select('name CA and resindex ' + str(c.contact_aa_resinds[0])).getResnames()[0]
    aa1 = c.full_pdb.select('name CA and resindex ' + str(c.contact_aa_resinds[1])).getResnames()[0]
    aaCaa = aa0 + '-' + aa1
    if aa0 < aa1:
        aaCaa = aa1 + '-' + aa0
    infos.append((c.full_pdb.getTitle(), num, aa0, aa1, aaCaa))

with open(workdir + 'AAcAAMetal_reps/_AAcAA_info.tsv', 'w') as f:
    f.write('Title\tNum\tAA0\tAA1\tAAcAA\n')
    for info in infos:
        f.write('\t'.join([str(x) for x in info]) + '\n')


