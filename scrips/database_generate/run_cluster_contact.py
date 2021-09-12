
import os
import sys
import prody as pr
from metalprot import ligand_database as ldb

'''
python /mnt/e/GitHub_Design/Metalprot/metalprot/scrips/run_cluster_contact.py
'''


workdir = "/mnt/e/DesignData/ligands/ZN_rcsb/"

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

cores = ldb.load_cores(workdir + '2_seq_cores_reps/')


'''
key = 'TriAtomContact'

for c in cores:
    c.generate_binary_contact(key = key)

    c.write_vdM(workdir + 'M8_TriAtomContact_reps/',  key = key)

#Align atom core

_pdbs = ldb.get_all_pbd_prody(workdir + 'M8_TriAtomContact_reps/')

ldb.run_cluster(_pdbs, workdir, 'M8_TriAtomContact_clusters/', rmsd = 0.1, metal_sel = metal_sel, len_sel = 4, align_sel = 'all', min_cluster_size = 2, tag = 'm7_')

'''


'''
key = 'AtomContact'

for c in cores:
    c.generate_atom_contact(key = key)
    for i in range(4, 8):
        _key = key + str(i)
        #print(_key)
        c.write_vdM(workdir + 'M8_AtomContact' + str(i) +  '_reps/',  key = _key)
'''

_pdbs = ldb.get_all_pbd_prody(workdir + 'M8_AtomContact4_reps/')

ldb.run_cluster(_pdbs, workdir, 'M8_AtomContact4_clusters/', rmsd = 0.2, metal_sel = metal_sel, len_sel = 4, align_sel = 'all', min_cluster_size = 2, tag = 'm7_')
