import os
import sys
import prody as pr
sys.path.append(r'/mnt/e/GitHub_Design/MetalDesign')
from metalprot import ligand_database as ldb


#------------------------------------------------------------------
'''
# Download metal containing pdbs.

workdir = "/mnt/e/DesignData/ligands/CA_rcsb/"

filename='all_rcsb.txt'
ldb.organize_rcsb_file(workdir)
ldb.download_pdb(workdir, filename, resolution = 2.5)
# cd /mnt/e/DesignData/ligands/FE_rcsb/
# gunzip *.gz
'''


#------------------------------------------------------------------

# Extract metal core.
# All rcsb pids are extracted by each NI core with sequence (9 aa)

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb/"
metal_sel = 'name ZN'

#pdbs = get_all_pbd_prody(workdir + 'all_rcsb/')
#cores = extract_all_core_seq(pdbs, metal_sel, extend = 4)

cores = ldb.extract_all_core_seq_from_path(workdir + 'all_rcsb/', metal_sel, extend = 4, extract_2ndshell = True)

ldb.superimpose_core_and_writepdb(cores, cores[0], metal_sel, workdir + 'Seq_cores_with_2ndshell/')


#------------------------------------------------------------------


# cluster based on core sequence. write summary and extract represent (the first on in each cluster)

pdbs = ldb.get_all_pbd_prody(workdir + 'Seq_cores_with_2ndshell/')

clusters = ldb.reduce_dup(pdbs, metal_sel)

ldb.write_dup_summary(workdir, pdbs, clusters)

ldb.extract_rep_and_writepdb(pdbs, clusters, metal_sel, workdir + 'Seq_cores_with_2ndshell_reps/')

#------------------------------------------------------------------