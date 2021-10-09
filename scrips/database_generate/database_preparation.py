import os
import sys
import prody as pr
import shutil
#sys.path.append(r'/mnt/e/GitHub_Design/MetalDesign')
from metalprot.database import database_extract as ldb


#------------------------------------------------------------------
'''
# Download metal containing pdbs.

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/"

filename='all_rcsb.txt'
ldb.organize_rcsb_file(workdir)
ldb.download_pdb(workdir, filename, resolution = 2.5)
# cd /mnt/e/DesignData/ligands/FE_rcsb/
# gunzip *.gz
'''


#------------------------------------------------------------------

# Extract metal core.
# All rcsb pids are extracted by each NI core with sequence (9 aa)

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/"
metal_sel = 'name ZN'

cores = ldb.extract_all_core_seq_from_path(workdir + 'all_rcsb/', metal_sel, extend = 4, extract_2ndshell = False)

ldb.superimpose_core_and_writepdb(cores, cores[0], metal_sel, workdir + '20210610/_Seq_cores/')


#------------------------------------------------------------------


# cluster based on core sequence. write summary and extract represent (the first on in each cluster)

pdbs = ldb.get_all_pbd_prody(workdir + '_Seq_core_filter/')

clusters = ldb.reduce_dup(pdbs, metal_sel)

ldb.write_dup_summary(workdir, pdbs, clusters)

ldb.extract_rep_and_writepdb(pdbs, clusters, metal_sel, workdir + '_Seq_core_date_reps/')

