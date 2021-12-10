'''
prepare metal-binding vdms. 1. Download pdbs, 2. Extract metal core, 3. filter the core of bfactor, 4. Remove duplication. 
Note that the filter is based on bfactor, occupancy and number of contact aa. 
 
'''
import os
import sys
import prody as pr
import shutil
from metalprot.database import database_extract

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/database_generate/database_preparation.py
'''

#------------------------------------------------------------------
'''
# Download metal containing pdbs.
'''
workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/"

filename='all_rcsb.txt'
database_extract.organize_rcsb_file(workdir)
database_extract.download_pdb(workdir, filename, resolution = 2.5)
# cd /mnt/e/DesignData/ligands/FE_rcsb/
# gunzip *.gz


#------------------------------------------------------------------
'''
# Extract metal core.
# All rcsb pids are extracted by each NI core with sequence (9 aa)
'''
workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/"
metal_sel = 'name ZN'

cores = database_extract.extract_all_core_seq_from_path(workdir + 'all_rcsb/', metal_sel, extend = 4, extract_2ndshell = False)

database_extract.superimpose_core_and_writepdb(cores, cores[0], metal_sel, workdir + '20210610/_Seq_cores/')


#------------------------------------------------------------------

'''
Here we want to filter the core we already extracted.
Filter is based on B factor, occupancy and number of contacting aa. 
'''
workdir += '20211207/'
database_extract.rm_core(workdir + '_Seq_cores/', workdir + '_Seq_core_filter/', min_contact_aa_num=2)

#------------------------------------------------------------------

# cluster based on core sequence. write summary and extract represent (the first on in each cluster)

pdbs = database_extract.get_all_pbd_prody(workdir + '_Seq_core_filter/')

clusters = database_extract.reduce_dup(pdbs, metal_sel)

database_extract.write_dup_summary(workdir, pdbs, clusters)

database_extract.extract_rep_and_writepdb(pdbs, clusters, metal_sel, workdir + '_Seq_core_reps/')

