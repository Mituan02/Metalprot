'''
Generating database with 3 binding amino acids. 
For the ShaB DL metal position prediction method.
'''

import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract as ldb
from metalprot.database import database_cluster
import itertools

### set up parameters

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/"

metal_sel = 'name NI MN ZN CO CU MG FE' 
align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE'

cores = ldb.load_cores(workdir + '_Seq_core_date_reps/')

# Align 2 separate aa core

for c in cores:   
    c.generate_AAdAAdAA_Metal( key = 'AAdAAdAA_Metal')
    c.write_vdM(workdir + 'AAdAAdAAMetal_reps/',  key = 'AAdAAdAA_Metal')

_pdbs = ldb.get_all_pbd_prody(workdir + 'AAdAAdAAMetal_reps/')
database_cluster.run_cluster(_pdbs, workdir, 'AAdAAdAAMetal_reps_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 13, align_sel = align_sel_backbone, min_cluster_size = 0, tag = '')



