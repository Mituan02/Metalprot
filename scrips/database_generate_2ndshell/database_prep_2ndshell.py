import os
import sys
import prody as pr
import shutil
#sys.path.append(r'/mnt/e/GitHub_Design/MetalDesign')
from metalprot.database import database_extract

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/"
metal_sel = 'name ZN'

# To get the core reps with 2nd shell, extract the core with 2nd shell first. Then match the name from the 1st shell core reps. 

cores_2nd = database_extract.extract_all_core_seq_from_path('/mnt/e/DesignData/ligands/ZN_rcsb/all_rcsb/', metal_sel, extend = 4, extract_2ndshell = True, _2ndshell_extend = 4)
os.makedirs(workdir + '20220116_2ndshell/_Seq_core_2ndshell/', exist_ok=True)
database_extract.superimpose_core_and_writepdb(cores_2nd, cores_2nd[0], metal_sel, workdir + '20220116_2ndshell/_Seq_core_2ndshell/')

# Check datesplit.py to change add the date.

### Copy based on 1st shell duplication.
core_reps = set()
for pdb_path in os.listdir('/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/_Seq_core_date_reps/'):
    if not pdb_path.endswith(".pdb"):
        continue
    core_reps.add(pdb_path)

outdir = workdir + '20220114_bb2ndshell/_Seq_cores_with_2ndshell_copyby1stshell/'
os.makedirs(outdir, exist_ok=True)

for pdb_path in os.listdir(workdir + '20220114_bb2ndshell/_Seq_core_2ndshell_date/'):
    #print(pdb_path)
    if not pdb_path.endswith(".pdb"):
        continue
    if pdb_path in core_reps:
        shutil.copy(workdir + '20220114_bb2ndshell/_Seq_core_2ndshell_date/' + pdb_path, outdir + pdb_path)


### Remove duplicate by 2nd shell.

metal_sel = 'name NI MN ZN CO CU MG FE' 

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/"

database_extract.rm_core(workdir + '_Seq_core_2ndshell_date/', workdir + '_Seq_core_2ndshell_date_filter/', min_contact_aa_num=2)

pdbs = database_extract.get_all_pbd_prody(workdir + '_Seq_core_2ndshell_date_filter/')

clusters = database_extract.reduce_dup(pdbs, metal_sel)

database_extract.write_dup_summary(workdir, pdbs, clusters)

database_extract.extract_rep_and_writepdb(pdbs, clusters, metal_sel, workdir + '_Seq_core_2ndshell_date_reps/')

