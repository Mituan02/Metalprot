import os
import sys
import prody as pr
import shutil
#sys.path.append(r'/mnt/e/GitHub_Design/MetalDesign')
from metalprot import ligand_database as ldb

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/"
metal_sel = 'name ZN'

cores = ldb.extract_all_core_seq_from_path(workdir + 'all_rcsb/', metal_sel, extend = 4, extract_2ndshell = False)

ldb.superimpose_core_and_writepdb(cores, cores[0], metal_sel, workdir + '20210610/_Seq_cores/')


# To get the core reps with 2nd shell, extract the core with 2nd shell first. Then match the name from the core reps w/o 2nd shell. 

cores_2nd = ldb.extract_all_core_seq_from_path(workdir + 'all_rcsb/', metal_sel, extend = 4, extract_2ndshell = True)

ldb.superimpose_core_and_writepdb(cores_2nd, cores_2nd[0], metal_sel, workdir + '20210608/_Seq_cores/')

core_reps = set()
for pdb_path in os.listdir(workdir + '20210608/_Seq_cores_reps/'):
    if not pdb_path.endswith(".pdb"):
        continue
    core_reps.add(pdb_path)

outdir = workdir + '20210608/_Seq_cores_with_2ndshell_reps/'
if not os.path.exists(outdir):
    os.mkdir(outdir)

for pdb_path in os.listdir(workdir + '20210608/_Seq_cores_with_2ndshell/'):
    if not pdb_path.endswith(".pdb"):
        continue
    if pdb_path in core_reps:
        shutil.copy(workdir + '20210608/_Seq_cores_with_2ndshell/' + pdb_path, outdir + pdb_path)