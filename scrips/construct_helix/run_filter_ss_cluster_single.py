'''
In this script, I intended to extract best helix histine vdms.
In ...\ZN_rcsb_datesplit\20211013\20211015_AAext3, we filter ss first by copy the vdms with abple 'XXAAAXX' to a new folder.
Then run the extraction and clustering. 
Then align the vdms on HIS-sc.
Then add vdms in the middle of 15mer helixs.
'''
import os
from unittest import TestResult
from numpy import append
import prody as pr
import shutil
import pickle
from metalprot.database import database_extract
from metalprot.basic import utils, prody_ext
from metalprot.database import database_cluster as ldb_clu


workdir = "/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20211013/20211015_AAext3/"

### Filter AAA
outdir = workdir + 'AAextMetal_HIS_reps_AAA/'
os.makedirs(outdir, exist_ok=True)

_pdbs = database_extract.get_all_pbd_prody(workdir + 'AAextMetal_HIS_reps/')


for pdb in _pdbs:
    abples, phipsi = utils.seq_get_ABPLE(pdb)
    if abples[2] == 'A' and abples[3] == 'A' and abples[4] == 'A':
        shutil.copy2(workdir + 'AAextMetal_HIS_reps/' + pdb.getTitle() + '.pdb', outdir + pdb.getTitle() + '.pdb')

### Extract Core
cores = database_extract.load_cores(outdir)

core_outdir = workdir + 'AAextMetal_HIS_reps_AAA_Core/'
os.makedirs(outdir, exist_ok=True)

for c in cores:   
    c.generate_AA_phipsi_Metal('HIS', key = 'AAA_H_reps')
    c.write_vdM(core_outdir, key = 'AAA_H_reps')


### Cluster

_vdm_pdbs = database_extract.get_all_pbd_prody(core_outdir)
#ldb_clu.run_cluster(_vdm_pdbs, workdir, 'AAA_H_cluster/', rmsd = 0.5, metal_sel = 'name NI MN ZN CO CU FE' , len_sel = 15, align_sel = 'heavy', min_cluster_size = 0, tag = 'AAA_', is_self_center=True)



ldb_clu.run_cluster(_vdm_pdbs, workdir, 'AAA_H_cluster_combine/', rmsd = 0.5, metal_sel = 'name NI MN ZN CO CU FE' , len_sel = 15, align_sel = 'heavy', min_cluster_size = 0, tag = 'AAA_', is_self_center=False)

'''
# Copy the centroid out and superimpose on one bb.
clu_comb_cent_outdir = workdir + 'AAA_H_cluster_combine_cent/'
os.makedirs(clu_comb_cent_outdir, exist_ok=True)


for xx in os.listdir(workdir + 'AAA_H_cluster_combine/'):
    if not os.path.isdir(workdir + 'AAA_H_cluster_combine/' + xx):
        continue
    for x in os.listdir(workdir + 'AAA_H_cluster_combine/' + xx):
        if not '.pdb' in x or not 'centroid' in x:
            continue
        if 'cluster_0' in x:
            pdb0 = pr.parsePDB(workdir + 'AAA_H_cluster_combine/' + xx + '/' + x)
        pdb = pr.parsePDB(workdir + 'AAA_H_cluster_combine/' + xx + '/' + x)
        pr.calcTransformation(pdb.select('resindex 1 and name N CA C'), pdb0.select('resindex 1 and name N CA C')).apply(pdb)
        pr.writePDB(clu_comb_cent_outdir + pdb.getTitle() + '.pdb', pdb)
        #shutil.copy(workdir + 'AAA_H_cluster_combine/' + xx + '/' + x, clu_comb_cent_outdir + x)
'''
    
### Extract vdm.
# Before running the following code, run database_selfcenter_pickle.py first to extract the vdM infomation.
# Then load the vdM info.
query_dir = workdir + 'pickle_AAA_H/'
with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)

vdm_dir = workdir + 'AAA_H_vdms/' 
os.makedirs(vdm_dir, exist_ok=True)

for v in all_querys:
    title = '_'.join(str(v.query.getTitle().split('_')[i]) for i in range(3)) + '_sc_' + str(round(v.score, 2))
    pr.writePDB(vdm_dir + title, v.query)

### Align all vdms on HIS sc ring (without metal).

_vdms = database_extract.get_all_pbd_prody(vdm_dir)

vdm_align_dir = workdir + 'AAA_H_vdms_alignSC/' 
os.makedirs(vdm_align_dir, exist_ok=True)
v_0 = _vdm_pdbs[0]
for v in _vdms:
    pr.calcTransformation(v.select('name CG ND1 CE1 CD2 NE2'), v_0.select('name CG ND1 CE1 CD2 NE2')).apply(v)
    pr.writePDB(vdm_align_dir + v.getTitle() + '.pdb', v)

# Attach 15mer. 

_15mer_path = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/15mer_ALA.pdb'
_15mer = pr.parsePDB(_15mer_path)

vdmHelix_dir = workdir + 'AAA_H_vdms_helix/' 
os.makedirs(vdmHelix_dir, exist_ok=True)

_vdms_alignSC = database_extract.get_all_pbd_prody(vdm_align_dir)
for v in _vdms_alignSC:
    _15mer_cp = _15mer.copy()
    pr.calcTransformation(_15mer_cp.select('resindex 7 and name N C CA'), v.select('resindex 1 and name N C CA')).apply(_15mer_cp)
    title = v.getTitle()
    resind_vdm_dict = {7: v} 
    ag = prody_ext.mutate_vdm_target_into_ag(_15mer_cp, resind_vdm_dict, title = title, write_metal=True)
    pr.writePDB(vdmHelix_dir + ag.getTitle() + '.pdb', ag)