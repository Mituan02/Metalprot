'''
After the design, we want to check how the vdMs different from the designed one. 
Are they still the orginal vdM and in the same vdM cluster?
'''

import os
import sys
from metalprot.apps.core import Core
import prody as pr
from metalprot import search_struct, extract_vdm
from metalprot import ligand_database as ldb

### Set working directories.

workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C4_rosetta/Rosetta_Output/'

#target_path = workdir + '3into4_helix_assembly_renum.pdb'

query_dir = '/mnt/e/DesignData/ligands/CU_NI/'


### Get query pdbs 

querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['cluster', 'M3-3_HIS'], score_cut = 0, clu_num_cut = 1)
print(len(querys))


### Get target VdMs

outdir = workdir + 'eva_out_CUNI_M3-3/'

if not os.path.exists(outdir):
    os.mkdir(outdir)

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

_pdb = ldb.get_all_pbd_prody(workdir)

cores = ldb.extract_all_core_seq_from_path(workdir, metal_sel, extend = 4, extract_2ndshell = False)

ldb.writepdb(cores, outdir + 'target_core/')

_pdbs = ldb.get_all_pbd_prody(outdir + 'target_core/')

all_aas = ['HIS', 'ASP', 'GLU', 'CYS']
for aa in all_aas:
    for p in _pdbs:
        core = Core(p)
        core.generate_AA_ext_Metal(aa, key = 'AAExtMetal_' + aa)
        core.write_vdM(outdir + 'target_core_aa/', key = 'AAExtMetal_' + aa)

### Match target 

def val_vdm(pdb, querys, align_sel, rmsd_cut = 1):
    pairs = []
    for query in querys:
        if len(pdb.select(align_sel)) == len(query.query.select(align_sel)):
            qr = query.copy()
            pr.calcTransformation(qr.query.select(align_sel), pdb.select(align_sel)).apply(qr.query)
            rmsd = pr.calcRMSD(pdb.select(align_sel), qr.query.select(align_sel))

            if rmsd <= rmsd_cut:
                pairs.append((pdb.getTitle(), qr, rmsd))
    return pairs


_pdbs = ldb.get_all_pbd_prody(outdir + 'target_core_aa/')

all_pairs = []
count = 0
for _pdb in _pdbs:
    pairs = val_vdm(_pdb, querys, 'name N C CA O NI')
    all_pairs.extend(pairs)
    for pair in pairs:
        pr.writePDB(outdir +  pair[0] + '_' + str(count)  + '_' + pair[1].query.getTitle(), pair[1].query)
        count += 1


with open(outdir + '_summary.txt', 'w') as f:
    f.write('rmsd\ttarget_vdm\tquery_vdm\tscore\tclust_num\tclust_tol_num\n')
    for pair in all_pairs:
        f.write(str(pair[2]) + '\t' + pair[0] + '\t' + pair[1].to_tab_string() + '\n')  


### Supperimpose all vdm members for the best match. 
outdir_mem = outdir + 'mems/'
if not os.path.exists(outdir_mem):
    os.mkdir(outdir_mem)

mem_querys = extract_vdm.get_vdm_mem(querys[1])

count = 0
for _pdb in _pdbs:

    pairs = val_vdm(_pdb, mem_querys, 'name N C CA O NI')
    for pair in pairs:
        pr.writePDB(outdir_mem +  pair[0] + '_' + str(count)  + '_' + pair[1].query.getTitle(), pair[1].query)
        count += 1
    with open(outdir_mem + '_summary.txt', 'w') as f:
        f.write('rmsd\ttarget_vdm\tquery_vdm\tscore\tclust_num\tclust_tol_num\n')
        for pair in pairs:
            f.write(str(pair[2]) + '\t' + pair[0] + '\t' + pair[1].to_tab_string() + '\n')  