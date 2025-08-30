import os
import prody as pr
import pickle


query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/20210923_bb/M1-1_AAMetalSc_HIS_cluster05/'
his_querys = []
for path in os.listdir(query_dir):
    if '.' in path:
        continue
    #print(lib_dir + path)
    for file in os.listdir(query_dir + path):
        if '.pdb' not in file:
            continue
        #print(lib_dir + path + '/' + file)
        pdb = pr.parsePDB(query_dir + path + '/' + file)
        his_querys.append((pdb, path))


lib_dir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/M1-1_AAMetalSc_HIS_cluster05/'

core_pdb_reps = []
for path in os.listdir(lib_dir):
    if '.' in path:
        continue
    #print(lib_dir + path)
    for file in os.listdir(lib_dir + path):
        if '.pdb' not in file:
            continue
        #print(lib_dir + path + '/' + file)
        pdb = pr.parsePDB(lib_dir + path + '/' + file)
        core_pdb_reps.append((pdb, path))

def porphyrin_superimpose_sel(pdb, metal_sel):
    metal_ind = pdb.select(metal_sel)[0].getIndex()
    contact_atom_resind = pdb.select('protein and within 2.83 of index ' + str(metal_ind))[0].getResindex()
    contact_aa_bb_inds = pdb.select('heavy and resindex ' + str(contact_atom_resind)).getIndices()
    sel = pdb.select('index ' + str(metal_ind) + ' ' + ' '.join([str(x) for x in contact_aa_bb_inds]))
    return sel


metal_sel = 'name NI MN ZN CO CU MG FE' 

#his_querys = [q for q in all_querys if 'HIS' in q.query.getTitle()]

rmsds = []

for c_c in core_pdb_reps:
    c= c_c[0]
   
    #core_sel = c.select('protein and ' + metal_sel)
    #core_sel = porphyrin_superimpose_sel(c, metal_sel)
    core_sel = c.select('protein and heavy or ' + metal_sel)

    min_rmsd = 5
    min_q = None

    for q in his_querys:

        #q_sel = q.query.select('resindex ' + str(q.contact_resind) + ' and heavy or ' + metal_sel)
        q_sel = q[0].select('heavy')

        if len(q_sel) != len(core_sel):
            print('Failed superimpose: ' + c.getTitle())
            continue

        tr = pr.calcTransformation(q_sel, core_sel)
        tr.apply(q_sel)

        rmsd = pr.calcRMSD(q_sel, core_sel)
        
        if rmsd < min_rmsd:
            min_rmsd = rmsd
            min_q = q[0].copy()

    rmsds.append((c_c, min_rmsd, min_q))

workdir = '/mnt/e/DesignData/ligands/porphyrin/database_compare/'

os.makedirs(workdir, exist_ok=True)

with open(workdir + 'database_compare.tsv', 'w') as f:
    for r in rmsds:
        f.write(r[0][0].getTitle() + '\t' + str(r[1]) + '\t' + r[2].getTitle() + '\n')

count = 0
for r in rmsds:
    #if 'centroid' in r[0][0].getTitle():
    pr.writePDB(workdir + str(count) + '_rmsd_' + str(round(r[1], 3)) +  r[0][0].getTitle(), r[0][0])
    pr.writePDB(workdir + str(count) + '_rmsd_' + str(round(r[1], 3)) +  r[2].getTitle(), r[2])
    count += 1