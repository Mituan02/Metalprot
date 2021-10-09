import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import search_struct, extract_vdm, ligand_database
from metalprot import hull
'''
python /mnt/e/GitHub_Design/Metalprot/scripts_other/metal_position_analysis.py
'''

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/'


# Get query pdbs, Add the cluster metal points into the query.hull_points

querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['AAMetalPhiPsi_CYS'], score_cut = 1, clu_num_cut = 0)

align_sel = 'name N CA C'
query_all_metal = querys[0].copy()
clu_rank = 0
for query in querys:
    inds = np.unique(query.query.getResindices())
    ind = inds[int(inds.shape[0]/2)-1]
    transform = pr.calcTransformation(query.query.select(align_sel + ' and resindex ' + str(ind)), query_all_metal.query.select(align_sel + ' and resindex ' + str(query_all_metal.contact_resind)))
    transform.apply(query.query)

    clu = extract_vdm.get_vdm_cluster(query)
    clu.realign_by_CCAN(target = query, align_sel=align_sel)
    query.cluster = clu
    query.extract_mem_metal_point()

    for c in query.cluster.querys:
        outdir = query_dir + 'AAMetalPhiPsi_CYS_realignCenter/' + str(clu_rank) + '/'
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        pr.writePDB(outdir + 'align_' + c.query.getTitle(), c.query)
    clu_rank +=1

#-----------------------------------------

query_all_metal = querys[0].copy()
points = []
count = 0
for _query in querys:
    query = _query.copy()
    inds = np.unique(query.query.getResindices())
    ind = inds[int(inds.shape[0]/2)-1]
    #transform = pr.calcTransformation(query.query.select(align_sel + ' and resindex ' + str(ind)), query_all_metal.query.select(align_sel + ' and resindex ' + str(query_all_metal.contact_resind)))
    #transform.apply(query.query)
    #transform.apply(query.hull_ag)
    #points.extend(query.hull_ag.getCoords())
    points = query.hull_ag.getCoords()
    #query_all_metal.hull_ag = hull.transfer2pdb(points)
    hull.write2pymol(points, query_dir + 'AAMetalPhiPsi_CYS_coords/', 'align_' + str(count) + '_' + query_all_metal.query.getTitle())
    count += 1
