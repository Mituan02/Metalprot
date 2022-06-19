'''
Sample the lig positions based on vdm cg alignment.
In the example here, the belinostat has a benzyl ring ing the middle, which is a binding hot spot.
The script  here is to sample the positions by the benzyl ring.
'''

import os
import sys
import prody as pr
import numpy as np
from metalprot.basic import constant, utils
from metalprot.combs import search_lig_indep
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import cdist

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/belinostat/pose_lig_by_vdm.py
'''

def _ligcoord_on_vdms(lig, df_vdm, vdm_cg_aa_atommap_dict, cg_id):
    '''
    The vdms are all aligned on ideal_aa.
    The ligands here are supperimposed on the ideal_aa to generate the ligand special poses.
    '''

    labels, vdm_coords = search_lig_indep.get_vdm_labels_coords_4old_vdm_db(df_vdm, vdm_cg_aa_atommap_dict[cg_id]['correspond_resname'], vdm_cg_aa_atommap_dict[cg_id]['correspond_names'])
    lig_sel_coord = search_lig_indep.lgd_sel_coord(lig, vdm_cg_aa_atommap_dict[cg_id]['lgd_sel'])
    lig_sel_coord_r = np.array(lig_sel_coord).reshape((-1, 3))
    ligs_coords = []
    for i in range(len(labels)):
        vdm_coord = vdm_coords[i].reshape((-1,3))
        lig_cp = lig.copy()
        pr.calcTransformation(lig_sel_coord_r, vdm_coord).apply(lig_cp)        
        ligs_coords.append(lig_cp.select('heavy').getCoords())
         
    return labels, np.array(ligs_coords)


def _lig_bb_clashing(target, lig_coords, dist_cut=2.5):
    '''
    return any coords that are not clashing with the target.
    '''
    lig_atom_num = lig_coords.shape[0]
    lig_coords_r = lig_coords.reshape((-1, 3))
    bb_coords = target.select('heavy').getCoords()
    nn = NearestNeighbors(radius= dist_cut).fit(bb_coords)
    x_in_y = nn.radius_neighbors_graph(lig_coords_r).toarray().astype(bool)

    clashing0 = np.any(x_in_y, axis=1)
    clashing1 = np.any(clashing0.reshape((lig_atom_num, -1)), axis=1)
    #np.sum(clashing) #The clashing is too less to be considered.
    return clashing1   


def _vdm_bb_clashing(target, chidres, df_vdm_filter, dist_cut = 2.5):
    '''
    
    '''
    df_vdm_filter_heavy = df_vdm_filter[(df_vdm_filter['chain'] == 'X') 
            & (df_vdm_filter['resnum'] == 10)
            & ~(df_vdm_filter['name'].isin(['N', 'C', 'CA', 'O']))
            & ~(df_vdm_filter['atom_type_label'].isin(['h_pol', 'h_alkyl', 'h_aro']))]

    bb_coords = target.select('heavy and not ( chid ' + chidres[0] + ' and resnum ' + str(chidres[1]) + ')').getCoords()

    nn = NearestNeighbors(radius= dist_cut).fit(bb_coords)
    x_in_y = nn.radius_neighbors_graph(df_vdm_filter_heavy[['c_x', 'c_y', 'c_z']]).toarray().astype(bool)
    clashing = np.any(x_in_y, axis=1)

    labels_clashing = df_vdm_filter_heavy[['CG', 'rota', 'probe_name']][clashing]

    return labels_clashing


def _lig_vdm_clashing(df_vdm_filter, labels_filter, lig_coords_filter, dist_cut = 2.5):
    '''
    
    '''
    kk = np.ones(labels_filter.shape[0], dtype=bool)
    vdms = []
    for i in range(labels_filter.shape[0]):
        x = labels_filter.iloc[i]
        v = df_vdm_filter[(df_vdm_filter['CG'] == x['CG']) 
        & (df_vdm_filter['rota'] == x['rota']) 
        & (df_vdm_filter['probe_name'] == x['probe_name'])
        & (df_vdm_filter['resnum'] == 10)
        & (df_vdm_filter['chain'] == 'X')
        & ~(df_vdm_filter['atom_type_label'].isin(['h_pol', 'h_alkyl', 'h_aro']))]
        
        dists = cdist(lig_coords_filter[i], v[['c_x', 'c_y', 'c_z']].to_numpy())
        if (dists < dist_cut).any():
            kk[i] = False
        else:
            vdms.append(v)
    return vdms, lig_coords_filter[kk]


def _ligcoord_on_vdm_on_bb(target, cg_id, aa, chidres, df_vdm_filter, labels, lig_coords):
    '''
    
    '''
    #Get transformation on bb.
    pos = target.select('name N CA C and chid ' + chidres[0] + ' and resnum ' + str(chidres[1]))
    tf = pr.calcTransformation(constant.ideal_ala_coords, pos.getCoords())
    shape = lig_coords.shape
    lig_coords_new = tf.apply(lig_coords.reshape((-1, 3))).reshape(shape)
    df_vdm_filter.loc[:, ['c_x', 'c_y', 'c_z']] = tf.apply(df_vdm_filter[['c_x', 'c_y', 'c_z']].to_numpy())

    #Clashing the ligs on bb
    clashing = _lig_bb_clashing(target, lig_coords_new, dist_cut=2.5)

    #Get vdm score filtered labels index.
    labels_prefilter = df_vdm_filter[['CG', 'rota', 'probe_name']].drop_duplicates().astype(str).T.agg(','.join)
    labels_prefilter_ind = labels.astype(str).T.agg(','.join).isin(labels_prefilter)

    #Get vdm sc clashing filtered labels index.
    labels_clashing = _vdm_bb_clashing(target, chidres, df_vdm_filter, dist_cut = 2.5)
    labels_clashing = labels_clashing.drop_duplicates().astype(str).T.agg(','.join)
    labels_clashfilter_ind = ~labels.astype(str).T.agg(','.join).isin(labels_clashing)

    #get filtered ligcoords
    filter_comb = np.invert(clashing) & labels_prefilter_ind & labels_clashfilter_ind
    labels_filter = labels[filter_comb]
    lig_coords_filter = lig_coords_new[filter_comb]

    #Filter by vdm lig clashing
    vdms_filter, lig_coords_filter = _lig_vdm_clashing(df_vdm_filter, labels_filter, lig_coords_filter, dist_cut = 2.5)
    
    return (cg_id, aa, chidres), (vdms_filter, lig_coords_filter)


def write_ligs(lig, cg_id, chidres, aa, vdms_filter, lig_coords_filter, outdir):
    '''
    
    '''
    for i in range(lig_coords_filter.shape[0]):
        v = vdms_filter[i]
        x = v[['CG', 'rota', 'probe_name']].iloc[0]
        ag = search_lig_indep.df2ag(v)

        title = chidres[0] + '_' + str(chidres[1]) + '_' + aa + '_' + cg_id + '_' + '_'.join([str(_x) for _x in x])
        coord = lig_coords_filter[i]
        lig_cp = lig.copy()
        lig_cp.setCoords(coord)

        pr.writePDB(outdir +  title + '_vdm', ag)
        pr.writePDB(outdir + title + '_lig', lig_cp)
    return 


def all_ligcoord_on_vdm_on_bb(target, lig, para, outdir):
    '''
    The input data structure:
        para.vdm_cg_aa_atommap_dict = {
            ('ph_0'):{
                'cg' : 'ph',
                'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
                'correspond_resname' : 'PHE',
                'correspond_names' : ['CG', 'CD1', 'CD2', 'CZ'],
                'aas' : 'FWY',
                'filter_hb' : False,
                'filter_cc' : True
            },
        }

    '''
    #lig_coord_on_vdms_pos_dict = {}
    abples, phipsi = utils.seq_get_ABPLE(target)

    for cg_id in para.vdm_cg_aa_atommap_dict.keys():
        for a in para.vdm_cg_aa_atommap_dict[cg_id]['aas']:
            aa = constant.inv_one_letter_code[a]
            # Load vdm
            df_vdm = search_lig_indep.load_old_vdm(para.path_to_database, para.vdm_cg_aa_atommap_dict[cg_id]['cg'], aa)
            # Lig coord superimposed on vdm
            labels, ligs_coords = _ligcoord_on_vdms(lig, df_vdm, para.vdm_cg_aa_atommap_dict, cg_id)

            for chidres in para.predefined_resnums:
                print('Searching: ' + cg_id + ' ' + aa + ' ' + chidres[0] + str(chidres[1]))

                pos = target.select('chid ' + chidres[0] + ' and resnum ' + str(chidres[1]))
                abple = abples[pos.getResindices()[0]]

                # Prefilter by score and ABPLE
                df_vdm_filter = search_lig_indep.filter_db(df_vdm, use_enriched = para.use_enriched, use_abple=para.use_abple, abple=abple)        
                
                # get the lig_coords aligned on bb and filter by vdm clashing and lig clashing. 
                cgid_aa_chidres, lablesfilter_ligcoordsfilter_tf = _ligcoord_on_vdm_on_bb(target, cg_id, aa, chidres, df_vdm_filter, labels, ligs_coords)
                if len(lablesfilter_ligcoordsfilter_tf[0]) > 0:
                    print('Find {} candidates.'.format(len(lablesfilter_ligcoordsfilter_tf[0])))
                    #lig_coord_on_vdms_pos_dict[cgid_aa_chidres] = lablesfilter_ligcoordsfilter_tf
                    write_ligs(lig, cg_id, chidres, aa, lablesfilter_ligcoordsfilter_tf[0], lablesfilter_ligcoordsfilter_tf[1], outdir)
                else:
                    print('Find None candidates.')
    #return lig_coord_on_vdms_pos_dict
    return


class Para:
    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

    resnums = [3, 7, 10, 14, 17, 18, 21, 24, 25, 
        51, 54, 58, 61, 65, 68, 69, 72, 77, 81, 84, 88, 91, 92, 95, 99, 
        125, 128, 132, 135, 139, 142, 146]
    #resnums = [61]
    predefined_resnums = [('A', r) for r in resnums]

    use_enriched = True
    use_abple=True

    rmsd = 0.6

    vdm_cg_aa_atommap_dict = {
        ('ph_0'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CG', 'CD1', 'CD2', 'CZ'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_1'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CG', 'CD2', 'CD1', 'CZ'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_2'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CD1', 'CG', 'CE1', 'CE2'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_3'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CD1', 'CE1', 'CG', 'CE2'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_4'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CE1', 'CD1', 'CZ', 'CD2'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_5'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CE1', 'CZ', 'CD1', 'CD2'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_6'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CZ', 'CE1', 'CE2', 'CG'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_7'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CZ', 'CE2', 'CE1', 'CG'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_8'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CE2', 'CZ', 'CD2', 'CD1'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_9'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CE2', 'CD2', 'CZ', 'CD1'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_10'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CD2', 'CE2', 'CG', 'CE1'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        },
        ('ph_11'):{
            'cg' : 'ph',
            'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
            'correspond_resname' : 'PHE',
            'correspond_names' : ['CD2', 'CG', 'CE2', 'CE1'],
            'aas' : 'FWY',
            'filter_hb' : False,
            'filter_cc' : True
        }
    }


def main():

    workdir = '/mnt/e/DesignData/Metalloenzyme/'

    lig = pr.parsePDB(workdir + 'ligs/50g.pdb')

    target = pr.parsePDB(workdir + 'targets/01_f63440_nick_ala.pdb')

    #Para
    para = Para()

    outdir = workdir + 'outdir/'
    os.makedirs(outdir, exist_ok = True)

    #run the search.
    all_ligcoord_on_vdm_on_bb(target, lig, para, outdir)

    #lig_coord_on_vdms_pos_dict = all_ligcoord_on_vdm_on_bb(target, lig, para, outdir)

    #write the results.
    # for key in lig_coord_on_vdms_pos_dict.keys():
    #     cg_id, aa, chidres = key 
    #     vdms_filter, lig_coords_filter = lig_coord_on_vdms_pos_dict[key]

    #     write_ligs(lig, cg_id, chidres, aa, vdms_filter, lig_coords_filter, outdir)

    return


if __name__=='__main__':
    main()


'''
import os
import sys
import prody as pr
import numpy as np
from metalprot.basic import constant, utils, transformation
from metalprot.combs import search_lig_indep, position_ligand
from sklearn.neighbors import NearestNeighbors

workdir = '/mnt/e/DesignData/Metalloenzyme/'
lig = pr.parsePDB(workdir + 'ligs/50g.pdb')
target = pr.parsePDB(workdir + 'targets/01_f63440_nick_ala.pdb')

cg ='ph'
cg_id = 'ph_0'
aa = 'PHE'
chidres = ('A', 61)
para = Para()
abple = 'A'

df_vdm = search_lig_indep.load_old_vdm(para.path_to_database, cg, aa)

df_vdm_filter = search_lig_indep.filter_db(df_vdm, use_enriched = para.use_enriched, use_abple=para.use_abple, abple=abple)        
                
labels, lig_coords = _ligcoord_on_vdms(lig, df_vdm, para.vdm_cg_aa_atommap_dict, cg_id)

pos = target.select('name N CA C and chid ' + chidres[0] + ' and resnum ' + str(chidres[1]))

    #Get transformation on bb.
    shape = lig_coords.shape
    tf = pr.calcTransformation(constant.ideal_ala_coords, pos.getCoords())
    lig_coords_new = tf.apply(lig_coords.reshape((-1, 3))).reshape(shape)
    df_vdm_filter.loc[:, ('c_x', 'c_y', 'c_z')] = tf.apply(df_vdm_filter[['c_x', 'c_y', 'c_z']].to_numpy())

    #Clashing the ligs on bb
    clashing = _lig_bb_clashing(target, lig_coords_new, dist_cut=2.5)

    #Get vdm score filtered labels index.
    labels_prefilter = df_vdm_filter[['CG', 'rota', 'probe_name']].drop_duplicates().astype(str).T.agg(','.join)
    labels_prefilter_ind = labels.astype(str).T.agg(','.join).isin(labels_prefilter)

    #Get vdm sc clashing filtered labels index.
    labels_clashing = _vdm_bb_clashing(target, chidres, df_vdm_filter, dist_cut = 2.5)
    labels_clashing = labels_clashing.drop_duplicates().astype(str).T.agg(','.join)
    labels_clashfilter_ind = ~labels.astype(str).T.agg(','.join).isin(labels_clashing)

    #get filtered ligcoords
    filter_comb = np.invert(clashing) & labels_prefilter_ind & labels_clashfilter_ind
    labels_filter = labels[filter_comb]
    lig_coords_filter = lig_coords_new[filter_comb]

    #Filter by vdm lig clashing
    vdms_filter, lig_coords_filter = _lig_vdm_clashing(df_vdm_filter, labels_filter, lig_coords_filter, dist_cut = 2.5)

outdir = workdir + 'test/'
os.makedirs(outdir, exist_ok = True)

write_ligs(lig, cg_id, chidres, lig_coords_filter, df_vdm, vdms_filter, outdir)
'''