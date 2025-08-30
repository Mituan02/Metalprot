
import os
import sys
import prody as pr
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors

from rvrsprot.basic import constant
from rvrsprot.combs import gvdm_helper
from rvrsprot.combs import __search_cg_vdms as search_cg_vdms

pd.set_option("display.max_columns", None)


vdm_cg_aa_atommap_dict = {
    (0, 0):{
        'cg' : 'coo',
        'lgd_sel' : ['O0'],
        'represent_name' : 'OD1',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD1'],
    },    
        (0, 1):{
        'cg' : 'coo',
        'lgd_sel' : ['O1'],
        'represent_name' : 'OD1',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD1'],
    },    
        (0, 2):{
        'cg' : 'coo',
        'lgd_sel' : ['O2'],
        'represent_name' : 'OD1',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD1'],
    },
    (0, 3):{
        'cg' : 'coo',
        'lgd_sel' : ['O0'],
        'represent_name' : 'OD1',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD1'],
    },    
        (0, 4):{
        'cg' : 'coo',
        'lgd_sel' : ['O1'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD2'],
    },    
        (0, 5):{
        'cg' : 'coo',
        'lgd_sel' : ['O2'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD2'],
    },       
    (1, 0):{
        'cg' : 'coo',
        'lgd_sel' : ['O0'],
        'represent_name' : 'OE1',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE1'],
    },    
        (1, 1):{
        'cg' : 'coo',
        'lgd_sel' : ['O1'],
        'represent_name' : 'OE1',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE1'],
    },    
        (1, 2):{
        'cg' : 'coo',
        'lgd_sel' : ['O2'],
        'represent_name' : 'OE1',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE1'],
    },
    (1, 3):{
        'cg' : 'coo',
        'lgd_sel' : ['O0'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE2'],
    },    
        (1, 4):{
        'cg' : 'coo',
        'lgd_sel' : ['O1'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE2'],
    },    
        (1, 5):{
        'cg' : 'coo',
        'lgd_sel' : ['O2'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE2'],
    },       
    
        (2, 1):{
        'cg' : 'coh',
        'lgd_sel' : ['O1'],
        'represent_name' : 'OG',
        'correspond_resname' : 'SER',
        'correspond_names' : ['OG'],
    },    
        (2, 2):{
        'cg' : 'coh',
        'lgd_sel' : ['O2'],
        'represent_name' : 'OG',
        'correspond_resname' : 'SER',
        'correspond_names' : ['OG'],
    },

        (3, 1):{
        'cg' : 'coh',
        'lgd_sel' : ['O1'],
        'represent_name' : 'OG1',
        'correspond_resname' : 'THR',
        'correspond_names' : ['OG1'],
    },    
        (3, 2):{
        'cg' : 'coh',
        'lgd_sel' : ['O2'],
        'represent_name' : 'OG1',
        'correspond_resname' : 'THR',
        'correspond_names' : ['OG1'],
    },
            (4, 0):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['O0'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['O'],
    },    
        (4, 1):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['O1'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['O'],
    },
            (5, 1):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['O0'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['O'],
    },    
        (5, 1):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['O1'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['O'],
    }
}


vdm_cg_aa_atommap_dict = {
    (0, 0):{
        'cg' : 'coo',
        'lgd_sel' : ['O0'],
        'represent_name' : 'OD1',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD1'],
    },    

    (0, 3):{
        'cg' : 'coo',
        'lgd_sel' : ['O0'],
        'represent_name' : 'OD1',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD1'],
    },    

    (1, 0):{
        'cg' : 'coo',
        'lgd_sel' : ['O0'],
        'represent_name' : 'OE1',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE1'],
    },    
 
    (1, 3):{
        'cg' : 'coo',
        'lgd_sel' : ['O0'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE2'],
    },    
 
    
    (2, 1):{
        'cg' : 'coh',
        'lgd_sel' : ['O1'],
        'represent_name' : 'OG',
        'correspond_resname' : 'SER',
        'correspond_names' : ['OG'],
    },    

    (3, 1):{
        'cg' : 'coh',
        'lgd_sel' : ['O1'],
        'represent_name' : 'OG1',
        'correspond_resname' : 'THR',
        'correspond_names' : ['OG1'],
    },    

    (4, 0):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['O0'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['O'],
    },    
    (5, 1):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['O0'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['O'],
    }
}


def search_vdm2(dfa, ligand_coords, cg_id, input_dict, labels_cgs, df_cgs, dist_ind_cgs, rmsd = 0.5):
    '''
    s: combs2.design._sample.Sample

    cg_id: (0, 0)  #The first is to record the cg in ligand, the second is to record the cg in vdms. 

    input_dict[cg_id] = {

        cg : 'coo'  #dfa = s.cg_dict['coo']

        lgd_sel : ['C9', 'O3', 'O4']  # correpond to coo: ASP ['CG', 'OD1', 'OD2'] 

        represent_name : 'OD2'

        correspond_resname : 'ASP'

        correspond_names : ['CG', 'OD1', 'OD2']
    }
    '''

    #dfa = cg_dict[input_dict[cg_id]['cg']]

    #ligand_coords = get_ligand_coords(ligands, input_dict[cg_id]['lgd_sel'])

    labels, vdm_coords = search_cg_vdms.get_vdm_labels_coords(dfa, input_dict[cg_id]['represent_name'], input_dict[cg_id]['correspond_resname'], input_dict[cg_id]['correspond_names'])

    if vdm_coords.shape[0] <=0:
        print('No vdM coords are generated.')
        return
    
    radius = rmsd

    dists, inds = search_cg_vdms.get_nearest_vdms_rmsd(vdm_coords, ligand_coords, radius = radius)

    labels_cgs[cg_id] = labels
    df_cgs[cg_id] = dfa
    dist_ind_cgs[cg_id] = (dists, inds)
    print('{} Found {} vdm'.format(cg_id, sum([len(x) for x in inds])))
    return 


def vdm_ligand_clash(vdm, ligand, clash_radius = 1.2):
    '''
    The vdm X sidechain may clash with the ligand.
    return True if clash.
    '''
    vdm_coords =[]
    for i in range(vdm[(vdm['chain'] == 'X')].shape[0]):
        k = vdm[(vdm['chain'] == 'X')].iloc[i]
        if k['name'][0] == 'H':
            continue
        vdm_coords.append(k[['c_x', 'c_y', 'c_z']].values)

    ligand_coords = ligand.select('heavy').getCoords()

    if len(vdm_coords) == 0:
        print(vdm[['CG', 'rota', 'probe_name', 'chain', 'name']])
        return False
    nbr = NearestNeighbors(radius=clash_radius).fit(np.array(vdm_coords))
    adj_matrix = nbr.radius_neighbors_graph(ligand_coords).astype(bool)

    if adj_matrix.sum() > 0:
        return True
    return False


def vdm_target_clash(vdm, target, resnum, clash_radius = 1.2):
    '''
    The vdm X sidechain may clash with the ligand.
    return True if clash.
    '''
    vdm_coords =[]
    for i in range(vdm[(vdm['chain'] == 'X')].shape[0]):
        k = vdm[(vdm['chain'] == 'X')].iloc[i]
        if k['name'][0] == 'H':
            continue
        vdm_coords.append(k[['c_x', 'c_y', 'c_z']].values)
    sel = 'name N C CA O and (not resnum ' + str(resnum - 1) + ' ' + str(resnum) + ' ' + str(resnum + 1) + ')'
    taget_coords = target.select(sel).getCoords()

    if len(vdm_coords) == 0:
        print(vdm[['CG', 'rota', 'probe_name', 'chain', 'name']])
        return False
    nbr = NearestNeighbors(radius=clash_radius).fit(np.array(vdm_coords))
    adj_matrix = nbr.radius_neighbors_graph(taget_coords).astype(bool)

    if adj_matrix.sum() > 0:
        return True
    return False


def construct_vdm_write2(outdir, ligand, target, resnum, labels_cgs, df_cgs, dist_ind_cgs, clash_radius = 1.2):
    '''
    dist_ind_cgs: dict. {cg: (dists, inds)}, where dists in shape: (len(ligands), )
    df_cgs: {cg: df}
    '''

    cgCombInfo = search_cg_vdms.CgCombInfo()
    cgCombInfo.ligand_id = 0

    for cg in dist_ind_cgs.keys():
        info = []

        labels = labels_cgs[cg]
        dfa = df_cgs[cg]

        #dists = dist_ind_cgs[cg][0]
        inds = dist_ind_cgs[cg][1]

        ids = set()
        for _x in range(len(inds)):
            for _y in range(len(inds[_x])):
                ids.add(inds[_x][_y])

        for id in ids:
            x = labels.iloc[id]
            v = dfa[(dfa['CG'] == x['CG']) & (dfa['rota'] == x['rota']) & (dfa['probe_name'] == x['probe_name'])]

            #Filter based on C score.
            if v['C_score_bb_ind'].iloc[0] < 0:
                continue

            
            if vdm_ligand_clash(v, ligand, clash_radius):
                print('clash')
                continue

            # if vdm_target_clash(v, target, resnum, clash_radius):
            #     print('clash')
            #     continue

            info.append((round(1, 2), round(v['C_score_bb_ind'].iloc[0], 2), v)) 
            prefix = 'K' + '-'.join(str(z) for z in cg) + '_C_' + str(round(v['C_score_bb_ind'].iloc[0], 2)) + '_'                              
            #combs2.design.functions.print_dataframe(v, outpath=outdir, tag = '_' + str(ind), prefix = prefix)
            _label = str(x['CG']) + '_' + str(x['rota']) + '_' + str(x['probe_name'])
            ag = gvdm_helper.df2ag(v)
            pr.writePDB(outdir + prefix + '_' + _label,  ag)

        cgCombInfo.vdm_cgs[cg] = info

    return cgCombInfo


def run_vdm_sample(target, resnum, aa, ligand, ligand_coords, path_to_vdm_database, vdm_cg_aa_atommap_dict, outdir):

    labels_cgs = {}
    df_cgs = {}
    dist_ind_cgs = {}

    cg_ids = vdm_cg_aa_atommap_dict.keys()
    for cg_id in cg_ids:
        #Load vdm database.
        #Here we specify to check ASP vdMs. Note this can be programed to be easier.
        df_vdm = gvdm_helper.load_old_vdm(path_to_vdm_database, vdm_cg_aa_atommap_dict[cg_id]['cg'], aa)
        
        #Transformation.
        pos = target.select('name N CA C and resnum ' + str(resnum))
        tf = pr.calcTransformation(constant.ideal_ala_coords, pos.getCoords())
        df_vdm.loc[:, ['c_x', 'c_y', 'c_z']] = tf.apply(df_vdm[['c_x', 'c_y', 'c_z']].to_numpy())

        #Test if the first vdm's coords are changed.
        #x = labels_cgs[cg_id].iloc[0]
        #v = df_vdm[(df_vdm['CG'] == x['CG']) & (df_vdm['rota'] == x['rota']) & (df_vdm['probe_name'] == x['probe_name'])]
        #ag = gvdm_helper.df2ag(v, 'test_position_of_1st_vdm')
        #pr.writePDB(outdir + ag.getTitle(), ag)

        #Search vdms.
        search_vdm2(df_vdm, ligand_coords, cg_id, vdm_cg_aa_atommap_dict, labels_cgs, df_cgs, dist_ind_cgs, rmsd = 0.3)
        

    CgCombInfoDict = construct_vdm_write2(outdir, ligand, target, resnum, labels_cgs, df_cgs, dist_ind_cgs, clash_radius = 1.5)

    #search_cg_vdms.write_summary(outdir, CgCombInfoDict, name = '_summary.tsv')
    return


def sample_vdMs(target, ligand_trim, ligand_coords, path_to_vdm_database, outdir):
    resnums = [66, 67, 69, 70, 73, 74, 77, 106, 107, 110, 113, 114]
    aas = ['SER', 'THR', 'TYR', 'ASP', 'GLU', 'ASN', 'GLN', 'ARG', 'LYS', 'HIS']
    #resnums = [110]
    #aas = ['SER', 'GLU']

    for i in range(len(resnums)):
        resnum = resnums[i]
        for aa in aas:
            _outdir = outdir + 'v_' + str(resnum) + '_' + aa + '_'
            
            run_vdm_sample(target, resnum, aa, ligand_trim, ligand_coords, path_to_vdm_database, vdm_cg_aa_atommap_dict, _outdir)
    return



#####################################################

workdir = '/Users/lonelu/DesignData/APEX/Design_210/'


ro1 = ['FE', 'NA']
ro1_range = range(-40, 40, 8) 
ro2 = ['FE', 'NE2']
ro2_range = range(0, 360, 20) 

target_path = workdir + 'kh210_af3_rf.pdb'
target = pr.parsePDB(target_path)
tar_sel = 'resnum 50 and name CG ND1 CE1 CD2 NE2'

lig_path = '/Users/lonelu/DesignData/APEX/Design_210/1a2f_hem_oo.pdb'
lig = pr.parsePDB(lig_path)

# all_ligs = generate_rotated_porphyrins(lig, ro1, ro2, ro1_range, ro2_range)

filter_ligs = sample_porphyrin_pos(target, tar_sel, lig, ro1, ro2, ro1_range, ro2_range, dist =1.8, lig_sel = 'chain B or chain C')

os.makedirs(workdir + 'filtered_ligs/', exist_ok=True)
for lg in filter_ligs:
    pr.writePDB(workdir + 'filtered_ligs/' +  lg.getTitle(), lg)


outdir = workdir + 'output_240823/'
os.makedirs(outdir, exist_ok = True)

ligand_coords = search_cg_vdms.get_ligand_coords(filter_ligs, lgd_sel = ['O0','O1','O2'])
ligand_coords = np.reshape(ligand_coords, (-1,3))

ligand_trim = pr.parsePDB('/Users/lonelu/DesignData/APEX/Design_210/hem_trim_af3.pdb')

path_to_vdm_database='/Users/lonelu/DesignData/Combs/vdMs/'
sample_vdMs(target, ligand_trim, ligand_coords, path_to_vdm_database, outdir)