'''
Searching 2ndshell using general vdM database.
Please check 'scrips/position_ligand/1ukr/run_search_2ndshell.py'
'''
from ..basic import constant, utils
from ..combs import search_lig_indep, gvdm_helper, search_lig_indep_wrap, search_lig_indep_2ndshell

class Para:
    task_type = 'search_2ndshell'
    #>>> vdM database filter.
    use_enriched = True
    use_abple=True

    rmsd = 0.75
    
    aa_cg_2ndshell_dict = {
            'D':['coo'],
            'E':['coo'],
            'H':['hid', 'hie']
        }

    vdm_cg_aa_atommap_dict = {
        ('coo_0'):{
            'cg' : 'coo',
            'lgd_sel' : ['CB', 'CG', 'OD1', 'OD2'],
            'represent_name' : 'OD2',
            'correspond_resname' : 'ASP',
            'correspond_names' : ['CB', 'CG', 'OD1', 'OD2'],
            'aas' : 'GASTYNQDEWKRHC',
            #'aas': 'H',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('coo_1'):{
            'cg' : 'coo',
            'lgd_sel' : ['CG', 'CD', 'OE1', 'OE2'],
            'represent_name' : 'OE2',
            'correspond_resname' : 'GLU',
            'correspond_names' : ['CG', 'CD', 'OE1', 'OE2'],
            'aas' : 'GASTYNQDEWKRHC',
            #'aas': 'H',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('hid_0'):{
            'cg' : 'hid',
            'lgd_sel' : ['ND1', 'CE1', 'CG', 'CD2'],
            'represent_name' : 'ND1',
            'correspond_resname' : 'HIS',
            'correspond_names' : ['ND1', 'CE1', 'CG', 'CD2'],
            'aas' : 'GASTYNQDEHC',
            #'aas' : 'Q',
            'filter_hb' : False,
            'filter_cc' : False
        },
        ('hie_0'):{
            'cg' : 'hie',
            'lgd_sel' : ['NE2', 'CE1', 'CD2', 'CG'],
            'represent_name' : 'NE2',
            'correspond_resname' : 'HIS',
            'correspond_names' : ['NE2', 'CE1', 'CD2', 'CG'],
            'aas' : 'GASTYNQDEHC',
            #'aas' : 'Q',
            'filter_hb' : False,
            'filter_cc' : False
        }
    }


def prepare_mvdm_ligs(neighbor_comb_dict, aa):
    '''
    ss.neighbor_comb_dict
        (34, 60, 64):
            ((34, 60, 64), (('H', 2051), ('H', 1902), ('H', 930)))
                {34: [<metalprot.basic.vdmer.VDM at 0x7fb3131dc970>],
                 60: [<metalprot.basic.vdmer.VDM at 0x7fb3131dca30>],
                 64: [<metalprot.basic.vdmer.VDM at 0x7fb3131dcaf0>]}

    aa: 'H', 'D' or 'E'

    contact_aa_lig_dict
        (34, 'H', 1687):
            (<AtomGroup: AAMetalPhiPsi_HIS_cluster_1687_mem_26_centroid_2017_5xzj_ZN_6_mem0 (15 atoms)>,
             [((34, 60, 64), (('H', 1687), ('H', 1063), ('H', 400)))])
    '''
    contact_aa_lig_dict = {}
    for win_comb in neighbor_comb_dict.keys():
        comb_dict = neighbor_comb_dict[win_comb]
        for comb_key in comb_dict.keys():
            poses = comb_key[0]
            aa_vdmIds = comb_key[1]

            query_dict = comb_dict[comb_key].query_dict
            for i in range(len(poses)):
                if aa_vdmIds[i][0] != aa:
                    continue

                pos = poses[i]
                _key = (pos, aa_vdmIds[i][0], aa_vdmIds[i][1]) #>>> The _key looks like: (34, 'H', 2051)      
                if _key in contact_aa_lig_dict.keys():
                    contact_aa_lig_dict[_key][1].append(comb_key)
                    continue

                mvdm = query_dict[pos][0].query
                contact_aa_lig_dict[_key]= (mvdm, [comb_key])
    return contact_aa_lig_dict


def run_search_2ndshell(target, neighbor_comb_dict, target_ind2chidres, resnum_filtered, path_to_database, para):
    '''
    lig_vdm_dict:{
        54:
            [(<AtomGroup: Lig-54_hid_A58_GLN_rmsd_1.67_v_0.3_198 (33 atoms)>,
                ('hid', 'GLN', ('A', 58), 'hid_0', 54, 1.6684779654166868, 198, 0.2699825, None, False, True)), 
                #cg, aa, chidres, cg_id, lig_id, _rmsd, vdm_id, score, v, contact_hb, contact_cc
                #please check search_lig_indep.search_lig_at_cg_aa_resnum()
          ]
    }
    '''
    for aa in ['H', 'D', 'E']:

        contact_aa_lig_dict = prepare_mvdm_ligs(neighbor_comb_dict, aa)

        lig_keys = list(contact_aa_lig_dict.keys())
        if len(lig_keys) <= 0:
            continue

        ligs = [contact_aa_lig_dict[k][0] for k in contact_aa_lig_dict.keys()]

        predefined_chidres = search_lig_indep_2ndshell.calc_all_chidres_around_pose(target, resnum_filtered, dist = 12)

        lig_vdm_dict = search_lig_indep_wrap.run_search(target, ligs, path_to_database, para, predefined_chidres, (aa, para.aa_cg_2ndshell_dict[aa]))

        for lig_id in lig_vdm_dict.keys():
            lig_key = lig_keys[lig_id]
            comb_keys = contact_aa_lig_dict[lig_key][1]
            for comb_key in comb_keys:
                for win_comb in neighbor_comb_dict.keys():
                    _the_chid_res = [target_ind2chidres[x] for x in comb_key[0]]
                    for _vdm in lig_vdm_dict[lig_id]:
                        if _vdm[1][2] in _the_chid_res:
                            continue
                        if comb_key not in neighbor_comb_dict[win_comb].keys():
                            continue
                        if lig_key not in neighbor_comb_dict[win_comb][comb_key].secondshell_dict.keys():
                            neighbor_comb_dict[win_comb][comb_key].secondshell_dict[lig_key] = [_vdm]
                        else:
                            neighbor_comb_dict[win_comb][comb_key].secondshell_dict[lig_key].append(_vdm)
                        
    return

'''
#<<< Development
from metalprot.combs import search_lig_indep, gvdm_helper, search_lig_indep_wrap, search_lig_indep_2ndshell
path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

para = Para()
aa = 'H'
contact_aa_lig_dict = prepare_mvdm_ligs(ss.neighbor_comb_dict, aa)
lig_keys = list(contact_aa_lig_dict.keys())
ligs = [contact_aa_lig_dict[k][0] for k in contact_aa_lig_dict.keys()]

_resnum_filtered = [ ('A', 61)] 
predefined_chidres = search_lig_indep_2ndshell.calc_all_chidres_around_pose(ss.target, _resnum_filtered, dist = 15)
predefined_chidres = [('A', 58)]
lig_vdm_dict = search_lig_indep_wrap.run_search(ss.target, ligs, path_to_database, para, predefined_chidres, para.aa_cg_2ndshell_dict[aa])

for lig_id in lig_vdm_dict.keys():
    lig_key = lig_keys[lig_id]
    comb_keys = contact_aa_lig_dict[lig_key][1]
    for win_comb in ss.neighbor_comb_dict.keys():
        for comb_key in comb_keys:
            _the_chid_res = [ss.target_ind2chidres[x] for x in comb_key[0]]
            for _vdm in lig_vdm_dict[lig_id]:
                if _vdm[1][2] in _the_chid_res:
                    continue
                if _vdm[1][2] not in ss.neighbor_comb_dict[win_comb][comb_key].secondshell_dict.keys():
                    ss.neighbor_comb_dict[win_comb][comb_key].secondshell_dict[_vdm[1][2]] = [_vdm]
                else:
                    ss.neighbor_comb_dict[win_comb][comb_key].secondshell_dict[_vdm[1][2]].append(_vdm)

#>>>
'''
