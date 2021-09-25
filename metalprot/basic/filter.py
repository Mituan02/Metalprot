import numpy as np
import prody as pr
from scipy.spatial.distance import cdist
import itertools
from sklearn.neighbors import NearestNeighbors

class Search_filter:
    def __init__(self, filter_abple = False, filter_phipsi = True, max_phipsi_val = 15, 
    filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
    after_search_filter = False, pair_angle_range = None, pair_aa_aa_dist_range = None, pair_metal_aa_dist_range = None,
    filter_qt_clash = False, write_filtered_result = False):

        self.filter_abple = filter_abple
        self.filter_phipsi = filter_phipsi
        self.max_phipsi_val = max_phipsi_val

        self.filter_vdm_score = filter_vdm_score
        self.min_vdm_score = min_vdm_score
        self.filter_vdm_count = filter_vdm_count
        self.min_vdm_clu_num = min_vdm_clu_num

        
        self.after_search_filter = after_search_filter       
        self.pair_angle_range = pair_angle_range # [90, 110]
        self.pair_aa_aa_dist_range = pair_aa_aa_dist_range # [2.8, 3.4]
        self.pair_metal_aa_dist_range = pair_metal_aa_dist_range # [2.0, 2.3]
        self.filter_qt_clash = filter_qt_clash
        self.write_filtered_result = write_filtered_result

    def para2string(self):
        parameters = "Filter parameters: \n"
        parameters += 'filter_abple: ' + str(self.filter_abple) + ' \n'
        parameters += 'filter_phipsi: ' + str(self.filter_phipsi) + ' \n'
        parameters += 'max_phipsi_val: ' + str(self.max_phipsi_val) + ' \n'

        parameters += 'filter_vdm_score: ' + str(self.filter_vdm_score) + ' \n'
        parameters += 'min_vdm_score: ' + str(self.min_vdm_score) + ' \n'
        parameters += 'filter_vdm_count: ' + str(self.filter_vdm_count) + ' \n'
        parameters += 'min_vdm_clu_num: ' + str(self.min_vdm_clu_num) + ' \n'

        parameters += 'after_search_filter: ' + str(self.after_search_filter) + ' \n'
        parameters += 'pair_angle_range: ' + str(self.pair_angle_range) + ' \n'
        parameters += 'pair_aa_aa_dist_range: ' + str(self.pair_aa_aa_dist_range) + ' \n'
        parameters += 'pair_metal_aa_dist_range: ' + str(self.pair_metal_aa_dist_range) + ' \n'
        parameters += 'filter_qt_clash: ' + str(self.filter_qt_clash) + ' \n'
        parameters += 'write_filtered_result: ' + str(self.write_filtered_result) + ' \n'
        return parameters


    @staticmethod
    def vdm_clash(vdms, target, clash_dist = 1.8, unsupperimposed = True, wins = None, align_sel = 'name N CA C'):
        '''
        clashing with sklearn.neighbors NearestNeighbors method.
        All sc except CB atom of vdm are used for clashing checking.
        All bb of target are used for clashing chekcing.
        
        If clash detected, return True.
        '''
        coords = []
        for i in range(len(vdms)):
            vdm = vdms[i]
            if unsupperimposed:
                win = wins[i]
                target_sel = 'resindex ' + str(win) + ' and ' + align_sel
                query_sel = 'resindex ' + str(vdm.contact_resind) + ' and '+ align_sel

                if len(vdm.query.select(query_sel)) != len(target.select(target_sel)):
                    print('supperimpose_target_bb not happening')
                    continue
                
                transform = pr.calcTransformation(vdm.query.select(query_sel), target.select(target_sel))
                transform.apply(vdm.query)
            
            vdm_sel = 'protein and heavy and sc and not name CB'

            coord = vdm.query.select(vdm_sel).getCoords()
            coords.append(coord)

        coord = target.select('protein and heavy and bb').getCoords()
        coords.append(coord)

        for i, j in itertools.combinations(range(len(coords)), 2):

            neigh_y = NearestNeighbors(radius= clash_dist)
            neigh_y.fit(coords[i])
            x_in_y = neigh_y.radius_neighbors(coords[j])
            x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])
            if x_has_y:
                return True
        
        return False






