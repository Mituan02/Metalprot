import numpy as np
import prody as pr
from prody.measure.transform import calcRMSD
from scipy.spatial.distance import cdist
import itertools
from sklearn.neighbors import NearestNeighbors
from .vdmer import pair_wise_geometry_matrix

class Search_filter:
    def __init__(self, filter_abple = False, filter_phipsi = True, max_phipsi_val = 15, 
    filter_vdm_score = False, min_vdm_score = 0, filter_vdm_count = False, min_vdm_clu_num = 20,
    after_search_filter_geometry = False, filter_based_geometry_structure = False, angle_tol = 5, aa_aa_tol = 0.5, aa_metal_tol = 0.2,
    pair_angle_range = None, pair_aa_aa_dist_range = None, pair_metal_aa_dist_range = None,
    after_search_filter_qt_clash = False, vdm_vdm_clash_dist = 2.7, vdm_bb_clash_dist = 2.2, 
    after_search_open_site_clash = True, open_site_dist = 3.0, write_filtered_result = False, 
    selfcenter_filter_member_phipsi = True):

        self.filter_abple = filter_abple
        self.filter_phipsi = filter_phipsi
        self.max_phipsi_val = max_phipsi_val

        self.filter_vdm_score = filter_vdm_score
        self.min_vdm_score = min_vdm_score
        self.filter_vdm_count = filter_vdm_count
        self.min_vdm_clu_num = min_vdm_clu_num

        
        self.after_search_filter_geometry = after_search_filter_geometry  
        self.filter_based_geometry_structure = filter_based_geometry_structure
        self.angle_tol = angle_tol 
        self.aa_aa_tol = aa_aa_tol
        self.aa_metal_tol = aa_metal_tol

        self.pair_angle_range = pair_angle_range # [90, 110]
        self.pair_aa_aa_dist_range = pair_aa_aa_dist_range # [2.8, 3.4]
        self.pair_metal_aa_dist_range = pair_metal_aa_dist_range # [2.0, 2.3]

        self.after_search_filter_qt_clash = after_search_filter_qt_clash 
        self.vdm_vdm_clash_dist = vdm_vdm_clash_dist
        self.vdm_bb_clash_dist = vdm_bb_clash_dist
        self.after_search_open_site_clash = after_search_open_site_clash
        self.open_site_dist = open_site_dist
        self.write_filtered_result = write_filtered_result

        self.selfcenter_filter_member_phipsi = selfcenter_filter_member_phipsi

    def para2string(self):
        parameters = "Filter parameters: \n"
        parameters += 'filter_abple: ' + str(self.filter_abple) + ' \n'
        parameters += 'filter_phipsi: ' + str(self.filter_phipsi) + ' \n'
        parameters += 'max_phipsi_val: ' + str(self.max_phipsi_val) + ' \n'

        parameters += 'filter_vdm_score: ' + str(self.filter_vdm_score) + ' \n'
        parameters += 'min_vdm_score: ' + str(self.min_vdm_score) + ' \n'
        parameters += 'filter_vdm_count: ' + str(self.filter_vdm_count) + ' \n'
        parameters += 'min_vdm_clu_num: ' + str(self.min_vdm_clu_num) + ' \n'

        parameters += 'after_search_filter_geometry: ' + str(self.after_search_filter_geometry) + ' \n'
        parameters += 'filter_based_geometry_structure: ' + str(self.filter_based_geometry_structure) + ' \n'
        parameters += 'pair_angle_range: ' + str(self.pair_angle_range) + ' \n'
        parameters += 'pair_aa_aa_dist_range: ' + str(self.pair_aa_aa_dist_range) + ' \n'
        parameters += 'pair_metal_aa_dist_range: ' + str(self.pair_metal_aa_dist_range) + ' \n'
        parameters += 'filter_qt_clash: ' + str(self.after_search_filter_qt_clash) + ' \n'
        parameters += 'vdm_vdm_clash_dist: ' + str(self.vdm_vdm_clash_dist) + ' \n'
        parameters += 'vdm_bb_clash_dist: ' + str(self.vdm_bb_clash_dist) + ' \n'
        parameters += 'after_search_open_site_clash: ' + str(self.after_search_open_site_clash) + ' \n'
        parameters += 'open_site_dist: ' + str(self.open_site_dist) + ' \n'
        parameters += 'write_filtered_result: ' + str(self.write_filtered_result) + ' \n'
        
        parameters += 'selfcenter_filter_member_phipsi: ' + str(self.selfcenter_filter_member_phipsi) + ' \n'
        return parameters


    @staticmethod
    def after_search_geo_pairwise_satisfied(combinfo, pair_angle_range = None, pair_aa_aa_dist_range = None, pair_metal_aa_dist_range = None):
        '''
        range = (75, 125) for Zn.
        if all pairwise angle is between the range. The geometry is satisfied.
        '''
        satisfied = True
        if pair_angle_range:
            for an in combinfo.angle_pair:
                if an < pair_angle_range[0] or an > pair_angle_range[1]:
                    combinfo.pair_angle_ok = -1
                    satisfied = False
                    break

        if pair_aa_aa_dist_range:           
            for ad in combinfo.aa_aa_pair:
                if ad < pair_aa_aa_dist_range[0] or ad > pair_aa_aa_dist_range[1]:
                    combinfo.pair_aa_aa_dist_ok = -1
                    satisfied = False
                    break

        if pair_metal_aa_dist_range:
            combinfo.pair_aa_metal_dist_ok = 1
            for amd in combinfo.metal_aa_pair:
                if amd < pair_metal_aa_dist_range[0] or amd > pair_metal_aa_dist_range[1]:
                    combinfo.pair_aa_metal_dist_ok = -1
                    satisfied = False
                    break
                    

        return satisfied


    @staticmethod
    def get_min_geo(geometry, geo_struct, metal_sel = 'name NI MN ZN CO CU MG FE' ):
        '''
        Metal must be the last atom in the prody object.
        '''
        aa_coords = geo_struct.select('not ' + metal_sel).getCoords()
        metal_coord = geo_struct.select(metal_sel).getCoords()[0]
        ct_len = len(aa_coords)

        min_rmsd = 0
        min_geo_struct = None
        for xs in itertools.permutations(range(ct_len), ct_len):
            _geo_struct = geo_struct.copy()
            coords = []
            for x in xs:
                coords.append(aa_coords[x])
            coords.append(metal_coord) 
            _geo_struct.setCoords(np.array(coords))
            try:
                pr.calcTransformation(_geo_struct.select('not oxygen'), geometry).apply(_geo_struct)
                
                rmsd = pr.calcRMSD(_geo_struct.select('not oxygen'), geometry)
                if not min_geo_struct:
                    min_geo_struct = _geo_struct
                    min_rmsd = rmsd
                elif rmsd < min_rmsd:
                    min_geo_struct = _geo_struct
                    min_rmsd = rmsd
            except Exception as e:
                print(f"{type(e).__name__} at line {e.__traceback__.tb_lineno} of {__file__}: {e}")
                
        return min_geo_struct, min_rmsd
        

    @staticmethod
    def after_search_geo_strcut_satisfied(comb_info, min_geo_struct, angle_tol, aa_aa_tol, aa_metal_tol):

        aa_aa_pair, metal_aa_pair, angle_pair = pair_wise_geometry_matrix(min_geo_struct)

        info_aa_aa_pair, info_metal_aa_pair, info_angle_pair = pair_wise_geometry_matrix(comb_info.geometry)

        satisfied = True
        
        comb_info.pair_aa_metal_dist_ok = 1
        for i in range(len(metal_aa_pair)):
            if info_metal_aa_pair[i] < metal_aa_pair[i] - aa_metal_tol or info_metal_aa_pair[i] > metal_aa_pair[i] + aa_metal_tol:
                comb_info.pair_aa_metal_dist_ok = -1
                satisfied = False
                break

        comb_info.pair_aa_aa_dist_ok = 1
        for i, j in itertools.combinations(range(aa_aa_pair.shape[0]), 2):
            if info_aa_aa_pair[i, j] < aa_aa_pair[i, j] - aa_aa_tol or info_aa_aa_pair[i, j] > aa_aa_pair[i, j] + aa_aa_tol:
                comb_info.pair_aa_aa_dist_ok = -1
                satisfied = False
                break

        comb_info.pair_angle_ok = 1
        for i, j in itertools.combinations(range(aa_aa_pair.shape[0]), 2):
            if info_angle_pair[i, j] < angle_pair[i, j] - angle_tol or info_angle_pair[i, j] > angle_pair[i, j] + angle_tol:
                comb_info.pair_angle_ok = -1
                satisfied = False
                break

        return satisfied


    @staticmethod
    def vdm_clash(vdms, target, vdm_vdm_clash_dist = 2.7, vdm_bb_clash_dist = 2.2, unsupperimposed = True, wins = None, align_sel = 'name N CA C'):
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

        for i, j in itertools.combinations(range(len(coords)), 2):

            neigh_y = NearestNeighbors(radius= vdm_vdm_clash_dist)
            neigh_y.fit(coords[i])
            x_in_y = neigh_y.radius_neighbors(coords[j])
            x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])
            if x_has_y:
                return True
        
        bb_coord = target.select('protein and heavy and bb').getCoords()
        for i in range(len(coords)):
            neigh_y = NearestNeighbors(radius= vdm_bb_clash_dist)
            neigh_y.fit(bb_coord)
            x_in_y = neigh_y.radius_neighbors(coords[i])
            x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])
            if x_has_y:
                return True
        return False


    @staticmethod
    def open_site_clashing(vdms, target, ideal_geo, open_site_dist = 3.0):
        '''
        The open site of ideal_geo must be Oxygen, the other atom could not be Oxygen.

        If clash detected, return True.
        '''
        ideal_geo_coord = [ideal_geo.select('oxygen')[0].getCoords()]

        coords = []
        for i in range(len(vdms)):
            vdm = vdms[i]
            vdm_sel = 'protein and heavy and sc and not name CB'
            coord = vdm.query.select(vdm_sel).getCoords()
            coords.extend(coord)
        bb_coord = target.select('protein and heavy and bb').getCoords()
        coords.extend(bb_coord)

        neigh_y = NearestNeighbors(radius= open_site_dist)
        neigh_y.fit(coords)
        x_in_y = neigh_y.radius_neighbors(ideal_geo_coord)
        x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])
        if x_has_y:
            return True
        return False




