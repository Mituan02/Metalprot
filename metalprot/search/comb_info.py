import prody as pr

from ..basic import quco
from ..basic import hull


class CombInfo:
    def __init__(self):
        #TO DO: Plan to remove
        self.comb = None
        self.query_dict = {}

        #scores
        self.totals = None
        self.scores = None
        self.fracScore = -10.00
        self.multiScore = -10.00  #Calc multiScore (By Bill: -ln(Na/SumNa * Nb/SumNb * Nc/SumNc))

        #Geometry
        self.geometry = None
        self.aa_aa_pair = None
        self.metal_aa_pair = None
        self.angle_pair = None

        
        self.volume = 0
        self.volPerMetal = 0
        self.diameter = 0

        #Querys
        self.centroid_dict = {}

        #evaluation property for evaluation search
        self.eval_mins = None
        self.eval_min_vdMs = None
        self.eval_is_origin = False

        #After search filter property
        self.pair_aa_aa_dist_ok = 0 #0: unchecked. -1: condition unsatisfied; 1: condition satisfied.
        self.pair_angle_ok = 0
        self.pair_metal_aa_dist_ok = 0
        self.vdm_no_clash = 0
        
        #For Search_selfcenter
        self.overlap_query_id_dict = None
        self.overlap_id_dict = None

        #After search filter result. If ture means any filter works.
        self.after_search_filtered = False


    def calc_geometry(self):
        all_coords = []
        metal_coords = []  

        for key in self.query_dict.keys():
            coords = []
            for _query in self.query_dict[key]:                                  
                coords.append(_query.get_contact_coord())
                metal_coords.append(_query.get_metal_coord())
            all_coords.append(pr.calcCenter(hull.transfer2pdb(coords)))         
        all_coords.append(pr.calcCenter(hull.transfer2pdb(metal_coords)))

        self.geometry = hull.transfer2pdb(all_coords, ['NI' if i == len(all_coords)-1 else 'N' for i in range(len(all_coords))])
        self.aa_aa_pair, self.metal_aa_pair, self.angle_pair  = quco.pair_wise_geometry(self.geometry)
        return
    
    def calc_centroid_geometry(self):
        all_coords = []
        metal_coords = []  

        for _query in self.centroid_dict.values():                                
            all_coords.append(_query.contact_ag)
            metal_coords.append(_query.candidates_metal_points)   
        all_coords.append(pr.calcCenter(hull.transfer2pdb(metal_coords)))

        self.geometry = hull.transfer2pdb(all_coords, ['NI' if i == len(all_coords)-1 else 'N' for i in range(len(all_coords))])
        self.aa_aa_pair, self.metal_aa_pair, self.angle_pair  = quco.pair_wise_geometry(self.geometry)
        return        

    def after_search_condition_satisfied(self, pair_angle_range = None, pair_aa_aa_dist_range = None, pair_metal_aa_dist_range = None):
        '''
        range = (75, 125) for Zn.
        if all pairwise angle is between the range. The geometry is satisfied.
        '''
        if pair_angle_range:
            for an in self.angle_pair:
                if an < pair_angle_range[0] or an > pair_angle_range[1]:
                    self.pair_angle_ok = -1
                    return False
                else:
                    self.pair_angle_ok = 1
        if pair_aa_aa_dist_range:           
            for ad in self.aa_aa_pair:
                if ad < pair_aa_aa_dist_range[0] or ad > pair_aa_aa_dist_range[1]:
                    self.pair_aa_aa_dist_ok = -1
                    return False
                else:
                    self.pair_aa_aa_dist_ok = 1

        if pair_metal_aa_dist_range:
            for amd in self.metal_aa_pair:
                if amd < pair_metal_aa_dist_range[0] or amd > pair_metal_aa_dist_range[1]:
                    self.pair_aa_metal_dist_ok = -1
                    return False
                else:
                    self.pair_aa_metal_dist_ok = 1

        return True