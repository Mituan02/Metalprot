import prody as pr
import numpy as np
from ..basic import vdmer
from ..basic import hull
from ..basic import constant


def supperimpose_ideal_geo(geometry):
    ideal_geometry = constant.tetrahydra_geo.copy()
    all_coords = geometry.getCoords()
    pr.calcTransformation(ideal_geometry.getCoords(), all_coords).apply(ideal_geometry)
    rmsd = pr.calcRMSD(ideal_geometry.getCoords(), all_coords)

    all_coords2 = np.array([all_coords[i] for i in [1, 0, 2, 3]])
    ideal_geometry2 = constant.tetrahydra_geo.copy()
    pr.calcTransformation(ideal_geometry2.getCoords(), all_coords2).apply(ideal_geometry2)
    rmsd2 = pr.calcRMSD(ideal_geometry2.getCoords(), all_coords2)

    if rmsd2 < rmsd:
        return ideal_geometry2, rmsd2
    return ideal_geometry, rmsd


class CombInfo:
    def __init__(self):
        #TO DO: Plan to remove
        self.comb = None
        self.query_dict = {}

        #scores
        self.overlaps = None
        self.scores = None
        self.cluScore = -10.00
        self.overlapScore = -10.00 

        #Geometry
        self.geometry = None
        self.geo_rmsd = -10.00
        self.aa_aa_pair = None
        self.metal_aa_pair = None
        self.angle_pair = None

        
        self.volume = 0
        self.volPerMetal = 0
        self.diameter = 0

        #Centroid
        self.centroid_dict = {}

        #2nd shell vdm
        self.secondshell_dict = {}

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
        self.overlap_ind_dict = None

        #After search filter result. If ture means the CombInfo pass all filter conditions.
        self.after_search_filtered = False

        #Tag
        self.tag = '' # tag if the CombInfo has the best OverlapScore or best rmsd or best ClusterScore


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
        self.aa_aa_pair, self.metal_aa_pair, self.angle_pair  = vdmer.pair_wise_geometry(self.geometry)
        return

    
    def calc_centroid_geometry(self):
        all_coords = []
        metal_coords = []  

        for _query in self.centroid_dict.values():                                
            all_coords.append(_query.get_contact_coord())
            metal_coords.append(_query.get_metal_coord())   
        all_coords.append(pr.calcCenter(hull.transfer2pdb(metal_coords)))

        self.geometry = hull.transfer2pdb(all_coords, ['NI' if i == len(all_coords)-1 else 'N' for i in range(len(all_coords))])      
        self.ideal_geometry = self.geometry
        
        self.aa_aa_pair, self.metal_aa_pair, self.angle_pair  = vdmer.pair_wise_geometry(self.geometry)
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