import prody as pr
import numpy as np
from ..basic import vdmer
from ..basic import prody_ext
from ..basic import constant


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
        self.ideal_geo = None
        self.geo_rmsd = 10.00
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
        self.secondshell_score = -10.00 
        self.secondshell_rmsd = 10.0
        #evaluation property for evaluation search
        self.eval_mins = None
        self.eval_min_vdMs = None
        self.eval_is_origin = False
        self.eval_contain_origin = False

        #After search filter property
        self.pair_aa_aa_dist_ok = 0 #0: unchecked. -1: condition unsatisfied; 1: condition satisfied.
        self.pair_angle_ok = 0
        self.pair_metal_aa_dist_ok = 0
        self.vdm_no_clash = 0
        self.open_site_clash = 0
        
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
            all_coords.append(pr.calcCenter(prody_ext.transfer2pdb(coords)))         
        all_coords.append(pr.calcCenter(prody_ext.transfer2pdb(metal_coords)))

        self.geometry = prody_ext.transfer2pdb(all_coords, ['NI' if i == len(all_coords)-1 else 'N' for i in range(len(all_coords))])
        self.aa_aa_pair, self.metal_aa_pair, self.angle_pair  = vdmer.pair_wise_geometry(self.geometry)
        return

    
    def calc_centroid_geometry(self):
        all_coords = []
        metal_coords = []  

        for _query in self.centroid_dict.values():                                
            all_coords.append(_query.get_contact_coord())
            metal_coords.append(_query.get_metal_coord())   
        all_coords.append(pr.calcCenter(prody_ext.transfer2pdb(metal_coords)))

        self.geometry = prody_ext.transfer2pdb(all_coords, ['NI' if i == len(all_coords)-1 else 'N' for i in range(len(all_coords))])      
        
        self.aa_aa_pair, self.metal_aa_pair, self.angle_pair  = vdmer.pair_wise_geometry(self.geometry)
        return     


    def calc_2ndshell(self):
        if len(list(self.secondshell_dict.keys())) == 0:
            return

        min_rmsd = 10.0
        for key in self.secondshell_dict:
            max_score = -10.0
            
            for ag, info2ds in self.secondshell_dict[key]:
                if info2ds[7] > max_score:
                    max_score = info2ds[7]
                if info2ds[5] < min_rmsd:
                    min_rmsd = info2ds[5]
            self.secondshell_score+= max_score
        self.secondshell_rmsd = min_rmsd
        return

