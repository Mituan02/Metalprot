import prody as pr
import itertools
import numpy as np
from numpy.core.fromnumeric import argmin
from .utils import get_ABPLE
from . import constant 

metal_sel = 'name NI MN ZN CO CU MG FE' 

def get_metal_contact_atoms(pdb):
    '''
    Extract list of contact atoms from pdb. 
    Here the pdb could be the core that contain multiple contact atoms.
    '''
    #Get metal, binding atom for each binding atom
    cts = []

    metal = pdb.select(metal_sel)[0]
    cts.append(metal)
    _contact_aas = pdb.select('protein and not carbon and not hydrogen and within 2.83 of resindex ' + str(metal.getResindex()))
    #For each aa, only select one contact atom. 
    resindices = np.unique(_contact_aas.getResindices())
    for rid in resindices:
        #TO DO: Not the first one, but the closest one for the  such ASP contains two contact atoms.
        _ct = _contact_aas.select('resindex ' + str(rid))
        if len(_ct) > 1:
            dists = [None]*len(_ct)
            for i in range(len(_ct)):
                dist = pr.calcDistance(metal, _ct[i])
                dists[i] = dist
            ct = _ct[argmin(dists)]
        else:
            ct = _ct[0]
        cts.append(ct)

    return cts

def get_contact_atom(pdb):
    '''
    Get contact atom of a vdM.
    '''
    metal = pdb.select(metal_sel)[0]
    _contact_aas = pdb.select('protein and not carbon and not hydrogen and within 2.83 of resindex ' + str(metal.getResindex()))
    if not _contact_aas:
        print('No contact atom: ' + pdb.getTitle())
        _contact_aas = pdb.select('protein and not carbon and not hydrogen and within 5 of resindex ' + str(metal.getResindex()))      
    if len(_contact_aas) > 1:
        dists = [None]*len(_contact_aas)
        for i in range(len(_contact_aas)):
            dist = pr.calcDistance(metal, _contact_aas[i])
            dists[i] = dist
        contact_aa = _contact_aas[argmin(dists)]
    else:
        contact_aa = _contact_aas[0]
    return contact_aa


def pair_wise_geometry(geometry_ag, metal_sel = 'name NI MN ZN CO CU MG FE'):
    '''
    cts: contact atoms
    '''
    metal = geometry_ag.select(metal_sel)[0]
    cts = geometry_ag.select('not ' + metal_sel)
    ct_len = len(cts)
    #print(cts.getNames())
    aa_aa_pair = []
    metal_aa_pair =[]
    angle_pair = []
    for i, j in itertools.combinations(range(ct_len), 2):   
        dist = pr.calcDistance(cts[i], cts[j])
        aa_aa_pair.append(dist)
        angle = pr.calcAngle(cts[i], metal, cts[j])
        angle_pair.append(angle)
    for i in range(ct_len):
        metal_aa_pair.append(pr.calcDistance(cts[i], metal))
        
    return aa_aa_pair, metal_aa_pair, angle_pair  


def pair_wise_geometry_matrix(geometry_ag, metal_sel = 'name NI MN ZN CO CU MG FE'):
    '''
    cts: contact atoms
    '''
    metal = geometry_ag.select(metal_sel)[0]
    cts = geometry_ag.select('not ' + metal_sel)
    ct_len = len(cts)
    #print(cts.getNames())
    aa_aa_pair = np.zeros((ct_len, ct_len), dtype=float)
    angle_pair = np.zeros((ct_len, ct_len), dtype=float)

    metal_aa_pair = np.zeros(ct_len, dtype=float)

    for i, j in itertools.combinations(range(ct_len), 2):   
        dist = pr.calcDistance(cts[i], cts[j])
        aa_aa_pair[i, j] = dist
        angle = pr.calcAngle(cts[i], metal, cts[j])
        angle_pair[i, j] = angle
    for i in range(ct_len):
        metal_aa_pair[i] = pr.calcDistance(cts[i], metal)
        
    return aa_aa_pair, metal_aa_pair, angle_pair  


def adjust_metal_bond_lenght(pdb_prody, bond_length_diff):
    '''
    To combine vdM with different metal, the metal bond distance is adjusted based on the database.
    It is tricky to adjust the distance between Metal-Oxygen(ASP/GLU) if considering the bidentate contact.
    '''
    return


class VDM:
    def __init__(self, query, id = -1, clu_rank = -1, score = 0, clu_num = 0, max_clu_num = 0, clu_total_num = 0, clu_member_ids= None, metal_atomgroup = None, metalcontact_atomgroup = None, sc_atomgroup = None, sc_atomgroup_ids = None, win = None, path = None):
        '''
        Note that the query or the metal_atomgroup are prody object, which may be transformed during the searching.
        '''
        self.query = query # The prody vdM pdb.
        
        self.id = id # Each vdM has a unique id, for index the vdM library. 
        self.clu_rank = clu_rank # The rank of cluster in the same type of vdM.

        # vdM info
        self.score = score
        self.clu_num = clu_num
        self.max_clu_num = max_clu_num
        self.clu_total_num = clu_total_num

        # cluster member info
        self.clu_member_ids = clu_member_ids
        self.metal_atomgroup = metal_atomgroup
        self.metalcontact_atomgroup = metalcontact_atomgroup
        self.sc_atomgroup = sc_atomgroup
        self.sc_atomgroup_ids = sc_atomgroup_ids
        self.candidate_inds = None 

        # search info
        self.win = win

        # where is the original vdm file.
        self.path = path

        #Calc contact_resind
        self.contact_resind = None
        self.aa_type = None
        self.phi = None
        self.psi = None
        self.set_vdm()


    def set_vdm(self):
        cen_inds = np.unique(self.query.getResindices())
        self.contact_resind = cen_inds[int(cen_inds.shape[0]/2)-1]
        aa_type = self.query.select('resindex ' + str(self.contact_resind)).getResnames()[0] # what the contacting amino acid.
        self.aa_type = constant.one_letter_code[aa_type]
        
        self.get_phi_psi()
        self.abple = get_ABPLE(self.query.select('name CA and resindex ' + str(self.contact_resind)).getResnames()[0], self.phi, self.psi)
        return


    def get_phi_psi(self):
        '''
        Get phi psi angle of the contact metal.
        '''
        indices = self.query.select('name N C CA').getIndices() 
        atoms = [self.query.select('index ' + str(i)) for i in indices]
        self.phi = pr.calcDihedral(atoms[0], atoms[1], atoms[2], atoms[3])[0]
        self.psi = pr.calcDihedral(atoms[1], atoms[2], atoms[3], atoms[4])[0]
        return 


    def get_cluster_key(self):
        return (self.aa_type, self.clu_rank)

    def get_metal_coord(self):
        return self.query.select(metal_sel)[0].getCoords()

    def get_metal_mem_coords(self):
        return self.metal_atomgroup.getCoords()
    
    def get_metalcontact_mem_coords(self):
        return self.metalcontact_atomgroup.getCoords()

    def get_sc_mem_coords_and_ids(self):
        '''
        The sc_atomgroup here is used for clashing filter. It doesn't contain CB.
        '''
        return self.sc_atomgroup.getCoords(), self.sc_atomgroup_ids

    def get_contact_coord(self):
        atm = get_contact_atom(self.query)           
        return atm.getCoords()

    def get_candidate_metal_coords(self):
        return self.metal_atomgroup.select('index ' + ' '.join([str(x) for x in self.candidate_inds])).getCoords()

    def get_win_str(self):
        return '-'.join([str(w) for w in self.win])

        
    def to_tab_string(self):
        query_info = self.query.getTitle() + '\t' + str(round(self.score, 2)) + '\t' + str(self.clu_num)  + '\t'+ str(self.clu_total_num)
        return query_info

    def copy(self):
        metal_atomgroup = None
        if self.metal_atomgroup:
            metal_atomgroup = self.metal_atomgroup.copy()
        metalcontact_atomgroup = None
        if self.metalcontact_atomgroup:
            metalcontact_atomgroup = self.metalcontact_atomgroup.copy()
        sc_atomgroup = None
        sc_atomgroup_ids = None
        if self.sc_atomgroup:
            sc_atomgroup = self.sc_atomgroup.copy()
            sc_atomgroup_ids = self.sc_atomgroup_ids.copy()
        return VDM(self.query.copy(), self.id, self.clu_rank, self.score, self.clu_num, self.max_clu_num, self.clu_total_num, self.clu_member_ids, metal_atomgroup, metalcontact_atomgroup, sc_atomgroup, sc_atomgroup_ids, self.win, self.path)
    
    def writepdb(self, outpath):
        pr.writePDB(outpath, self.query)




    