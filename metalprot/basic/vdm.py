import prody as pr
import itertools
import numpy as np
from numpy.core.fromnumeric import argmin
from prody.atomic import pointer
from .hull import transfer2pdb, write2pymol
from .utils import get_ABPLE

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

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


def pair_wise_geometry(geometry_ag):
    '''
    cts: contact atoms
    '''
    metal = geometry_ag.select('name NI')[0]
    cts = geometry_ag.select('name N')
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


class VDM:
    def __init__(self, query, id = -1, clu_rank = -1, score = 0, clu_num = 0, clu_total_num = 0, clu_member_ids= None, metal_atomgroup = None, win = None, path = None):
        '''
        Note that the query or the metal_atomgroup are prody object, which may be transformed during the searching.
        '''
        self.query = query # The prody vdM pdb.
        
        self.id = id # Each vdM has a unique id, for index the vdM library. 
        self.clu_rank = clu_rank # The rank of cluster in the same type of vdM.

        # vdM info
        self.score = score
        self.clu_num = clu_num
        self.clu_total_num = clu_total_num

        # cluster member info
        self.clu_member_ids = clu_member_ids
        self.metal_atomgroup = metal_atomgroup
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
        self.aa_type = self.query.select('resindex ' + str(self.contact_resind)).getResnames()[0] # what the contacting amino acid.
        
        self.get_phi_psi()
        self.abple = get_ABPLE(self.query.select('name CA and resindex ' + str(self.contact_resind)).getResnames()[0], self.phi, self.psi)


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
        return VDM(self.query.copy(), self.id, self.clu_rank, self.score, self.clu_num, self.clu_total_num, self.clu_member_ids, metal_atomgroup, self.win, self.path)
    
    def writepdb(self, outpath):
        pr.writePDB(outpath, self.query)