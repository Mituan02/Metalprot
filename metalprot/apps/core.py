from logging import NullHandler
import os
import prody as pr
from scipy.spatial.distance import cdist
from dataclasses import dataclass
import shutil
import sys
import numpy as np
import matplotlib.pyplot as plt
from . import cluster
from . import transformation
from . import ligand_database

# According to the prody atom flag (http://prody.csb.pitt.edu/manual/reference/atomic/flags.html#flags), NI MN ZN CO CU MG FE CA are not all flag as ion. 
# Note that calcium is read as CA, which is same as alpha carbon in prody selection. 
# So to select calcium, we need to add ion before it.

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

class Core:
    def __init__(self, full_pdb):
        self.full_pdb = full_pdb

        self.metal = full_pdb.select(metal_sel)[0]
        self.metal_resind = self.metal.getResindex()

        self.contact_aas = self.full_pdb.select('protein and within 2.83 of resindex ' + str(self.metal.getResindex()))
        self.contact_aa_resinds = np.unique(self.contact_aas.getResindices())    

        # Inside the atomGroupDict. We will generate different type of pdb for clustering into vdM.
        self.atomGroupDict = dict()

    def add2atomGroupDict(self, key, sel_pdb_prody):
        if not key in self.atomGroupDict.keys():
            self.atomGroupDict[key] = []
        self.atomGroupDict[key].append(sel_pdb_prody)

        
    def generate_AA_Metal(self, AA = 'HIS', key = 'AAMetal_HIS'):
        aa_sel = 'resname ' + AA

        for resind in self.contact_aa_resinds:
            if not self.full_pdb.select(aa_sel + ' and resindex ' + str(resind)):
                continue
            sel_pdb_prody = self.full_pdb.select('resindex ' + str(resind) + ' '+ str(self.metal_resind))
            
            self.add2atomGroupDict(key, sel_pdb_prody)            
        return 
        
    
    def generate_AA_ext_Metal(self, AA = 'HIS', extention = 3, key = 'AAExtMetal_HIS_ext3'):
        aa_sel = 'resname ' + AA

        for resind in self.contact_aa_resinds:
            if not self.full_pdb.select(aa_sel + ' and resindex ' + str(resind)):
                continue

            ext_inds = ligand_database.extend_res_indices([resind], self.full_pdb, extend =extention)
            if len(ext_inds) != 2*extention + 1:
                continue
            #print(self.full_pdb.getTitle() + '+' + '-'.join([str(x) for x in ext_inds]))
            sel_pdb_prody = self.full_pdb.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(self.metal_resind))

            self.add2atomGroupDict(key, sel_pdb_prody)
        return

    
    def generate_AA_phipsi_Metal(self, AA = 'HIS', key = 'AAMetalPhiPsi_HIS'):
        aa_sel = 'resname ' + AA

        for resind in self.contact_aa_resinds:
            if not self.full_pdb.select(aa_sel + ' and resindex ' + str(resind)):
                continue
            
            ext_inds = ligand_database.extend_res_indices([resind], self.full_pdb, extend =1)
            if len(ext_inds) != 3:
                continue
            #print(self.full_pdb.getTitle() + '+' + '-'.join([str(x) for x in ext_inds]))
            atom_inds = []
            atom_inds.extend(self.full_pdb.select('resindex ' + str(resind-1)).select('name C O').getIndices())
            atom_inds.extend(self.full_pdb.select('resindex ' + str(resind)).getIndices())
            atom_inds.extend(self.full_pdb.select('resindex ' + str(resind+1)).select('name N CA').getIndices())
            sel_pdb_prody = self.full_pdb.select('index ' + ' '.join([str(x) for x in atom_inds]) + ' '+ str(self.metal.getIndex()))
            
            self.add2atomGroupDict(key, sel_pdb_prody)
        return


    def generate_AAcAA_Metal(self, filter_aas = False, aas = ['HIS', 'HIS'], extention = 3, extention_out = 0, key = 'AAcAA_Metal_ext0'):
        exts = []
        for resind in self.contact_aa_resinds:
            ext_inds = ligand_database.extend_res_indices([resind], self.full_pdb, extend =extention)
            exts.append(ext_inds)
        
        pairs = []
        pairs_ind = []
        for i in range(len(self.contact_aa_resinds) - 1):
            for j in range(i+1, len(self.contact_aa_resinds)):
                overlap = list(set(exts[i]) & set(exts[j]))
                if len(overlap) ==0: continue

                if filter_aas:
                    filtered = True
                    if aas and len(aas)==2:
                        if self.full_pdb.select('resindex ' + str(self.contact_aa_resinds[i])).getResnames()[0] == aas[0] and self.full_pdb.select('resindex ' + str(self.contact_aa_resinds[j])).getResnames()[0] == aas[1]:
                            filtered= False
                        elif self.full_pdb.select('resindex ' + str(self.contact_aa_resinds[i])).getResnames()[0] == aas[1] and self.full_pdb.select('resindex ' + str(self.contact_aa_resinds[j])).getResnames()[0] == aas[0]:
                            filtered= False
                    if filtered: continue

                pairs.append((i, j))
                if self.contact_aa_resinds[i] < self.contact_aa_resinds[j]:
                    pairs_ind.append((self.contact_aa_resinds[i], self.contact_aa_resinds[j]))
                else:
                    pairs_ind.append((self.contact_aa_resinds[j], self.contact_aa_resinds[i]))
        
        for v in range(len(pairs)):
            i = pairs[v][0]
            j = pairs[v][1]
            ext_inds = list(set(exts[i]) | set(exts[j]))
            ext_inds = [x for x in ext_inds if x >= pairs_ind[v][0]-extention_out and x <= pairs_ind[v][1]+extention_out]
            sel_pdb_prody = self.full_pdb.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(self.metal_resind))
            
            self.add2atomGroupDict(key, sel_pdb_prody)
        return


    def generate_AAdAA_Metal(self, filter_aas = False, aas = ['HIS', 'HIS'], key = 'AAdAA_Metal'):
        pairs_ind = []
        for i in range(len(self.contact_aa_resinds) - 1):
            for j in range(i+1, len(self.contact_aa_resinds)):

                if filter_aas:
                    filtered = True
                    if aas and len(aas)==2:
                        if self.full_pdb.select('resindex ' + str(self.contact_aa_resinds[i])).getResnames()[0] == aas[0] and self.full_pdb.select('resindex ' + str(self.contact_aa_resinds[j])).getResnames()[0] == aas[1]:
                            filtered= False
                        elif self.full_pdb.select('resindex ' + str(self.contact_aa_resinds[i])).getResnames()[0] == aas[1] and self.full_pdb.select('resindex ' + str(self.contact_aa_resinds[j])).getResnames()[0] == aas[0]:
                            filtered= False
                    if filtered: continue

                pairs_ind.append((self.contact_aa_resinds[i], self.contact_aa_resinds[j]))
                pairs_ind.append((self.contact_aa_resinds[j], self.contact_aa_resinds[i]))  # In the database, we don't know the order of the two aa.
        
        for p in pairs_ind:
            sel_pdb_prody = self.full_pdb.select('resindex ' + str(p[0]) + ' ' +  str(p[1]) + ' '+ str(self.metal_resind))
            
            self.add2atomGroupDict(key, sel_pdb_prody)
        return


    def get_inds_from_resind(self, pdb_prody, resind, aa = 'resname HIS'):
        if aa == 'resname HIS':
            inds = pdb_prody.select('name ND1 NE2 and resindex ' + str(resind)).getIndices()
        elif aa == 'resname ASP':
            inds = pdb_prody.select('name OD1 OD2 and resindex ' + str(resind)).getIndices()
        elif aa == 'resname GLU':
            inds = pdb_prody.select('name OE1 OE2 and resindex ' + str(resind)).getIndices()
        else:
            inds = []
        return inds


    def generate_AA_2ndShell_Metal(self, key = 'AA2ndShellMetal'):
        for resind in self.contact_aa_resinds:      
            #inds = get_inds_from_resind(pdb_prody, resind, aa)
            _2nshell_resinds = self.get_2ndshell_indices([resind], self.full_pdb, self.metal_resind.getIndex())
            if len(_2nshell_resinds) > 0:
                for _2resind in _2nshell_resinds:
                    #print(self.full_pdb.getTitle() + '+' + '-'.join([str(x) for x in _2nshell_resinds]))
                    sel_pdb_prody = self.full_pdb.select('resindex ' + str(resind) + ' ' + str(_2resind) + ' ' + str(self.metal_resind))
                  
                    self.add2atomGroupDict(key, sel_pdb_prody)
        return


    def generate_AA_kMetal(self, AA = 'HIS', k = 5, key = 'AAMetal_HIS', key_out = 'AA5Metal_HIS'):        
        if not key in self.atomGroupDict:
            return
        for ag in self.atomGroupDict[key]:
            sel_pdb_prody = ag.select('protein').toAtomGroup() 
            for i in range(k):
                sel_pdb_prody = sel_pdb_prody + self.metal.toAtomGroup()

            if not key_out in self.atomGroupDict.keys():
                self.atomGroupDict[key_out] = []
            self.atomGroupDict[key_out].append(sel_pdb_prody)
        return


    def write_vdM(self, outdir, key):
        if not key in self.atomGroupDict.keys():
            return

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        count = 0
        for ag in self.atomGroupDict[key]:
            pr.writePDB(outdir + self.full_pdb.getTitle() + '_' + key + '_mem' + str(count) + '.pdb', ag)
            count+=1
            

