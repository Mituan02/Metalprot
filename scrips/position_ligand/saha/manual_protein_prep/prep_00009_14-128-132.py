'''
Prepare the lig-helix for ligand search and rosetta design.
'''
import prody as pr
import numpy as np


#>>> mvdms. Due to the lack of hydrogen in current mvdm library, there is a bug here for Hydrogen position.
#>>> HIS/ASP sc can be transformed. 
#>>> GLU can not be transformed if including CB, one way is not including CB. Another way is to find a similar one and only change the coord of sc heavy atoms.


workdir = '/mnt/e/DesignData/Metalloprotein/Design_00009_54-128-132_H-H-H/'

mvdm = pr.parsePDB(workdir + 'W_A-54-A-128-A-132_H-H-H_2633-300-372_BestOPscore_BestCluscore_Best2ShScore_Best2ShRmsd.pdb')

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

target = pr.parsePDB(workdir + 'bb_prep.pdb')

#gvdm_A62 = pr.parsePDB(workdir + 'Lig-318_phenol_A62_ASP_rmsd_0.74_v_2.8_1293.pdb.gz')
#gvdm_A65 = pr.parsePDB(workdir + 'Lig-318_ph_A65_PHE_rmsd_0.34_v_2.3_19174.pdb.gz') 
#gvdm_A84 = pr.parsePDB(workdir + 'Lig-318_bb_cco_A84_SER_rmsd_0.61_v_2.6_4224.pdb.gz')

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
his_order = ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']

#>>> hid
target.select('resnum 57 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '1HD', '2HD', '1HE']
coord_a = np.array([target.select('resnum 57 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid B and resnum 475 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 57 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').copy()
tf.apply(ag)
target.select('resnum 57 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').setCoords(ag.getCoords())
target.select('resnum 57 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1'])

#>>> hid
target.select('resnum 61 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '1HD', '2HD', '1HE']
coord_a = np.array([target.select('resnum 61 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid C and resnum 500 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 61 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').copy()
tf.apply(ag)
target.select('resnum 61 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').setCoords(ag.getCoords())
target.select('resnum 61 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1'])


#>>> hid
target.select('resnum 130 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '1HD', '2HD', '1HE']
coord_a = np.array([target.select('resnum 130 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid A and resnum 36 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 130 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').copy()
tf.apply(ag)
target.select('resnum 130 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').setCoords(ag.getCoords())
target.select('resnum 130 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1'])



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> mvdm2ndshell.


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> gvdms


pr.writePDB(workdir + target.getTitle() + '_fix.pdb', target)