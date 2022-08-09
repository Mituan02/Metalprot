'''
Prepare the lig-helix for ligand search and rosetta design.
'''

#>>> mvdms. Due to the lack of hydrogen in current mvdm library, there is a bug here for Hydrogen position.
#>>> HIS/ASP sc can be transformed. 
#>>> GLU can not be transformed if including CB, one way is not including CB. Another way is to find a similar one and only change the coord of sc heavy atoms.


import prody as pr
import numpy as np

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
workdir = '/mnt/e/DesignData/Metalloenzyme/SAHA_Vorinostat/run_design_cgs3/SAHA_Rosetta20220730_2/'

target = pr.parsePDB(workdir + 'bb_loop.pdb')

mvdm = pr.parsePDB(workdir + 'W_A-54-A-128-A-132_H-E-H_846-611-370_BestOPscore_BestGeo_BestCluscore_Best2ShScore_Best2ShRmsd.pdb')

gvdm_A17 = pr.parsePDB(workdir + 'W_A-54-A-128-A-132_H-E-H_846-611-370_A54_hid_SER_A_17_0.62_1.33.pdb')
gvdm_A125 = pr.parsePDB(workdir + 'W_A-54-A-128-A-132_H-E-H_846-611-370_A128_coo_LYS_A_125_0.27_1.52.pdb')
gvdm_A135 = pr.parsePDB(workdir + 'W_A-54-A-128-A-132_H-E-H_846-611-370_A132_hid_GLU_A_135_0.62_3.02.pdb')


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
his_order = ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']

#>>> hid
target.select('resnum 130 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '1HD', '2HD', '1HE']
coord_a = np.array([target.select('resnum 130 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid A and resnum 63 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 130 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').copy()
tf.apply(ag)
target.select('resnum 130 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').setCoords(ag.getCoords())
target.select('resnum 130 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1'])

#>>> hid
target.select('resnum 61 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '1HD', '2HD', '1HE']
coord_a = np.array([target.select('resnum 61 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid C and resnum 301 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 61 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').copy()
tf.apply(ag)
target.select('resnum 61 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').setCoords(ag.getCoords())
target.select('resnum 61 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1'])

#>>> GLU
target.select('resnum 57 and sc').getNames()
mvdm.select('chid B and resnum 37 and sc').getNames() 
target.select('resnum 57 and sc').setCoords(mvdm.select('chid B and resnum 37 and sc').getCoords())
target.select('resnum 57 and sc').setNames(['CB', 'CG', 'CD', 'OE1', 'OE2', 'HA', 'HB2', 'HB3', 'HG2', 'HG3'])

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> mvdm2ndshell.

#>>> GLU
target.select('resnum 64 and sc').getNames()
gvdm_A135.select('chid X and resnum 10 and sc').getNames() 
target.select('resnum 64 and sc').setCoords(gvdm_A135.select('chid X and resnum 10 and sc').getCoords())
target.select('resnum 64 and sc').setNames(gvdm_A135.select('chid X and resnum 10 and sc').getNames() )

#>>> LYS
target.select('resnum 54 and sc').getNames()
gvdm_A125.select('chid X and resnum 10 and sc').getNames() 
target.select('resnum 54 and sc').setCoords(gvdm_A125.select('chid X and resnum 10 and sc').getCoords())
target.select('resnum 54 and sc').setNames(gvdm_A125.select('chid X and resnum 10 and sc').getNames())


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> gvdms

pr.writePDB(workdir + target.getTitle() + '_fix.pdb', target)