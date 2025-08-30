'''
Prepare the lig-helix for ligand search and rosetta design.
'''
import prody as pr
import numpy as np


#>>> mvdms. Due to the lack of hydrogen in current mvdm library, there is a bug here for Hydrogen position.
#>>> HIS/ASP sc can be transformed. 
#>>> GLU can not be transformed if including CB, one way is not including CB. Another way is to find a similar one and only change the coord of sc heavy atoms.


workdir = '/mnt/e/DesignData/Metalloenzyme/1ukr/68-77-79_DHH/'

target = pr.parsePDB(workdir + 'a_1ukr_68-77-79_DHH_lig20.pdb')

mvdm = pr.parsePDB(workdir + 'W_A-68-A-77-A-79_D-H-H_294-2334-1658_Best2ShScore.pdb')

mvdm2nd_A70 = pr.parsePDB(workdir + 'W_A-68-A-77-A-79_D-H-H_294-2334-1658_A77_hid_ASP_A_70_0.25_2.8.pdb')

mvdm2nd_A64 = pr.parsePDB(workdir + 'W_A-68-A-77-A-79_D-H-H_294-2334-1658_A79_hie_GLU_A_64_0.58_0.87.pdb')

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

his_order = ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']
#>>> HIE

# target.select('resnum 77 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '2HD', '1HE', '2HE']
# coord_a = np.array([target.select('resnum 60 and name ' + x).getCoords()[0] for x in his_order])
# coord_b = np.array([mvdm.select('chid C and resnum 77 and name ' + x).getCoords()[0] for x in his_order])
# tf  = pr.calcTransformation(coord_a, coord_b)
# ag = target.select('resnum 77 and name CB CG CD2 ND1 CE1 NE2 2HD 1HE 2HE').copy()
# tf.apply(ag)
# target.select('resnum 77 and name CB CG CD2 ND1 CE1 NE2 2HD 1HE 2HE').setCoords(ag.getCoords())
# target.select('resnum 77 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD2', 'HE1', 'HE2'])

#>>> HID
target.select('resnum 77 and sc').getNames() 
coord_a = np.array([target.select('resnum 77 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid B and resnum 82 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 77 and name CB CG CD2 ND1 CE1 NE2 HD1 HD2 HE1').copy()
tf.apply(ag)
target.select('resnum 77 and name CB CG CD2 ND1 CE1 NE2 HD1 HD2 HE1').setCoords(ag.getCoords())

#>>> HID
target.select('resnum 79 and sc').getNames() 
coord_a = np.array([target.select('resnum 79 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid C and resnum 377 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 79 and name CB CG CD2 ND1 CE1 NE2 HD1 HD2 HE1').copy()
tf.apply(ag)
target.select('resnum 79 and name CB CG CD2 ND1 CE1 NE2 HD1 HD2 HE1').setCoords(ag.getCoords())


#>>> ASP
target.select('resnum 68 and sc').getNames() #['CB', 'CG', 'OD1', 'OD2', 'HA', '2HB', '3HB']
mvdm.select('chid A and resnum 36 and heavy and sc').getNames()
target.select('resnum 68 and sc and heavy').setCoords(mvdm.select('chid A and resnum 36 and heavy and sc').getCoords())
target.select('resnum 68 and sc').setNames(['CB', 'CG', 'OD1', 'OD2', 'HA', 'HB1', 'HB2'])

#>>> mvdm2ndshell.
target.select('resnum 70 and sc').getNames()
#['CB', 'CG', 'OD1', 'OD2', 'HA', 'HB2', 'HB3']
mvdm2nd_A70.select('chid X and resnum 10').getNames()
coord_b = np.array([mvdm2nd_A70.select('chid X and resnum 10 and name ' + x).getCoords()[0] for x in ['CB', 'CG', 'OD1', 'OD2', 'HA', 'HB2', 'HB3']])
target.select('resnum 70 and sc').setCoords(coord_b)
target.select('resnum 70 and sc').setNames(['CB', 'CG', 'OD1', 'OD2', 'HA', 'HB2', 'HB3'])


target.select('resnum 64 and sc').getNames()
#['CB', 'CG', 'CD', 'OE1', 'OE2', 'HA', 'HB2', 'HB3', 'HG2', 'HG3']
coord_b = np.array([mvdm2nd_A64.select('chid X and resnum 10 and name ' + x).getCoords()[0] for x in ['CB', 'CG', 'CD', 'OE1', 'OE2', 'HA', 'HB2', 'HB3', 'HG2', 'HG3']])
target.select('resnum 64 and sc').setCoords(coord_b)
target.select('resnum 64 and sc').setNames(['CB', 'CG', 'CD', 'OE1', 'OE2', 'HA', 'HB2', 'HB3', 'HG2', 'HG3'])
#>>> gvdms


pr.writePDB(workdir + target.getTitle() + '_fix.pdb', target)