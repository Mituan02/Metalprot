'''
Prepare the lig-helix for ligand search and rosetta design.
'''
import prody as pr
import numpy as np


#>>> mvdms.
#>>> HIS/ASP sc can be transformed. 
#>>> GLU can not be transformed if including CB, one way is not including CB. Another way is to find a similar one and only change the coord of sc heavy atoms.


workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_14-135-139_HHE/'

mvdm = pr.parsePDB(workdir + 'structure_prep/W_A-14-A-135-A-139_H-H-E_351-334-451_BestOPscore_BestCluscore_Best2ShScore_Best2ShRmsd.pdb')

gvdm = pr.parsePDB(workdir + 'structure_prep/W_A-14-A-135-A-139_H-H-E_351-334-451_A14_hie_THR_A_136_0.53_2.44.pdb')


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

target = pr.parsePDB(workdir + 'helix6_preloop_mvdm.pdb')

his_order = ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']
coord_a = np.array([target.select('resnum 135 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid B and resnum 119 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 135 and sc').copy()
tf.apply(ag)
target.select('resnum 135 and sc').setCoords(ag.getCoords())
target.select('resnum 135 and sc').getNames()
target.select('resnum 135 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD1','HD2', 'HE1'])

coord_a = np.array([target.select('resnum 14 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid A and resnum 77 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 14 and sc').copy()
tf.apply(ag)
target.select('resnum 14 and sc').setCoords(ag.getCoords())
target.select('resnum 14 and sc').getNames()
target.select('resnum 14 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD2', 'HE1', 'HE2'])

target.select('resnum 139 and sc and heavy').setCoords(mvdm.select('chid C and resnum 77 and heavy and sc').getCoords())
target.select('resnum 139 and sc').setNames(['CB', 'CG', 'CD', 'OE1', 'OE2', 'HA', 'HB2', 'HB3', 'HG2', 'HG3'])
#>>> gvdms.

target.select('resnum 136 and sc').getNames()
#['CB', 'CG2', 'OG1', 'HA', 'HB', '1HG', '1HG2', '2HG2', '3HG2']
coord_b = np.array([gvdm.select('chid X and resnum 10 and name ' + x).getCoords()[0] for x in ['CB',  'CG2', 'OG1', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23']])
target.select('resnum 136 and sc').setCoords(coord_b)
target.select('resnum 136 and sc').setNames(['CB',  'CG2', 'OG1', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23'])

pr.writePDB(workdir + target.getTitle() + '_fix.pdb', target)



