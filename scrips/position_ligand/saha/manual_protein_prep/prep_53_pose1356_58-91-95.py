'''
Prepare the lig-helix for ligand search and rosetta design.
'''
import prody as pr
import numpy as np


#>>> mvdms. Due to the lack of hydrogen in current mvdm library, there is a bug here for Hydrogen position.
#>>> HIS/ASP sc can be transformed. 
#>>> GLU can not be transformed if including CB, one way is not including CB. Another way is to find a similar one and only change the coord of sc heavy atoms.


workdir = '/mnt/e/DesignData/Metalloenzyme/SAHA_Vorinostat/run_design_cgs3/SAHA_Rosetta20220728/'

mvdm = pr.parsePDB(workdir + 'W_A-58-A-91-A-95_H-H-E_2079-2414-463_BestOPscore_BestCluscore_Best2ShScore.pdb')

gvdm_A62 = pr.parsePDB(workdir + 'W_A-58-A-91-A-95_H-H-E_2079-2414-463_A58_hie_TYR_A_62_0.48_0.64.pdb')
gvdm_A99 = pr.parsePDB(workdir + 'W_A-58-A-91-A-95_H-H-E_2079-2414-463_A95_coo_LYS_A_99_0.3_2.62.pdb')


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

target = pr.parsePDB(workdir + '00009.f63440efff7e.allbb_ala_min_ala_0001.pdb')

#gvdm_A62 = pr.parsePDB(workdir + 'Lig-318_phenol_A62_ASP_rmsd_0.74_v_2.8_1293.pdb.gz')
#gvdm_A65 = pr.parsePDB(workdir + 'Lig-318_ph_A65_PHE_rmsd_0.34_v_2.3_19174.pdb.gz') 
#gvdm_A84 = pr.parsePDB(workdir + 'Lig-318_bb_cco_A84_SER_rmsd_0.61_v_2.6_4224.pdb.gz')

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
his_order = ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']

#>>> hid
target.select('resnum 58 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '1HD', '2HD', '1HE']
coord_a = np.array([target.select('resnum 58 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid A and resnum 44 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 58 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').copy()
tf.apply(ag)
target.select('resnum 58 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').setCoords(ag.getCoords())
target.select('resnum 58 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1'])

#>>> hid
target.select('resnum 91 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '1HD', '2HD', '1HE']
coord_a = np.array([target.select('resnum 91 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid B and resnum 243 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 91 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').copy()
tf.apply(ag)
target.select('resnum 91 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').setCoords(ag.getCoords())
target.select('resnum 91 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1'])

#>>> GLU
target.select('resnum 95 and sc').getNames()
mvdm.select('chid C and resnum 152 and heavy and sc').getNames() #['CB', 'CG', 'CD', 'OE1', 'OE2']
target.select('resnum 95 and sc and heavy').setCoords(mvdm.select('chid C and resnum 152 and heavy and sc').getCoords())
target.select('resnum 95 and sc').setNames(['CB', 'CG', 'OD1', 'OD2', 'HA', 'HB1', 'HB2'])


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> mvdm2ndshell.

target.select('resnum 91 and sc').getNames()
#['CB', 'CG2', 'OG1', 'HA', 'HB', '1HG', '1HG2', '2HG2', '3HG2']
coord_b = np.array([gvdm.select('chid X and resnum 10 and name ' + x).getCoords()[0] for x in ['CB',  'CG2', 'OG1', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23']])
target.select('resnum 91 and sc').setCoords(coord_b)
target.select('resnum 91 and sc').setNames(['CB',  'CG2', 'OG1', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23'])

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> gvdms

target.select('resnum 131 and sc').getNames()
#['CB', 'CG', 'OD1', 'OD2', 'HA', '2HB', '3HB']
gvdm_A62.select('chid X and sc and resnum 10').getNames()
coord_b = gvdm_A62.select('chid X and sc and resnum 10').getCoords()
target.select('resnum 131 and sc').setCoords(coord_b)
target.select('resnum 131 and sc').setNames(gvdm_A62.select('chid X and sc and resnum 10').getNames())


target.select('resnum 134 and sc').getNames()
#['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HA', '2HB', '3HB', '1HD', '2HD', '1HE', '2HE', 'HZ']
gvdm_A65.select('chid X and sc and resnum 10').getNames()
coord_b = gvdm_A65.select('chid X and sc and resnum 10').getCoords()
target.select('resnum 134 and sc').setCoords(coord_b)
target.select('resnum 134 and sc').setNames(gvdm_A65.select('chid X and sc and resnum 10').getNames())

target.select('resnum 16 and sc').getNames()
#['CB', 'OG', 'HA', '2HB', '3HB', 'HG']
gvdm_A84.select('chid X and sc and resnum 10').getNames()
coord_b = gvdm_A84.select('chid X and sc and resnum 10').getCoords()
target.select('resnum 16 and sc').setCoords(coord_b)
target.select('resnum 16 and sc').setNames(gvdm_A84.select('chid X and sc and resnum 10').getNames())

pr.writePDB(workdir + target.getTitle() + '_fix.pdb', target)