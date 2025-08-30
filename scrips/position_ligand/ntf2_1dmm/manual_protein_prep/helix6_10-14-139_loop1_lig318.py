'''
Prepare the lig-helix for ligand search and rosetta design.
'''
import prody as pr
import numpy as np


#>>> mvdms. Due to the lack of hydrogen in current mvdm library, there is a bug here for Hydrogen position.
#>>> HIS/ASP sc can be transformed. 
#>>> GLU can not be transformed if including CB, one way is not including CB. Another way is to find a similar one and only change the coord of sc heavy atoms.


workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_10-14-139_HDH/Lig_318/loop_0/'

mvdm = pr.parsePDB(workdir + 'W_A-10-A-14-A-139_H-D-H_316-397-352_BestOPscore_BestGeo_BestCluscore_Best2ShScore_Best2ShRmsd.pdb')

gvdm = pr.parsePDB(workdir + 'W_A-10-A-14-A-139_H-D-H_316-397-352_A139_hie_THR_A_11_0.5_2.44.pdb')


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

target = pr.parsePDB(workdir + 'helix6_10-14-139_lig318_prep.pdb')

gvdm_A62 = pr.parsePDB(workdir + 'Lig-318_phenol_A62_ASP_rmsd_0.74_v_2.8_1293.pdb.gz')
gvdm_A65 = pr.parsePDB(workdir + 'Lig-318_ph_A65_PHE_rmsd_0.34_v_2.3_19174.pdb.gz') 
gvdm_A84 = pr.parsePDB(workdir + 'Lig-318_bb_cco_A84_SER_rmsd_0.61_v_2.6_4224.pdb.gz')


target.select('resnum 60 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '2HD', '1HE', '2HE']
his_order = ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']
coord_a = np.array([target.select('resnum 60 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid C and resnum 77 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 60 and name CB CG CD2 ND1 CE1 NE2 2HD 1HE 2HE').copy()
tf.apply(ag)
target.select('resnum 60 and name CB CG CD2 ND1 CE1 NE2 2HD 1HE 2HE').setCoords(ag.getCoords())
target.select('resnum 60 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD2', 'HE1', 'HE2'])


target.select('resnum 90 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '1HD', '2HD', '1HE']
coord_a = np.array([target.select('resnum 90 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid A and resnum 198 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 90 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').copy()
tf.apply(ag)
target.select('resnum 90 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').setCoords(ag.getCoords())
target.select('resnum 90 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1'])

target.select('resnum 94 and sc').getNames() #['CB', 'CG', 'OD1', 'OD2', 'HA', '2HB', '3HB']
mvdm.select('chid B and resnum 19 and heavy and sc').getNames()
target.select('resnum 94 and sc and heavy').setCoords(mvdm.select('chid B and resnum 19 and heavy and sc').getCoords())
target.select('resnum 94 and sc').setNames(['CB', 'CG', 'OD1', 'OD2', 'HA', 'HB1', 'HB2'])

#>>> mvdm2ndshell.
target.select('resnum 91 and sc').getNames()
#['CB', 'CG2', 'OG1', 'HA', 'HB', '1HG', '1HG2', '2HG2', '3HG2']
coord_b = np.array([gvdm.select('chid X and resnum 10 and name ' + x).getCoords()[0] for x in ['CB',  'CG2', 'OG1', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23']])
target.select('resnum 91 and sc').setCoords(coord_b)
target.select('resnum 91 and sc').setNames(['CB',  'CG2', 'OG1', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23'])

#>>> gvdms
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