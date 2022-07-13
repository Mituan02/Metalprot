import prody as pr
import numpy as np


#>>> mvdms.
#>>> HIS/ASP sc can be transformed. 
#>>> GLU can not be transformed if including CB, one way is not including CB. Another way is to find a similar one and only change the coord of sc heavy atoms.


workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_14-135-139_HHE/'

mvdm = pr.parsePDB(workdir + 'structure_prep/W_A-14-A-135-A-139_H-H-E_351-334-451_BestOPscore_BestCluscore_Best2ShScore_Best2ShRmsd.pdb')

gvdm = pr.parsePDB(workdir + 'structure_prep/W_A-14-A-135-A-139_H-H-E_351-334-451_A14_hie_THR_A_136_0.53_2.44.pdb')



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

target = pr.parsePDB(workdir + 'helix6_loop2/helix6_loop2_prep.pdb')

gvdm_A58 = pr.parsePDB(workdir + 'helix6_loop1/Lig-157_ph_A58_PHE_rmsd_0.33_v_1.1_25862.pdb.gz')
gvdm_A65 = pr.parsePDB(workdir + 'helix6_loop1/Lig-157_coo_A65_ASN_rmsd_0.67_v_1.0_3738.pdb.gz') 
gvdm_A85 = pr.parsePDB(workdir + 'helix6_loop1/Lig-157_phenol_A85_ASN_rmsd_0.69_v_1.1_2414.pdb.gz')

target.select('resnum 78 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '1HD','2HD', '1HE']
his_order = ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']
coord_a = np.array([target.select('resnum 78 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid B and resnum 119 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 78 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').copy()
tf.apply(ag)
target.select('resnum 78 and name CB CG CD2 ND1 CE1 NE2 1HD 2HD 1HE').setCoords(ag.getCoords())
target.select('resnum 78 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD1','HD2', 'HE1'])

target.select('resnum 45 and sc').getNames() #['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', '2HB', '3HB', '2HD', '1HE', '2HE']
coord_a = np.array([target.select('resnum 45 and name ' + x).getCoords()[0] for x in his_order])
coord_b = np.array([mvdm.select('chid A and resnum 77 and name ' + x).getCoords()[0] for x in his_order])
tf  = pr.calcTransformation(coord_a, coord_b)
ag = target.select('resnum 45 and name CB CG CD2 ND1 CE1 NE2 2HD 1HE 2HE').copy()
tf.apply(ag)
target.select('resnum 45 and name CB CG CD2 ND1 CE1 NE2 2HD 1HE 2HE').setCoords(ag.getCoords())
target.select('resnum 45 and sc').setNames(['CB', 'CG', 'CD2', 'ND1', 'CE1', 'NE2', 'HA', 'HB2', 'HB3', 'HD2', 'HE1', 'HE2'])

target.select('resnum 82 and sc and heavy').setCoords(mvdm.select('chid C and resnum 77 and heavy and sc').getCoords())
target.select('resnum 82 and sc').getNames() #['CB', 'CG', 'CD', 'OE1', 'OE2', 'HA', '2HB', '3HB', '2HG', '3HG']
target.select('resnum 82 and sc').setNames(['CB', 'CG', 'CD', 'OE1', 'OE2', 'HA', 'HB2', 'HB3', 'HG2', 'HG3'])

#>>> mvdm2ndshell.
target.select('resnum 79 and sc').getNames()
#['CB', 'CG2', 'OG1', 'HA', 'HB', '1HG', '1HG2', '2HG2', '3HG2']
coord_b = np.array([gvdm.select('chid X and resnum 10 and name ' + x).getCoords()[0] for x in ['CB',  'CG2', 'OG1', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23']])
target.select('resnum 79 and sc').setCoords(coord_b)
target.select('resnum 79 and sc').setNames(['CB',  'CG2', 'OG1', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23'])

#>>> gvdms
target.select('resnum 12 and sc').getNames()
#'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HA', '2HB', '3HB', '1HD', '2HD', '1HE', '2HE', 'HZ']
gvdm_A58.select('chid X and sc and resnum 10').getNames()
coord_b = gvdm_A58.select('chid X and sc and resnum 10').getCoords()
target.select('resnum 12 and sc').setCoords(coord_b)
target.select('resnum 12 and sc').setNames(gvdm_A58.select('chid X and sc and resnum 10').getNames())

target.select('resnum 19 and sc').getNames()
#['CB', 'CG2', 'OG1', 'HA', 'HB', '1HG', '1HG2', '2HG2', '3HG2']
coord_b = np.array([gvdm_A65.select('chid X and sc and resnum 10 and name ' + x).getCoords()[0] for x in ['CB', 'CG', 'ND2', 'OD1', 'HA', 'HB2', 'HB3', 'HD21', 'HD22']])
target.select('resnum 19 and sc').setCoords(coord_b)
target.select('resnum 19 and sc').setNames(['CB', 'CG', 'ND2', 'OD1', 'HA', 'HB2', 'HB3', 'HD21', 'HD22'])

target.select('resnum 105 and sc').getNames()
#['CB', 'CG2', 'OG1', 'HA', 'HB', '1HG', '1HG2', '2HG2', '3HG2']
coord_b = np.array([gvdm_A85.select('chid X and sc and resnum 10 and name ' + x).getCoords()[0] for x in ['CB', 'CG', 'ND2', 'OD1', 'HA', 'HB2', 'HB3', 'HD21', 'HD22']])
target.select('resnum 105 and sc').setCoords(coord_b)
target.select('resnum 105 and sc').setNames(['CB', 'CG', 'ND2', 'OD1', 'HA', 'HB2', 'HB3', 'HD21', 'HD22'])

pr.writePDB(workdir + target.getTitle() + '_fix.pdb', target)