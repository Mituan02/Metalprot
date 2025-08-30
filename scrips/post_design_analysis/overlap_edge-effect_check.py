import os
from re import T
import sys
import prody as pr
import numpy as np
import pickle

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_category/pickle_noCYS/'

with open(query_dir + 'all_vdms.pkl', 'rb') as f:
    all_vdms = pickle.load(f)

workdir = '/mnt/e/DesignData/ligands/LigandBB/6dwv/AAMetalPhiPsi_HIS_cluster05/'

w4 = pr.parsePDB(workdir + 'win_4_6dwv_core.pdb')
w6 = pr.parsePDB(workdir + 'win_6_6dwv_core.pdb')
w15 = pr.parsePDB(workdir + 'win_15_6dwv_core.pdb')


for v in all_vdms: 
    if v.get_cluster_key()==('HIS', 10) or v.get_cluster_key() == ('HIS', 20): 
        selection = 'name N CA C and resindex 1'
        pr.calcTransformation(v.query.select(selection), w4.select(selection)).apply(v.query)
        pr.writePDB(workdir + v.query.getTitle() + '.pdb', v.query)
    if v.get_cluster_key()==('HIS', 4): 
        selection = 'name N CA C and resindex 1'
        pr.calcTransformation(v.query.select(selection), w6.select(selection)).apply(v.query)
        pr.writePDB(workdir + v.query.getTitle() + '.pdb', v.query)
    if v.get_cluster_key()==('HIS', 21): 
        selection = 'name N CA C and resindex 1'
        pr.calcTransformation(v.query.select(selection), w15.select(selection)).apply(v.query)
        pr.writePDB(workdir + v.query.getTitle() + '.pdb', v.query)


query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/AAMetalPhiPsi_HIS_cluster05/'

workdir = '/mnt/e/DesignData/ligands/LigandBB/6dwv/AAMetalPhiPsi_HIS_cluster05_selfcenter/'
os.makedirs(workdir, exist_ok=True)

for file in os.listdir(query_dir + '1713/'):
    pdb = pr.parsePDB(query_dir + '1713/' + file)
    pr.calcTransformation(pdb.select('heavy'), w4.select('heavy')).apply(pdb)
    pr.writePDB(workdir + file, pdb)

for file in os.listdir(query_dir + '1286/'):
    pdb = pr.parsePDB(query_dir + '1286/' + file)
    pr.calcTransformation(pdb.select('heavy'), w6.select('heavy')).apply(pdb)
    pr.writePDB(workdir + file, pdb)

for file in os.listdir(query_dir + '2305/'):
    pdb = pr.parsePDB(query_dir + '2305/' + file)
    pr.calcTransformation(pdb.select('heavy'), w15.select('heavy')).apply(pdb)
    pr.writePDB(workdir + file, pdb)