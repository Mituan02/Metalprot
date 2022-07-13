# Notes/scripts to Design helix bundles around the ligand TTS.

# 1. Search mvdMs in helix bundles. 'run_selfcenter_search.py'.
     Mannually pair-fit ideal geometry, then pair-fit ligand tts.
# 2. Search lig vdMs in helix bundles. 'run_search_lig_indep.py'
# 3. Search loops. 'run_loop_ss.py'
# 4. Generate resfile. Use rosetta to generate all ala file. 
#    Then 'DesignScript/topology/topology_usage.py'. 
#    Manually fix the output resfile.
# 5. Create others files need for rosetta and run rosetta. 


## Run rosetta to generate all ala pdb.

cd /wynton/home/degradolab/lonelu/DesignData/Metalloenzyme/HelixFe_TTS/20220625_design/
~/rosetta_bin_linux_2020.08.61146_bundle/main/source/bin/fixbb.static.linuxgccrelease @options_all_ala_fixbb.txt -ignore_unrecognized_res  -ignore_zero_occupancy

cd /mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_10-14-139_HDH/Lig_318/design_0/
/mnt/d/Rosetta/main/source/bin/fixbb.static.linuxgccrelease @options_all_ala_fixbb.txt -ignore_unrecognized_res -ignore_zero_occupancy false

## Run rosetta to generate .params file.

cd /mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_10-14-139_HDH/Lig_318/design_0_rosetta/
/mnt/d/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n TTS -p TTS --conformers-in-one-file tts.sdf

## Tips.
option -qsar:grid_dir is not set.  Use this flag to specify a directory to store scoring grids.  This will save you a huge amount of time.
When build helix-ligand, remember to 'rebuild' before save.
Rename the ligand name to 'TTS' in the protein .pdb file, add FE1-O3 CONNECT. Each time, generate new TTS.params file.
When manual mutate vdMs on bb, try to change the atom names. 
Adjust the '.cst' file. 
Adjust the '.resfile'. 
Pre make the 'output' folder.

## Account.
Run AF2 colab: lonelu.af@gmail.com, lonelu.af2@gmail.com

## Dialog
In the design process, I use nick_helix6 to start, because it contains several good metal binding sites. Eventually, all the nick_helix_bundles are same on bottom 2helix when considering the binding of the ligand.

First design on position helix6_10-14-139_HDH, considering there is an 2ndshell for one of the M-His. Design followed: 1 find mvdm, 2 find 2ndshell, 3 find gvdm, 4 find loop, 5 construct. 

Initial design I just use pymol's mutation to mutate vdms, which generate shifted coords and wrong atom names. Still I could get some reasonable designs (E:/DesignData/Metalloenzyme/HelixFe_TTS/helix6_10-14-139_HDH/Lig_318/design_0_rosetta/output_0_sel/helix6_10-14-139_loop1_lig318_821, helix6_10-14-139_loop1_lig318_646) I used two different loop topologies. It turned out the longer one is better. The rosetta design protocols are copied from the previous ntf2 design protocols, with .comp and .charge are copied from SamM's tutorial.

I fixed the coords and atom names of the input and rerun on the old designs. While most of them could not maintain the metal geometry. While rerun on the previous reasonable rosetta generated designs could maintain the geometry.

I designed on helix6_14-135-139_HHE with correct configuration, and could get reasonable designs (E:/DesignData/Metalloenzyme/HelixFe_TTS/helix6_14-135-139_HHE/helix6_loop1_rosetta/output_2_sel/helix6_14-135-139_loop1_33). While rerun on the reasonable rosetta generated designs could maintain the geometry.

It turned out all designs contain too many GLU and ASN. I learned that the charge is better to be around -5 and the ratio of each aa in helix bundles. 
I rerun on the reasonable rosetta generated designs (helix6_10-14-139_loop1_lig318_821, helix6_14-135-139_loop1_33) with the new files. Here I also commented out the dihedral angle constrain on GLU and ASP.   
     
Eventually, the following two looks promising: helix6_10-14-139_loop1_lig318_821_811 and 
helix6_14-135-139_loop1_33.

I will order two of the designs.