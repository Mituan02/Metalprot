The file describes how to run prebuild Metalprot on wynton.

The program is located at: /wynton/home/degradolab/lonelu/GitHub_Design/Metalprot

The databases, scripts are located at: /wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/data/

Using 'run_selfcenter_one_structre.py' as an example:
1. using development mode in wynton. > ssh dev1
2. get into the folder. cd /wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/
3. Activate the conda environment. > conda activate env_conda
4. run line14 'python /wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/data/run_selfcenter_one_structure.py'.

Inside the scripts:
1. query_dir (line_17) is the metal-vdMer database. (line_19 - 28 load the database)
2. workdir (line_34) is the working directory contains the target protein backbone. 
3. outdir is the where the output will be.
4. target_path is the target protein backbone.
5. win_filter is the positions on protien backbone where vdMers will be searched.
6. rmsd_cuts is the metal-metal distance cutoff to find a combination binding core.
7. selfcenter_rmsd is the density sphere radius.
8. num_iters is the number of vdMers for each binding core.
9. _filter is the filters implemented. In the _filter, their function agree with names:
   a. filter_abple
   b. filter_phipsi, max_phipsi_val
   c. filter_vdm_score, min_vdm_score 
   d. filter_vdm_count, min_vdm_clu_num
   e. after_search_filter: This one will filter the search results after nearest neighbor search.
      pair_angle_range
      pair_aa_aa_dist_range
      pair_metal_aa_dist_range (If set to None, it won't filter anything)
      filter_qt_clash: this one will filter vdMer-backbone clash or vdMer-vdMer clash.
   f. write_filtered_result: do you want to keep after_search_filtered results.
   g. selfcenter_filter_member_phipsi: filter member vdmers for density calculation.
10. Search_selfcenter is the class for searching. In the Search_selfcenter:
   a. validateOriginStruct = True will use the protein backbone's sequence. Turn it to false to allow mutations.
   b. parallel = True. will use parallel for different combination of positions.


To search your own structure, copy and rename the script 'run_selfcenter_one_structre.py' in your folder. change the workdir and other parameters accordingly. 
To search multiple sctuctures, use script 'run_selfcenter_multiple_structres.py'.
Please let me know if you have any questions. 



