'''
After search the ~500 benchmark protein, the next step is to check the different score methods.

The f_median = np.prod([overlaps])/np.prod([median clu of each aa_type])*total_clu/volume.

Two things need to be checked: 
1. Distribution of f_median for all protein in different radius.
2. Type of binding core. HIS-HIS-HIS vs HIS-HIS-GLU ... etc

The function here just combine info from all benchmark cores.
'''

import os
import re
import pandas as pd

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/_Seq_core_date_3contact/20211028_eval_selfcenter/'

outfile = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/_Seq_core_date_3contact/_all_log_info.tsv'


aa_vdm_info_dict = {
            'HIS':[2662, 285, 92.08, 60],
            'GLU':[829, 50, 16.32, 11],
            'ASP':[896, 90, 26.12, 22]
        }

df_all = pd.DataFrame()

for dir in os.listdir(workdir):
    log = workdir + 'output_selfcenter_eval1028_1982_4tln_ZN_1.pdb_20211028-162254/_log.txt' 
    log = workdir + dir + '/_log.txt'
    if not os.path.exists(log):
        continue
    title = dir.split('.')[0]
    try:
            # with open(log, 'r') as ff:
            #     lines = ff.readlines()
            #     key = lines[1].split('\t')[0]
            #     core_types = list(filter(None, re.split("[ (,)']", key)))
            #     core_type = '-'.join(sorted([core_types[3],core_types[5],core_types[7]]))
            #     for line in lines[1:]:
            #         f.write(title + '\t' + core_type + '\t' + line)
        df = pd.read_csv(log, sep='\t', lineterminator='\n', index_col=False)
        extract_core_types = list(filter(None, re.split("[ (,)']", df["key"][0]))) 
        core_types = [extract_core_types[3],extract_core_types[5],extract_core_types[7]]

        max_fmedian = df["f_median"].max()
        if max_fmedian == 0:
            print(title)
            max_fmedian = 1
        df['n_f_median'] = [x/max_fmedian for x in df["f_median"]]

        totals = [aa_vdm_info_dict[c][0] for c in core_types]
        df["total1"] = totals[0]
        df["total2"] = totals[1]
        df["total3"] = totals[2]
        maxs = [aa_vdm_info_dict[c][1] for c in core_types]
        df["max1"] = maxs[0]
        df["max2"] = maxs[1]
        df["max3"] = maxs[2]
        avgs = [aa_vdm_info_dict[c][2] for c in core_types]
        df["avg1"] = avgs[0]
        df["avg2"] = avgs[1]
        df["avg3"] = avgs[2]
        medians = [aa_vdm_info_dict[c][3] for c in core_types]
        df["median1"] = medians[0]
        df["median2"] = medians[1]
        df["median3"] = medians[2]
        
        df["clu2total"] = df["clu1"]*df["clu2"]*df["clu3"]/df["total1"]/df["total2"]/df["total3"]
        df["clu2max"] = df["clu1"]*df["clu2"]*df["clu3"]/df["max1"]/df["max2"]/df["max3"]


        df["title"] = title
        df["core_type"] = '-'.join(sorted(core_types))

        df_all = df_all.append(df)
    except:
        print(title)

df_all.to_csv(outfile, sep='\t')


