'''
The 2ndshell vdm counting.
'''



'''
The method is to extract 2nd shell vdm by heavy atom distance cut off 3.4, then count the types of extracted 2nd shell vdms.  
'''
import os

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/"

for _folder1 in os.listdir(workdir):
    if not os.path.isdir(workdir + _folder1) and not '20220116_selfcenter' in _folder1:
        continue
    
    for _folder2 in os.listdir(workdir + _folder1 + '/'):
        if not os.path.isdir(workdir + _folder1 + '/' + _folder2) or 'AA2ndS' not in _folder2 or 'reps' not in _folder2:
            continue
        filesFromSame1stShellVdM = {}
        for file in os.listdir(workdir + _folder1  + '/' + _folder2 + '/'):
            splits = file.split('_')
            refile = splits[0] + '_' + splits[1] + '_' + splits[2] + '_' + splits[3] + '_' + splits[5]
            if refile in filesFromSame1stShellVdM.keys():
                filesFromSame1stShellVdM[refile] += 1
            else:
                filesFromSame1stShellVdM[refile] = 1
        with open(workdir + 'reason/' +  _folder1 + '_' + _folder2 + '_summary.tsv', 'w') as f:
            for key in filesFromSame1stShellVdM.keys():
                f.write(key + '\t' + str(filesFromSame1stShellVdM[key]) + '\n')
        print(_folder1 + '-' + _folder2 + '-' +  str(len(filesFromSame1stShellVdM.keys())) + '-' + str(sum(filesFromSame1stShellVdM.values())))


'''
New method is to use Probe to calc the bonding of protein, then count the 'hb' which are 2nd shell. 
'''
import os

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_2ndshell_date_reps_probe2ndshell/"

for _folder1 in os.listdir(workdir):
    if not os.path.isdir(workdir + _folder1):
        continue
    total = 0
    total_bb = 0
    for file in os.listdir(workdir + _folder1 + '/'):
        total += 1
        if '_bb_' in file:
            total_bb += 1
    print(_folder1 + ' total: ' + str(total) + ' bb: ' + str(total_bb))    



