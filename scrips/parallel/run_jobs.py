import os
import shutil
import subprocess


def run(data_set_path, script_name, script_arguments, num_jobs, time='48:00:00', mem_free_GB=4, scratch_space_GB=1,
        architecture='linux-x64', hold_jid=None, keep_job_output_path=True):

    job_output_path = os.path.join(data_set_path, "job_outputs")
    
    if not keep_job_output_path and os.path.exists(job_output_path): # Clear the job output path
        shutil.rmtree(job_output_path)

    if not os.path.exists(job_output_path):
        os.mkdir(job_output_path)

    hj = ['-hold_jid', hold_jid] if hold_jid else []

    qsub_command = ['qsub',
                    '-cwd'] + hj \
                    + ['-N', script_name.split('/')[-1],
                    '-t', '1-{0}'.format(num_jobs),
                    '-l', 'h_rt={0}'.format(time),
                    '-l', 'mem_free={0}G'.format(mem_free_GB),
                    '-l', 'scratch={0}G'.format(scratch_space_GB),
                    '-l', 'arch=linux-x64',
                    '-o', job_output_path,
                    '-e', job_output_path,
                    './job_scripts/run_SGE_job.sh',
                    script_name,
                    data_set_path] \
                    + script_arguments \
                    + [num_jobs]

    print(qsub_command)
    
#    subprocess.check_call(qsub_command)



if __name__=='__main__':

    #workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/data/zn_eval_bench_mark/2013_2014/'
    workdir = workdir = '/mnt/e/DesignData/ligands/LigandBB/zn_eval_bench_mark/2013_2014/'
    # target_files = []
    # for target_file in os.listdir(workdir):
    #     if target_file.endswith('.pdb'):
    #         summaryfile = workdir + 'output_dist45_phipsi30_' + target_file.split('.')[0] + '/_summary.tsv'
    #         if os.path.isfile(summaryfile):
    #             continue
    #         target_files.append(target_file)

    
    run(workdir, 'run_parallel_eval_search.py', [workdir, ''], 1)




