def mp_handler(pairs, ss_dis_dict, scratch_dir, pdb_folder, result_file, np_pool):
    '''
    Multi process handler. Start a pool of processes and queues them
    on a number of cores. Runs the mp_worker and prints the output sequentially.
    '''

    if not args.pbs_range:
        # Print output header:
        print_header(result_file)
    else:
        run_range = args.pbs_range.split('-')
        run_range = list(map(int, run_range))

    # Generate a list of info that is needed for each pair
    # to run indepedently by the worker process:
    pair_info_list = list()
    n_comparisons = 0
    for pair_number, pair_tuple in enumerate(pairs):
        if pair_number < run_range[0] or pair_number > run_range[1]:
            continue
        pair_info_list.append([pair_number, pair_tuple, ss_dis_dict, scratch_dir, pdb_folder])
        n_comparisons += len(pair_tuple[0].split('-')) * len(pair_tuple[1].split('-'))
    n_pairs = len(pair_info_list)

    # Start the pool with X cores:
    pool = multiprocessing.Pool(np_pool)

    # Continue printing the results from each process in the pool:
    with open(result_file, 'a') as fh:
        completed = 0
        started = 0
        count = 0
        for result in pool.imap_unordered(mp_worker, pair_info_list):
            print_str, s, c = result
            completed += c
            started += s
            # Print stats:
            count += 1
            if count % 10 == 0:
                print('### {} out of {} pair tuples have been run'.format(count, n_pairs))
                print('### Currently {} out of {} comparisons have run, {} with success.'.format(started, n_comparisons, completed))
            if print_str:
                fh.write(print_str)

    pool.close()
    pool.join()




def pbs_submission(nsubs, njobs, flags):
    jobs_per_sub = round(njobs / nsubs + 0.5)
    script_path = run_dir + '/' + os.path.abspath(__file__).split('/')[-1]
    for i in range(nsubs):
        log_err = 'pair_log' + str(i) + '.err'
        log_out = 'pair_log' + str(i) + '.out'
        run_range = str(i * jobs_per_sub) + '-' + str((i + 1) * jobs_per_sub - 1)

        s1 = '#!/bin/sh\n\
### Note: No commands may be executed until after the #PBS lines\n\
### Account information\n\
#PBS -W group_list=cu_10020 -A cu_10020\n\
### Job name (comment out the next line to get the name of the script used as the job name)\n\
#PBS -N krdav_job\n\
### Output files (comment outa the next 2 lines to get the job name used instead)\n' + '#PBS -e ' + log_err + '\n' + '#PBS -o ' + log_out + '\n' + '### Email: no (n)\n\
#PBS -M n\n\
### Make the job rerunable (y)\n\
#PBS -r y\n\
### Number of nodes\n\
#PBS -l nodes=1:ppn=' + str(args.np) + ':thinnode\n\
### Requesting time - 12 hours - overwrites **long** queue setting\n\
#PBS -l walltime=24:00:00\n\
\n\
echo This is the STDOUT stream from a PBS Torque submission script.\n\
# Go to the directory from where the job was submitted (initial directory is $HOME)\n\
echo Working directory is $PBS_O_WORKDIR\n\
cd $PBS_O_WORKDIR\n\
\n\
# Load user Bash settings:\n\
source /home/people/krdav/.bash_profile\n\
\n\
echo Modules and user Bash settings was loaded.\n\
\n\
module unload anaconda2/4.0.0\n\
module load anaconda3/4.0.0\n\
\n\
echo  Now the user defined script is run. After the ---- line, the STDOUT stream from the script is pasted.\n\
echo -----------------------------------------------------------------------------------------------------\n\
# Run the desired script:\n\
python ' + script_path + ' ' + flags + ' -pbs_range ' + run_range + '\n'

        qsub_name = run_dir + '/sub' + str(i) + '.qsub'
        with open(qsub_name, 'w') as fh:
            print(s1, file=fh)

        cmd = 'qsub {}'.format(qsub_name)
        os.system(cmd)