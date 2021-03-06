import glob
import sys
import os
import shutil
import time
import argparse


# Build commandline parser
parser = argparse.ArgumentParser(description="Run the ddg_monomer in Rosetta.")

# Arguments
parser.add_argument(
    "-db_home_dir",
    "--db_home_dir",
    type=str,
    dest="db_home_dir",
    metavar="PATH",
    help="Absolute path to the run directory. Current directory by default.",
    required=False)
parser.add_argument(
    "-db_split_dir",
    "--db_split_dir",
    type=str,
    dest="db_split_dir",
    metavar="PATH",
    help="Absolute path to the split directory containing all the two character folders, pdb style. By default the directory \"split\" in db_home_dir.",
    required=False)
parser.add_argument(
    "-folders_for_ddg",
    "--folders_for_ddg",
    type=str,
    dest="folders_for_ddg",
    metavar="FILE",
    help="File with a list of paths for input folders in the database for which to make ddg calculations.",
    required=False)
parser.add_argument(
    "-a",
    "--check_all_folders",
    type=str,
    dest="check_all_folders",
    metavar="SWITCH",
    help="Should all the created folders in the database be cheched for ddG calculations?",
    required=False)
parser.add_argument(
    "-rosetta_db",
    "--rosetta_db",
    type=str,
    dest="rosetta_db",
    metavar="PATH",
    help="Absolute path to the Rosetta database.",
    required=False)
parser.add_argument(
    "-rosetta_ddg_app",
    "--rosetta_ddg_app",
    type=str,
    dest="rosetta_ddg_app",
    metavar="PATH",
    help="Absolute path to the Rosetta dgg application.",
    required=False)
parser.add_argument(
    "-v",
    "--verbose",
    type=int,
    dest="verbose",
    metavar="VERBOSE",
    help="Should the script be verbose?",
    required=False)
parser.add_argument(
    "-ve",
    "--verbose_error",
    type=int,
    dest="verbose_error",
    metavar="SWITCH",
    help="Should the script be verbose on error only? On by default.",
    required=False)
parser.add_argument(
    "-restart_failed",
    "--restart_failed",
    type=int,
    dest="restart_failed",
    metavar="SWITCH",
    help="Try to restart those jobs that failed. Default 0.",
    required=False)

# Set default arguments
parser.set_defaults(
    rosetta_db='/services/tools/rosetta/2016.10/main/database',
    rosetta_ddg_app='/services/tools/rosetta/2016.10/main/source/bin/ddg_monomer.default.linuxgccrelease',
    verbose_error=1,
    restart_failed=0)

# Put arguments into args object:
args = parser.parse_args()


# Determine what to do with the input arguments:
if not args.db_home_dir:
    args.db_home_dir = os.getcwd()
else:
    args.db_home_dir = args.db_home_dir.rstrip('/')
    os.chdir(args.db_home_dir)

if not args.db_split_dir:
    args.db_split_dir = args.db_home_dir + '/split'
else:
    args.db_split_dir = args.db_split_dir.rstrip('/')

run_dir = os.getcwd()
if args.folders_for_ddg and not args.check_all_folders:
    pass
elif not args.folders_for_ddg and args.check_all_folders:
    args.folders_for_ddg = run_dir + '/folders_for_ddg.txt'
else:
    parser.error('No action requested, add --folders_for_ddg or --check_all_folders')

if args.folders_for_ddg[0] != '/':
    args.folders_for_ddg = run_dir + '/' + args.folders_for_ddg


# Global variables:
const_flags_ddg = '-resfile resfile.txt -ddg:weight_file soft_rep_design -ddg::iterations 1 -ddg::dump_pdbs false -ignore_unrecognized_res -ddg::local_opt_only true -ddg::suppress_checkpointing true -in::file::fullatom -ddg::mean true -ddg::min false -mute all -ddg::output_silent false'
np = 1

residue_type_3to1_map = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "MSE": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y"
}

AAletters = 'ACDEFGHIKLMNPQRSTVWY'

D_isomer_AA = ["DAL", "DCY", "DAP", "DGU", "DPH",
               "DGL", "DHI", "DIL", "DLY", "DLE",
               "DME", "DMS", "DAN", "DPR", "DGN",
               "DAR", "DSE", "DTH", "DVA", "DTR",
               "DTY"]

# Don't try to add other ptm's like acetylation etc.
# It is not worth it.
whitelist = list(residue_type_3to1_map.keys())
whitelist.append('MSE')
ptm_residues = ["SEP", "TPO", "PTR"]
whitelist.extend(D_isomer_AA)
whitelist.extend(ptm_residues)


def pbs_submit_cmd(np, cmd_flags, run_dir, idx):
    '''Submit a job to the queue.'''
    os.chdir(run_dir)
    log_err = run_dir + '/sub' + str(idx) + '_log.err'
    log_out = run_dir + '/sub' + str(idx) + '_log.out'

    qsub_string = '#!/bin/sh\n\
### Note: No commands may be executed until after the #PBS lines\n\
### Account information\n\
#PBS -W group_list=pr_46690 -A pr_46690\n\
### Job name (comment out the next line to get the name of the script used as the job name)\n\
#PBS -N krdav_job\n\
### Output files (comment outa the next 2 lines to get the job name used instead)\n' + '#PBS -e ' + log_err + '\n' + '#PBS -o ' + log_out + '\n' + '### Email: no (n)\n\
#PBS -M n\n\
### Make the job rerunable (y)\n\
#PBS -r y\n\
### Number of nodes\n\
#PBS -l nodes=1:ppn=' + str(np) + ':thinnode\n\
#PBS -l walltime=100:00:00\n\
echo This is the STDOUT stream from a PBS Torque submission script.\n\
# Go to the directory from where the job was submitted (initial directory is $HOME)\n\
echo Working directory is $PBS_O_WORKDIR\n\
cd $PBS_O_WORKDIR\n\
\n\
# Load user Bash settings:\n\
source /home/people/krdav/.bash_profile\n\
module unload anaconda2/4.0.0\n\
module load anaconda3/4.0.0\n\
\n\
echo  Now the user defined script is run. After the ---- line, the STDOUT stream from the script is pasted.\n\
echo -----------------------------------------------------------------------------------------------------\n\
# Run the desired script:\n\
' + cmd_flags + '\n'

    qsub_path = run_dir + '/sub' + str(idx) + '.qsub'
    with open(qsub_path, 'w') as fh_out:
        print(qsub_string, file=fh_out)

    cmd = 'qsub {}'.format(qsub_path)
    os.system(cmd)
    time.sleep(1)  # Added a small time delay to let the queueing system work


def submit_ddg(pdb_file_path, idx):
    '''Prepare for and submit a ddg_monomer run to the queue.'''
    # Hacky way of getting the chain/pair name:
    name = pdb_file_path.split('/')[-1][12:-9]
    ddg_rundir = '/'.join(pdb_file_path.split('/')[0:-1]) + '/' + name + '_ddg_rundir'
    # Make the rundir:
    os.mkdir(ddg_rundir)
    dst_pdb_file_path = ddg_rundir + '/' + name
    # Copy the relevant pdb file:
    shutil.copy(pdb_file_path, dst_pdb_file_path)
    # Copy the resfile to the ddg_rundir:
    resfile_path = '/'.join(pdb_file_path.split('/')[0:-1]) + '/' + 'resfile.txt'
    dst_resfile_path = ddg_rundir + '/' + 'resfile.txt'
    shutil.copy(resfile_path, dst_resfile_path)
    # Create the run command:
    ddg_cmd = args.rosetta_ddg_app + ' ' + const_flags_ddg + ' -database ' + args.rosetta_db + ' -in:file:s ' + dst_pdb_file_path
    ddg_cmd += '\nrm mutant_traj* wt_traj*'  # Delete Rosetta annoying checkpoint files
    # Submit the job:
    print('Submitting for:', pdb_file_path)
    pbs_submit_cmd(np, ddg_cmd, ddg_rundir, idx)


def resubmit_ddg(pdb_file_path, idx):
    '''Resubmit a failed/terminated ddg monomer run to the queue.'''
    # Hacky way of getting the chain/pair name:
    name = pdb_file_path.split('/')[-1][12:-9]
    ddg_rundir = '/'.join(pdb_file_path.split('/')[0:-1]) + '/' + name + '_ddg_rundir'
    # Prepare for resubmission:
    prep_restart_job(ddg_rundir, name)
    dst_pdb_file_path = ddg_rundir + '/' + name
    # Create the run command:
    ddg_cmd = args.rosetta_ddg_app + ' ' + const_flags_ddg + ' -database ' + args.rosetta_db + ' -in:file:s ' + dst_pdb_file_path
    # Submit the job:
    print('Resubmitting for:', pdb_file_path)
    pbs_submit_cmd(np, ddg_cmd, ddg_rundir, idx)


def ddg_success(ddg_file):
    '''Investigate whether a job has already run and if so,
    get the success of the previously run ddg job.'''
    # Hacky way of getting the chain name from the cst_min pdb output name:
    name = ddg_file.split('/')[-1][12:-9]
    folder = '/'.join(ddg_file.split('/')[0:-1]) + '/' + name + '_ddg_rundir'
    ddg_outfile = folder + '/' + 'ddg_predictions.out'
    # If either the ddg_rundir or the ddg_predictions.out file is missing,
    # there cannot have been a successful ddg run:
    if not os.path.exists(folder):  # Code 1: No folder, or never started and the minimized structure is there
        if args.verbose:
            print('The ddg folder could not be found. Proceeding to submit this as a job. Folder:\n', folder)
        return(1)
    elif not os.path.exists(ddg_outfile):  # Code 1: No folder, or never started
        if args.verbose:
            print('The ddg_predictions.out file could not be found. Proceeding to submit this as a job. Folder:\n', folder)
        shutil.rmtree(folder)
        return(1)
    elif os.path.exists(ddg_outfile):  # The job has either fininshed succesfully or somthing went wrong
        # Glob to find error and out logs:
        out_log_glob_string = folder + '/' + 'sub*_log.out'
        out_log_glob = glob.glob(out_log_glob_string)
        err_log_glob_string = folder + '/' + 'sub*_log.err'
        err_log_glob = glob.glob(err_log_glob_string)
        if len(out_log_glob) > 1 and len(err_log_glob) > 1:  # Code 2: Unforseen error
            if args.verbose:
                print('This is weird... There are more than one ddG submission STDOUT log. Maybe clean folder first? Folder:\n', folder)
            # shutil.rmtree(folder)
            return(2)
        elif len(out_log_glob) == 0 and len(err_log_glob) == 0:  # Code 3: The job is still running
            if args.verbose:
                print('No log found, job must be in progress:\n', folder)
            # Job in progress log not yet created:
            return(3)
        else:  # Check if the logs indicate success
            out_log = out_log_glob[0]
            err_log = err_log_glob[0]
        # Parse the error log:
        with open(err_log) as fh:
            lines = fh.readlines()
        if not lines:  # If the error log is empty
            pass
        elif lines[-1].startswith('=>> PBS: job killed'):  # Code 4: The job was killed for exceed run time limits
            if args.verbose:
                print('Job was killed by the queueing system because of too much run time:\n', err_log)
            return(4)
        elif 'Aborted' in lines[-1]:  # Code 2: Unforseen error
            if args.verbose:
                print('Rosetta have thrown an error an aborted:\n', folder)
            return(2)
        # Parse the out log:
        with open(out_log) as fh:
            lines = fh.readlines()
            # Apparently the ddg:monomer app ends with this on success:
        if lines[-1].startswith(' Total weighted score'):  # Code 5: The job ended with success
            if args.verbose:
                print('Job has succesfully ended:\n', folder)
            return(5)
        else:  # Code 2: Unforseen error
            if args.verbose:
                print('The submission log indicates that the run did not end succesfully:\n', out_log)
            # shutil.rmtree(folder)
            return(2)


# Code 1: No folder, or never started
### Response: Proceed as normal (True)
# Code 2: Unforseen error
### Response: Report a problem, and don't try to run it (False)
# Code 3: The job is still running
### Response: Let it run until it finishes (False)
# Code 4: The job was killed for exceed run time limits
### Response: Restart the job from where it ended (True)
###### Action: Delete old logs, create a new resfile
# Code 5: The job ended with success
### Response: Nothing (False)
def ddg_choice(ddg_response, pdb_file_path, idx, folders_for_ddg2):
    '''Determine what to do following the return code of the ddg success.'''
    if ddg_response == 1:
        submit_ddg(pdb_file_path, idx)
        folders_for_ddg2.append(pdb_file_path)
    elif ddg_response == 2:
        folders_for_ddg2.append(pdb_file_path)
    elif ddg_response == 3:
        folders_for_ddg2.append(pdb_file_path)
    elif ddg_response == 4:
        folders_for_ddg2.append(pdb_file_path)
        resubmit_ddg(pdb_file_path, idx)
    elif ddg_response == 5:
        pass
    else:
        if args.verbose:
            print('Unforseen error. ddg_success could not be determined.')
    return(folders_for_ddg2)


### How an update resfile looks like:
# ALLAA
# start
# 1 - last_res - 1 A NATAA
def prep_restart_job(folder, name):
    '''Pepare for resubmission of an early terminated job.'''
    # Glob to find the submission logs:
    out_log_glob_string = folder + '/' + 'sub*_log.out'
    out_log_glob = glob.glob(out_log_glob_string)
    err_log_glob_string = folder + '/' + 'sub*_log.err'
    err_log_glob = glob.glob(err_log_glob_string)
    out_log = out_log_glob[0]
    err_log = err_log_glob[0]
    # Remove the submission logs:
    os.remove(out_log)
    os.remove(err_log)
    # If monomer the chain name is intacts, if dimer the chain name is renamed to A:
    if len(name) == 1:
        chain_name = name
    else:
        chain_name = 'A'
    # Read the predictions.out to find where the job ended:
    finished_res = dict()
    ddg_outfile = folder + '/' + 'ddg_predictions.out'
    with open(ddg_outfile) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 1:  # Empty lines
                pass
            elif cols[1] == 'description':  # Header line
                pass
            else:
                resnumb = int(cols[1][1:-1])
                finished_res[resnumb] = 1
    last_res = sorted(finished_res.keys())[-1]

    # Print the new resfile:
    resfile = folder + '/' + 'resfile.txt'
    with open(resfile, 'w') as fh_out:
        print('ALLAA\nstart', file=fh_out)
        # Make the resfile start at the second last residue (last_res - 1):
        print('1 - {} {} NATAA'.format(last_res - 1, chain_name), file=fh_out)


def cst_min_success(folder):
    '''Investigate whether the cst_min run was a success or not.'''
    # Glob to find the out log:
    out_log_glob_string = folder + '/' + 'sub*_log.out'
    out_log_glob = glob.glob(out_log_glob_string)
    if len(out_log_glob) > 1:
        if args.verbose:
            print('This is weird... There are more than one min_cst submission STDOUT log. Maybe clean folder first? Folder:\n', folder)
        return(False)
    elif len(out_log_glob) == 0:
        return(False)
    else:
        out_log = out_log_glob[0]
    # Parse the log:
    with open(out_log) as fh:
        lines = fh.readlines()
        # Apparently the cst_min app ends with this on success:
    if lines[-1].startswith('running another iteration of minimization'):
        return(True)
    else:
        if args.verbose:
            print('The cst_min application did not end succesfully in folder:\n', folder)
        return(False)


if __name__ == "__main__":
    # Either check all folders or a list provided:
    if not args.check_all_folders:
        with open(args.folders_for_ddg) as fh:
            folder_list = fh.read().splitlines()
    else:
        # If all folders should be checked then glob is used:
        folder_list_glob = args.db_split_dir + '/*/*/'
        folder_list = glob.glob(folder_list_glob)
        folder_list = [f.rstrip('/') for f in folder_list]

    job_idx = 0  # Index for the log, not strictly necessary
    folders_for_ddg2 = list()  # Folders that are still not finished
    for idx, folder in enumerate(folder_list):
        # If the cst_min failed skip it:
        if not cst_min_success(folder):
            continue
        # Find the cst_min structures:
        ddg_file_glob = folder + '/' + '*.pdb'
        files_for_ddg = glob.glob(ddg_file_glob)
        for ddg_file in files_for_ddg:
            ddg_response = ddg_success(ddg_file)
            folders_for_ddg2 = ddg_choice(ddg_response, ddg_file, job_idx, folders_for_ddg2)
            job_idx += 1

    # Write an updated "folders to update" list:
    folders_for_ddg2 = list({'/'.join(f.split('/')[:-1]): 1 for f in folders_for_ddg2}.keys())
    with open(args.folders_for_ddg, 'w') as fh_out:
        print('\n'.join(folders_for_ddg2), file=fh_out)
