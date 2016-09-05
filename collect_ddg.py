
import glob
import sys
import os
import shutil
import re
import gzip
import pickle
import datetime
import time
import math
import contextlib
import multiprocessing
import argparse

sys.path.insert(0, '/home/projects/cu_10020/apps/python3-site-packages/lib/python/')
from Bio.PDB import *


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
    "SEP": "S",  # Phosphoserine
    "THR": "T",
    "TPO": "T",  # Phosphothreonine
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "PTR": "Y",   # Phosphotyrosine
    "DAL": "A",   # Dextro from here
    "DCY": "C",
    "DAP": "D",
    "DGU": "E",
    "DPH": "F",
    "DGL": "G",
    "DHI": "H",
    "DIL": "I",
    "DLY": "K",
    "DLE": "L",
    "DME": "M",
    "DMS": "M",
    "DAN": "N",
    "DPR": "P",
    "DGN": "Q",
    "DAR": "R",
    "DSE": "S",
    "DTH": "T",
    "DVA": "V",
    "DTR": "W",
    "DTY": "Y"
}


AAletters = 'ACDEFGHIKLMNPQRSTVWY'

D_isomer_AA = ["DAL", "DCY", "DAP", "DGU", "DPH",
               "DGL", "DHI", "DIL", "DLY", "DLE",
               "DME", "DMS", "DAN", "DPR", "DGN",
               "DAR", "DSE", "DTH", "DVA", "DTR",
               "DTY"]

# Don't try to add other ptm's like acetylation etc.
# It is not worth it.
ptm_residues = ["SEP", "TPO", "PTR"]


whitelist = ['MSE']

whitelist.extend(D_isomer_AA)
whitelist.extend(ptm_residues)


# Currently nothing in blacklist
# The badly modified residues are normally just cut out of the crystal
blacklist = []

# Remove from crystal structure:
greylist = ["HOH"]


def ddg_success(ddg_file):
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
### Don't do anything (False)
# Code 2: Unforseen error
### Response: Report a problem, and don't try to run it (False)
# Code 3: The job is still running
### Response: Let it run until it finishes (False)
# Code 4: The job was killed for exceed run time limits
### Response: Don't do anything (False)
# Code 5: The job ended with success
### Response: Collect ddG values (True)
def ddg_choice(ddg_response):
    if ddg_response == 1:
        return(False)
    elif ddg_response == 2:
        return(False)
    elif ddg_response == 3:
        return(False)
    elif ddg_response == 4:
        return(False)
    elif ddg_response == 5:
        return(True)
    else:
        if args.verbose:
            print('Unforseen error. ddg_success could not be determined.')
        return(False)


def cst_min_success(folder):
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
    with open(out_log) as fh:
        lines = fh.readlines()
        # Apparently the cst_min app ends with this on success:
    if lines[-1].startswith('running another iteration of minimization'):
        return(True)
    else:
        if args.verbose:
            print('The cst_min application did not end succesfully in folder:\n', folder)
        return(False)


def write_ddG_results(folder, ddG_dict):
    resname = folder + '/ddG_results.txt'
    with open(resname, 'w') as fh_out:
        for name, res_list in sorted(ddG_dict):
            if len(name) == 1:
                print('# Monomeric stability for chain', name, file=fh_out)
                print('# ChainID Residue_number InsertionCode\tFromTo\tREU', file=fh_out)
                tab_list = ['\t'.join(t[1:3]) for t in res_list]
                print('\n'.join(tab_list), file=fh_out)
            elif len(name) == 2:
                print('# Complex stability for chain complex', name, file=fh_out)
                print('# ChainID Residue_number InsertionCode\tFromTo\tREU', file=fh_out)
                tab_list = ['\t'.join(t[1:3]) for t in res_list]
                print('\n'.join(tab_list), file=fh_out)


def check_mapping(mapping_dict):
    # From 'PRO A   1 '
    # To 'PRO A   1 '
    new_mapping = dict()
    for k, v in mapping_dict:
        if k[0:3] != v[0:3]:
            print(k, v)
            sys.exit()
        new_mapping[k[3:]] = v[3:]
    return(new_mapping)


def collect_ddg(ddg_file, mapping_dict, ddG_dict):
    name = ddg_file.split('/')[-1][12:-9]
    ddg_rundir = '/'.join(ddg_file.split('/')[0:-1]) + '/' + name + '_ddg_rundir'
    ddg_outfile = ddg_rundir + '/' + 'ddg_predictions.out'
    ori_AAseq = get_AA_string(ddg_rundir, name)
    os.mkdir(ddg_rundir)
    mapping_dict = check_mapping(mapping_dict)
    # Create a new dict for the chain:
    ddG_dict[name] = list()
    offset = 0
    with open(ddg_outfile) as fh:
        lines = fh.readlines()
        for line in lines:
            if offset > 20:
                print("Offset is very high. Many residues have been skipped. Terminating:", ddg_outfile)
                sys.exit()
            line = line.strip()
            if not line:
                continue
            elif line.startswith('ddG: description'):
                continue
            rows = line.split()
            AA = rows[1][1]
            res_count = int(rows[1][1:-1])
            # Check for skipped residue:
            while ori_AAseq[res_count - 1 + offset] != AA or offset > 20:
                offset += 1
            res_idx_string = 'A{:>4} '.format(res_count + offset)  # Notice chain A and the blank insertion code
            ori_res_idx_string = mapping_dict[res_idx_string]
            from_to = rows[1][1] + '->' + rows[1][-1]
            res_tuple = (ori_res_idx_string, from_to, rows[2])
            ddG_dict[name].append(res_tuple)
    return(ddG_dict)


def get_AA_string(ddg_rundir, name):
    fnam = ddg_rundir + '/' + name
    with open(fnam) as fh:
        lines = fh.readlines()
        prev_resnumb = '   0'
        AA_string = ''
        for line in lines:
            resnumb = line[22:26]
            if resnumb != prev_resnumb:
                res_name = line[17:20]
                AA = residue_type_3to1_map[res_name]
                AA_string += AA
    return(AA_string)


if __name__ == "__main__":
    if not args.check_all_folders:
        with open(args.folders_for_ddg) as fh:
            folder_list = fh.read().splitlines()
    else:
        folder_list_glob = args.db_split_dir + '/*/*/'
        folder_list = glob.glob(folder_list_glob)
        folder_list = [f.rstrip('/') for f in folder_list]

    for idx, folder in enumerate(folder_list):
        if not cst_min_success(folder):
            continue
        mapping_path = folder + '/' + 'residue_mapping_dict.p'
        mapping_dict = pickle.load(open(mapping_path, "rb"))
        ddg_file_glob = folder + '/' + '*.pdb'
        files_for_ddg = glob.glob(ddg_file_glob)
        ddG_dict = dict()
        for ddg_file in files_for_ddg:
            ddg_response = ddg_success(ddg_file)
            if not ddg_choice(ddg_response):
                continue
            ddG_dict = collect_ddg(ddg_file, mapping_dict, ddG_dict)
        write_ddG_results(folder, ddG_dict)











