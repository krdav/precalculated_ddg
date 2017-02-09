import glob
import sys
import os
import shutil
import pickle
import multiprocessing
import argparse
import numpy as np
import statistics
import json
import copy

# Build commandline parser
parser = argparse.ArgumentParser(description="Collect the results from the Rosetta ddg_monomer app.")

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
    type=int,
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

# Set default arguments
parser.set_defaults(
    check_all_folders=1)

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


def ddg_success(ddg_file):
    '''Investigate the success of the ddG run of a given protein chain.'''
    # Hacky way of getting the chain name from the cst_min pdb output name:
    name = ddg_file.split('/')[-1][12:-9]  # Chain name
    folder = '/'.join(ddg_file.split('/')[0:-1]) + '/' + name + '_ddg_rundir'  # ddG rundir for the chain
    ddg_outfile = folder + '/' + 'ddg_predictions.out'  # ddG predction results for the chain
    # If either the ddg_rundir or the ddg_predictions.out file is missing,
    # there cannot have been a successful ddg run:
    if not os.path.exists(folder):  # Code 1: No folder, or never started and the minimized structure is there
        if args.verbose:
            print('The ddg folder could not be found. Run the ddG application submission script. Folder:\n', folder)
        return(1)
    elif not os.path.exists(ddg_outfile):  # Code 1: No folder, or never started
        if args.verbose:
            print('The ddg_predictions.out file could not be found. Run the ddG application submission script again. Folder:\n', folder)
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
### Don't do anything (False)
# Code 2: Unforseen error
### Response: Don't do anything (False)
# Code 3: The job is still running
### Response: Let it run until it finishes (False)
# Code 4: The job was killed for exceed run time limits
### Response: Don't do anything (False)
# Code 5: The job ended with success
### Response: Collect ddG values (True)
def ddg_choice(ddg_response):
    '''Determine how to proceed with the ddG value collection based on the return code.'''
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
    '''Investigate if the cst_min application ran succesfully.'''
    # Glob to find error and out logs:
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


def write_ddG_results(folder, ddG_dict, ddG_dict_json, prot_name):
    '''Write all the results in ddG_dict to a file in the main folder for the given protein.'''
    ddG_dict_json = {prot_name: ddG_dict_json}
    resname = folder + '/ddG_results.txt'  # Name of the results file
    resname_json = folder + '/ddG_results.json'  # Name of the results file
    resname_pickle = folder + '/ddG_results.p'  # Name of the results file
    pickle.dump(ddG_dict_json, open(resname_pickle, "wb"))
    with open(resname_json, 'w') as outfile:
        json.dump(ddG_dict_json, outfile, sort_keys=True, indent=4)
    with open(resname, 'w') as fh_out:
        print('# ddG data for:', prot_name, file=fh_out)
        # Sort first according to the alphabet then according to length of the name.
        # Effectively printing the monomeric results first and the dimeric second.
        for name in sorted(sorted(ddG_dict), key=len):
            res_list = ddG_dict[name]
            if len(name) == 1:  # Monomeric stability
                print('# Monomeric stability for chain', name, file=fh_out)
                print('# ChainID Residue_number InsertionCode\tFromTo\tREU', file=fh_out)
                tab_list = ['\t'.join(t[0:]) for t in res_list]
                print('\n'.join(tab_list), file=fh_out)
            elif len(name) == 2:  # Dimer stability
                print('# Complex stability for chain complex', name, file=fh_out)
                print('# ChainID Residue_number InsertionCode\tFromTo\tREU', file=fh_out)
                tab_list = ['\t'.join(t[0:]) for t in res_list]
                print('\n'.join(tab_list), file=fh_out)


def check_mapping(mapping_dict):
    '''This is just a function to edit the mapping dict and remove the residue type.
    The reason for doing this is that the mapping between single letter codes in the
    ddg_predictions.txt are ambiguous with respect to the three letter code
    because of dextro -and modified residues.'''
    ### Shorten the keys:
    # From 'PRO A   1 '
    # To '   1 '
    new_mapping = dict()
    for k, v in mapping_dict.items():
        # Make sure that the residues are the same between the mappnings:
        if k[0:4] != v[0:4]:
            print(k, v)
            sys.exit()
        new_mapping[k[5:]] = v[4:]
    return(new_mapping)


def collect_ddg(ddg_file, mapping_dict, ddG_dict, ddG_dict_json):
    '''Collect the ddG results for a protein and write the results into a file,
     mapping the ddG values back to the original protein residue numbering.'''
    # Hacky way of getting the chain name from the cst_min pdb output name:
    name = ddg_file.split('/')[-1][12:-9]  # Chain name
    ddg_rundir = '/'.join(ddg_file.split('/')[0:-1]) + '/' + name + '_ddg_rundir'  # Chain rundir
    ddg_outfile = ddg_rundir + '/' + 'ddg_predictions.out'  # Prediction results from the ddg_monomer application
    # Get the original AA seq to test if any residues have been skipped by the ddg_monomer application:
    ori_AAseq = get_AA_string(ddg_rundir, name)
    # print(name, ori_AAseq)
    # print()
    # Correct the mapping dictionary to reflect usage of dextro -and modified residues:
    chain_map = check_mapping(mapping_dict[name])
    # Create a new list for the chain:
    ddG_dict[name] = list()
    ddG_dict_json[name] = dict()
    offset = 0     # Offset if ddg_monomer have skipped some residues
    res_count = 0  # Residue count according to the ddg_prediction.txt
    del_flag = 0   # Flag the chains where mapping fails
    with open(ddg_outfile) as fh:
        lines = fh.readlines()
        for line in lines:
            line = line.strip()
            if not line:
                continue
            elif line.startswith('ddG: description'):
                continue
            # Max number of skipped residues:
            if offset > 10 or (res_count - 1 + offset) > len(ori_AAseq):
                print("Offset is very high. Many residues have been skipped. Terminating:", ddg_outfile)
                del_flag = 1
                break
            # Extract results for the columns of the ddg_predictions.txt:
            cols = line.split()
            AA = cols[1][0]
            AA_mut = cols[1][-1]
            res_count = int(cols[1][1:-1])
            # Check for skipped residue:
            ori_AA = ori_AAseq[res_count - 1 + offset]
            while ori_AA != AA:  # Start looping if the chain is broken
                print(res_count)
                print('obs:', AA)
                print('ori', ori_AAseq[res_count - 1 + offset])
                print()
                offset += 1
                # Break either if the offset gets to high or the matching amino acid is found:
                if offset > 10 or (res_count - 1 + offset) >= len(ori_AAseq):
                    del_flag = 1
                    break
                elif ori_AAseq[res_count - 1 + offset] == AA:
                    break
            if del_flag:
                break
            # Create the new residue mapping:
            res_idx_string = '{:>4} '.format(res_count + offset)  # Notice the blank insertion code
            # Look-up the original residue mapping:
            ori_res_idx_string = chain_map[res_idx_string]
            # Create a "from to" string:
            from_to = AA + '->' + AA_mut
            # Put this into a tuple with the ddG value and append it the results dictionary:
            res_tuple = (ori_res_idx_string, from_to, cols[2])
            ddG_dict[name].append(res_tuple)

            chain_name = ori_res_idx_string[0]
            resn = AA + ori_res_idx_string[1:].strip()
            if chain_name not in ddG_dict_json[name]:
                ddG_dict_json[name][chain_name] = dict()
            if resn not in ddG_dict_json[name][chain_name]:
                ddG_dict_json[name][chain_name][resn] = dict()
            # Special case for e.g. the histidine tautomer:
            if AA_mut in ddG_dict_json[name][chain_name][resn] and\
               ddG_dict_json[name][chain_name][resn][AA_mut] < float(cols[2]):
                pass
            else:
                ddG_dict_json[name][chain_name][resn][AA_mut] = float(cols[2])

        # If faulty run, delete the chain:
        if del_flag:
            del ddG_dict[name]
            del ddG_dict_json[name]
    return(ddG_dict, ddG_dict_json)


def get_AA_string(ddg_rundir, name):
    '''Get the original amino acid sequence for the given chain. Monomer or complex.'''
    fnam = ddg_rundir + '/' + name  # Path of the renamed chain we are looking at
    with open(fnam) as fh:
        lines = fh.readlines()
        prev_resnumb = '   0'
        AA_string = ''
        for line in lines:
            if not line.startswith('ATOM'):
                continue
            # When observing a new resnumber add the amino acid:
            resnumb = line[22:27]
            if resnumb != prev_resnumb:
                prev_resnumb = resnumb
                res_name = line[17:20]
                AA = residue_type_3to1_map[res_name]
                AA_string += AA
    return(AA_string)


def median(lst):
    return(np.median(np.array(lst)))


def collect_REU(ddg_file, ddG_dict, ddG_dict_json):
    name = ddg_file.split('/')[-1][12:-9]  # Chain name
    if name not in ddG_dict or name not in ddG_dict_json:
        return(ddG_dict, ddG_dict_json)

    ddg_rundir = '/'.join(ddg_file.split('/')[0:-1]) + '/' + name + '_ddg_rundir'  # Chain rundir
    # Glob to find the output log:
    ddg_log_glob = ddg_rundir + '/' + '*_log.out'
    ddg_logfile = glob.glob(ddg_log_glob)
    # There can only be one log:
    assert(len(ddg_logfile) == 1)
    ddg_logfile = ddg_logfile[0]
    energy_list = list()
    # Now read the log file:
    with open(ddg_logfile) as fh:
        lines = fh.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith('Total weighted score:'):
                REU = line.split()[-1]
                REU = float(REU)
                energy_list.append(REU)

    median_REU = statistics.median_grouped(energy_list)
    res_tuple = ('Total_energy', 'X->X', str(median_REU))
    ddG_dict[name].append(res_tuple)
    ddG_dict_json[name]['energy_total'] = median_REU
    return(ddG_dict, ddG_dict_json)


def run_parallel(folder):
    # If the cst_min failed then return:
    if not cst_min_success(folder):
        return(False)

    prot_name = folder.split('/')[-1]
    # Use the residue mapping created by the run_cst app:
    mapping_path = folder + '/' + 'residue_mapping_dict.p'
    mapping_dict = pickle.load(open(mapping_path, "rb"))
    # Glob to find all the cst_min minimized structures:
    ddg_file_glob = folder + '/' + '*.pdb'
    files_for_ddg = glob.glob(ddg_file_glob)
    # If previous runs have already generated results delete them:
    resname = folder + '/ddG_results.txt'
    if os.path.isfile(resname):
        os.remove(resname)
    # Run through all the monomers/dimers and fill up the ddG_dict:
    ddG_dict = dict()
    ddG_dict_json = dict()
    for ddg_file in files_for_ddg:
        # Check if the ddG run was succesfull and if not skip this chain:
        ddg_response = ddg_success(ddg_file)
        if ddg_choice(ddg_response) is False:
            continue
        # Add the results to the ddG_dict:
        ddG_dict, ddG_dict_json = collect_ddg(ddg_file, mapping_dict, ddG_dict, ddG_dict_json)
        # Add the total REU of each protein monomer and dimer:
        ddG_dict, ddG_dict_json = collect_REU(ddg_file, ddG_dict, ddG_dict_json)
    # Finally write the result to a file in the protein folder:
    write_ddG_results(folder, ddG_dict, ddG_dict_json, prot_name)


def merge_results(folder_list):
    master_dict = dict()
    for folder in folder_list:
        resname_pickle = folder + '/ddG_results.p'  # Name of the results file
        if not os.path.isfile(resname_pickle):
            continue
        ddG_dict_json = pickle.load(open(resname_pickle, "rb"))
        master_dict.update(ddG_dict_json)
    all_resfile = run_dir + '/all_results.json'
    with open(all_resfile, 'w') as outfile:
        json.dump(master_dict, outfile, sort_keys=True, indent=4)


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

    # Paralellize the process:
    pool = multiprocessing.Pool(28)
    output = pool.map(run_parallel, folder_list)
    merge_results(folder_list)
