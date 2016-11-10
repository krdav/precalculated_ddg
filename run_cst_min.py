import sys
import os
import shutil
import gzip
import pickle
import time
import argparse
import glob

sys.path.insert(0, '/home/projects/pr_46690/apps/python3-site-packages/lib/python/')
from Bio.PDB import *


# Build commandline parser
parser = argparse.ArgumentParser(description="Run the carbon alpha constraint minimization in Rosetta.")

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
    "-prot_list_file",
    "--prot_list_file",
    type=str,
    dest="prot_list_file",
    metavar="FILE",
    help="File with a list of paths for input to the database.",
    required=True)
parser.add_argument(
    "-rosetta_db",
    "--rosetta_db",
    type=str,
    dest="rosetta_db",
    metavar="PATH",
    help="Absolute path to the Rosetta database.",
    required=False)
parser.add_argument(
    "-rosetta_min_cst_app",
    "--rosetta_min_cst_app",
    type=str,
    dest="rosetta_min_cst_app",
    metavar="PATH",
    help="Absolute path to the Rosetta cst_min application.",
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
parser.add_argument(
    "-skip_RW",
    "--skip_RW",
    type=int,
    dest="skip_RW",
    metavar="SWITCH",
    help="Skip the read/write step and go directly to the \"folders_for_update.txt\" run. Default 0.",
    required=False)

# Set default arguments
parser.set_defaults(
    rosetta_db='/services/tools/rosetta/2016.10/main/database',
    rosetta_min_cst_app='/services/tools/rosetta/2016.10/main/source/bin/minimize_with_cst.default.linuxgccrelease',
    verbose_error=1,
    restart_failed=0,
    skip_RW=0)

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
if args.prot_list_file[0] != '/':
    args.prot_list_file = run_dir + '/' + args.prot_list_file


# Global variables:
const_flags_min_cst = '-in:file:fullatom -ignore_unrecognized_res -fa_max_dis 9.0 -ddg::harmonic_ca_tether 0.5 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst_0.5 -ddg::sc_min_only false'
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


def QC_prot_check(prot_path):
    '''Check the quality of the protein in the list.'''
    min_res = 10    # Minimum resolution
    resolution = 0  # Resolution good/bad switch
    xray = 0        # Method good/bad switch
    canonical = 1   # Residue good/bad switch
    with gzip.open(prot_path) as fh:
        lines = fh.readlines()
        for line in lines:
            # Decode from binary to unicode:
            line = str(line, 'utf-8')
            # Find the line that tells if the PBD file is from X-RAY diffraction:
            if line.startswith('EXPDTA    X-RAY DIFFRACTION'):
                xray = 1
            # Find the resolution line:
            elif line.startswith('REMARK   2 RESOLUTION.'):
                # Validate that the resolution is under the cutoff:
                res_str = line[22:29]
                try:
                    res_float = float(res_str)
                    if res_float < min_res:
                        resolution = 1
                except:
                    pass
            # Exlude amino acids on the blacklist:
            elif line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                canonical = 0
                break
        # A PDB ID is kept if both xray and resolution flags checks out:
        keep = xray * resolution * canonical  # All switches must be on to keep
    return(keep)


class ReadWriteProtein:
    '''Class to read and write the input PDB file.'''
    def _fetch_chain_lines(self, prot_path):
        '''Read the PDB file into a dictionary based on each chain having its own key.'''
        chain_lines = dict()
        with gzip.open(prot_path) as fh:
            lines = fh.readlines()
            for line in lines:
                # Decode from binary to unicode:
                line = str(line, 'utf-8')
                line = line.rstrip('\n')
                # Find residues and ligands:
                if (line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM') and line[17:20] in whitelist:  # Require residues to be known by Rosetta
                    chain_name = line[21]
                    # Make the dictionary by chain name:
                    if chain_name not in chain_lines:
                        chain_lines[chain_name] = list()
                        chain_lines[chain_name].append(line)
                    else:
                        chain_lines[chain_name].append(line)
        return(chain_lines)

    def _make_split_key(self, prot_path):
        '''Function to make the two character split key from either a normal PDB filename (easy case) or a homology model filename.'''
        path_list = prot_path.split('/')
        if 'pdb' in path_list[-1]:
            split_key = path_list[-1][4:6]
        else:
            pass  # DO SOMETHING HERE IF HOMOLOGY MODELS ARE TO BE USED
        return(split_key)

    def _make_prot_name(self, prot_path):
        '''Make a name for the protein.'''
        path_list = prot_path.split('/')
        if 'pdb' in path_list[-1]:
            prot_name = path_list[-1][3:7]
        else:
            pass  # DO SOMETHING HERE IF HOMOLOGY MODELS ARE TO BE USED
        return(prot_name)

    def _write_protein(self, new_prot_folder, prot_name, chain_lines):
        '''Copy the protein as it is from the provided compressed file.'''
        outname = new_prot_folder + '/' + prot_name
        with open(outname, 'w') as fh_out:
            for chain_name, lines in chain_lines.items():
                print('\n'.join(lines), file=fh_out)

    def __init__(self, prot_path, db_split_dir):
        self.prot_path = prot_path
        self.chain_lines = self._fetch_chain_lines(self.prot_path)
        self.split_key = self._make_split_key(self.prot_path)
        self.prot_name = self._make_prot_name(self.prot_path)
        self.split_folder = db_split_dir + '/' + self.split_key
        self.mapping_dict = dict()
        # Make the split folder if not present already:
        if not os.path.exists(self.split_folder):
            os.mkdir(self.split_folder)
        # Make the protein folder in not already present:
        self.new_prot_folder = self.split_folder + '/' + self.prot_name
        if not os.path.exists(self.new_prot_folder):
            os.mkdir(self.new_prot_folder)
        # Write the native protein to the folder:
        self._write_protein(self.new_prot_folder, self.prot_name, self.chain_lines)

    def write_monomers(self):
        '''Write a monomeric protein chain to a file. Renumber the residue so they have consecutive numbering.
        Store the mapping so the original residue number can be recovered later.'''
        filenames = list()
        # Run through all the chains:
        for chain_name, lines in self.chain_lines.items():
            outname = self.new_prot_folder + '/' + chain_name
            filenames.append(outname)
            self.mapping_dict[chain_name] = dict()
            with open(outname, 'w') as fh_out:
                atom_idx = 0
                res_idx = 0
                prev_resnumb = '   X '
                for line in lines:
                    # Convert line to list:
                    lline = list(line)
                    atom_idx += 1
                    resnumb = line[22:27]
                    # The new residue number is incremented upon finding a new residue number in the original file:
                    if resnumb != prev_resnumb:
                        prev_resnumb = resnumb
                        res_idx += 1
                    lline[6:11] = list('{:>5}'.format(atom_idx))  # New atom index
                    res_idx_string = '{:>4} '.format(res_idx)     # New residue index, notice the blank insertion code
                    lline[22:27] = list(res_idx_string)
                    # Recreate the renumbered line and print it:
                    new_line = ''.join(lline)
                    print(new_line, file=fh_out)
                    # Save the mapping:
                    self.mapping_dict[chain_name][new_line[17:27]] = line[17:27]  # Notice that the list slice is including both residue type, chain name, residue number and insertion code
        # Dump the mapping:
        mapping_path = self.new_prot_folder + '/' + 'residue_mapping_dict.p'
        pickle.dump(self.mapping_dict, open(mapping_path, "wb"))
        return(filenames)

    def _find_dimers(self, new_prot_folder, prot_name, chain_lines):
        '''Look through all atoms in a PDB file to find the chains that have atoms closer than within a cutoff distance.'''
        parser = PDBParser()  # Use the Biopython PDB module to find the distances
        new_prot_path = new_prot_folder + '/' + prot_name
        pdb_obj = parser.get_structure(prot_name, new_prot_path)[0]  # Notice the first model in the structures object is fetch right away
        dist_cut = 20        # Interactions cannot be further away than 20A backbone atom to backbone atom
        # R -> E = 7.1A + 4.9A = 12A + bond max 15A to give a bit of freedom go to 20A
        interaction_cut = 3  # Maximum distance that defines an interaction
        dimers = list()
        chains = list(chain_lines.keys())
        chains_copy = chains[:]
        for chain1 in chains:
            chain1_obj = pdb_obj[chain1]
            chains_copy.pop(0)  # Remove one each time to avoid self -and double comparison
            for chain2 in chains_copy:
                chain2_obj = pdb_obj[chain2]
                stop_comparing = 0
                for res1 in chain1_obj:
                    if stop_comparing: break  # Break early if an interaction have been found
                    for res2 in chain2_obj:
                        if stop_comparing: break  # Break early if an interaction have been found
                        # Find a backbone atom that can be used to find the distance.
                        # CA is prefered:
                        if 'CA' in res1 and 'CA' in res2:
                            res_comp = 'CA'
                        elif 'N' in res1 and 'N' in res2:
                            res_comp = 'N'
                        elif 'C' in res1 and 'C' in res2:
                            res_comp = 'C'
                        # If no backbone atoms are found no interactions are annotated:
                        else:
                            continue

                        idist = res1[res_comp] - res2[res_comp]
                        # Only when the backbone distance is smaller than the cutoff,
                        # continue and find the distance between all atoms:
                        if idist < dist_cut:
                            for atom1 in res1:      # All atoms for the residue in the chain for comparison
                                if stop_comparing: break  # Break early if an interaction have been found
                                for atom2 in res2:  # All atoms for the residue in the chosen chain
                                    idist = atom1 - atom2
                                    # Is this distance short enough to imply an interation?:
                                    if idist <= interaction_cut:
                                        dimers.append(chain1 + chain2)
                                        stop_comparing = 1
                                        break
        # Sort pairs before returning:
        dimers = [''.join(sorted(d)) for d in dimers]
        return(dimers)

    def write_dimers(self):
        '''Write a dimeric protein chain to a file. Renumber the residue so they have consecutive numbering.
        Store the mapping so the original residue number can be recovered later.'''
        # Find the dimers based on a distance criteria:
        self.dimers = self._find_dimers(self.new_prot_folder, self.prot_name, self.chain_lines)
        filenames = list()  # List of chain filenames
        # The pair must be two characters chain_name1chain_name2
        # and the order of the characters must be alphabetic!
        for pair in self.dimers:
            self.mapping_dict[pair] = dict()
            outname = self.new_prot_folder + '/' + pair
            filenames.append(outname)
            # Fetch all the lines from both chains:
            lines1 = self.chain_lines[pair[0]]
            lines2 = self.chain_lines[pair[1]]
            all_lines = lines1 + lines2
            # Now write the dimeric PDB file:
            with open(outname, 'w') as fh_out:
                atom_idx = 0
                res_idx = 0
                prev_resnumb = '   X '
                new_chain_name = 'A'  # Rename the new concatenated chain, A
                for line in all_lines:
                    # Convert line to list:
                    lline = list(line)
                    atom_idx += 1
                    resnumb = line[22:27]
                    # The new residue number is incremented upon finding a new residue number in the original file:
                    if resnumb != prev_resnumb:
                        prev_resnumb = resnumb
                        res_idx += 1
                    lline[6:11] = list('{:>5}'.format(atom_idx))  # New atom index
                    lline[21] = new_chain_name
                    res_idx_string = '{:>4} '.format(res_idx)     # New residue index, notice the blank insertion code
                    lline[22:27] = list(res_idx_string)
                    # Recreate the renumbered line and print it:
                    new_line = ''.join(lline)
                    print(new_line, file=fh_out)
                    # Save the mapping:
                    self.mapping_dict[pair][new_line[17:27]] = line[17:27]  # Notice that the list slice is including both residue type, chain name, residue number and insertion code

        # Dump the mapping:
        mapping_path = self.new_prot_folder + '/' + 'residue_mapping_dict.p'
        pickle.dump(self.mapping_dict, open(mapping_path, "wb"))
        return(filenames)


def clean_list(l):
    '''Removes duplicates and empty string values while keeping the order of the list.'''
    seen = set()
    seen_add = seen.add
    seen_add('')
    return([x for x in l if not (x in seen or seen_add(x))])


def pbs_submit_cmd(np, cmd_flags, run_dir, idx):
    '''Submit the cst_min minimization as a job to the queueing system.'''
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


def make_split_key(prot_path):
    '''Function to make the two character split key from either a normal PDB filename (easy case) or a homology model filename.'''
    path_list = prot_path.split('/')
    if 'pdb' in path_list[-1]:
        split_key = path_list[-1][4:6]
    else:
        pass  # DO SOMETHING HERE IF HOMOLOGY MODELS ARE TO BE USED
    return(split_key)


def make_prot_name(prot_path):
    '''Make the protein name from the original path of the compressed input PDB.'''
    path_list = prot_path.split('/')
    if 'pdb' in path_list[-1]:
        prot_name = path_list[-1][3:7]
    else:
        pass  # DO SOMETHING HERE IF HOMOLOGY MODELS ARE TO BE USED
    return(prot_name)


def write_resfile(new_prot_folder):
    '''Write a resfile specifying that all residue are repackable.'''
    s = 'ALLAA'
    outname = new_prot_folder + '/resfile.txt'
    with open(outname, 'w') as fh_out:
        print(s, file=fh_out)


def write_update_folder(db_home_dir, update_sack):
    '''Empty the sack of paths to proteins to run.'''
    outname = db_home_dir + '/' + 'folders_for_update.txt'
    with open(outname, 'w') as fh_out:
        print('\n'.join(update_sack), file=fh_out)


def cst_min_success(prot_path, db_split_dir):
    '''Investigate whether the cst_min run was a success or not.'''
    split_key = make_split_key(prot_path)
    prot_name = make_prot_name(prot_path)
    new_prot_folder = db_split_dir + '/' + split_key + '/' + prot_name
    # If the folder does not exist return code 1:
    if not os.path.exists(new_prot_folder):
        return(1)
    else:
        pass
    # Glob to find the error -and out logs:
    out_log_glob_string = new_prot_folder + '/' + 'sub*_log.out'
    out_log_glob = glob.glob(out_log_glob_string)
    err_log_glob_string = new_prot_folder + '/' + 'sub*_log.err'
    err_log_glob = glob.glob(err_log_glob_string)
    if len(out_log_glob) > 1 and len(err_log_glob) > 1:  # Code 2: Unforseen error
        if args.verbose or args.verbose_error:
            print('This is weird... There are more than one ddG submission STDOUT log. Maybe clean folder first? Folder:\n', new_prot_folder)
        return(2)
    elif len(out_log_glob) == 0 and len(err_log_glob) == 0:  # Code 3: The job is still running
        if args.verbose:
            print('No log found, job must be in progress:\n', new_prot_folder)
        # Job in progress log not yet created:
        return(3)
    else:  # Check if the logs indicate success
        out_log = out_log_glob[0]
        err_log = err_log_glob[0]
    # Open the error log and parse it:
    with open(err_log) as fh:
        lines = fh.readlines()
    if not lines:  # If the error log is empty
        pass
    elif lines[-1].startswith('=>> PBS: job killed'):  # Code 4: The job was killed for exceed run time limits
        if args.verbose or args.verbose_error:
            print('Job was killed by the queueing system because of too much run time:\n', err_log)
        return(4)
    elif 'Aborted' in lines[-1]:  # Code 2: Unforseen error
        if args.verbose or args.verbose_error:
            print('Rosetta have thrown an error an aborted:\n', new_prot_folder)
        return(2)
    # Open the output log and parse it:
    with open(out_log) as fh:
        lines = fh.readlines()
        # Apparently the cst:min app ends with this on success:
    if lines[-1].startswith('running another iteration of minimization') or lines[-1].startswith('core.optimization.LineMinimizer: Inaccurate G!'):  # Code 5: The job ended with success
        if args.verbose:
            print('Job has succesfully ended:\n', new_prot_folder)
        return(5)
    elif 'already exists! skipping' in lines[-1]:
        if args.verbose:
            print('Job was submitted multiple times:\n', new_prot_folder)
        return(6)
    else:  # Code 2: Unforseen error
        if args.verbose or args.verbose_error:
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
### Response: Do nothing but inform the user (False)
# Code 5: The job ended with success
### Response: Nothing (False)
# Code 6: The job was submitted multiple times resulting in overlapping output names
### Response: Nothing (False)
def min_cst_choice(response, db_split_dir, prot_path):
    '''Determine what the response should be based on the cst_run return code.'''
    split_key = make_split_key(prot_path)
    prot_name = make_prot_name(prot_path)
    new_prot_folder = db_split_dir + '/' + split_key + '/' + prot_name

    if response == 1:
        return(True)
    elif response == 2 and args.restart_failed:  # Some kind of fail, try to restart
        shutil.rmtree(new_prot_folder)
        return(True)
    elif response == 2 and not args.restart_failed:  # Some kind of fail
        return(False)
    elif response == 3:  # Still running
        return(False)
    elif response == 4:  # Kill because of too much time
        return(False)
    elif response == 5:  # Success
        return(False)
    elif response == 6:  # Probably success
        return(False)
    else:
        if args.verbose or args.verbose_error:
            print('Unforseen error. ddg_success could not be determined.')
        return(False)


if __name__ == "__main__":
    # Skip reading and writing of the protein,
    # and go directly to submission of cst_min:
    if not args.skip_RW:
        # Read the input list:
        with open(args.prot_list_file) as fh:
            prot_list = fh.read().splitlines()
        prot_list = clean_list(prot_list)
        update_sack = list()  # Sack of paths to proteins that requires runnning cst_min
        for prot_path in prot_list:
            # Skip bad proteins:
            if not QC_prot_check(prot_path):
                print('Following protein did not meet the QC requiemnents:\n{}'.format(prot_path))
                continue
            # Skip runs that are already running or failing for some reason:
            cst_min_response = cst_min_success(prot_path, args.db_split_dir)
            cst_process = min_cst_choice(cst_min_response, args.db_split_dir, prot_path)
            if not cst_process:
                continue
            # Read protein into an object:
            RWprotein_obj = ReadWriteProtein(prot_path, args.db_split_dir)
            # Get info from the read/write protein object:
            split_key = RWprotein_obj.split_key
            prot_name = RWprotein_obj.prot_name
            new_prot_folder = RWprotein_obj.new_prot_folder
            monomer_filenanes = RWprotein_obj.write_monomers()
            dimer_filenanes = RWprotein_obj.write_dimers()
            # Concatenate all the monomer/dimer filename together and write them to the protein folder.
            # This list will be used as input to the cst_min later.
            all_filenames = monomer_filenanes + dimer_filenanes
            filenames_outname = new_prot_folder + '/' + 'filenames_for_cst_min.txt'
            with open(filenames_outname, 'w') as fh_out:
                print('\n'.join(all_filenames), file=fh_out)
            # Write a resfile to the folder with the newly created PDB files:
            write_resfile(new_prot_folder)
            # Put the folder in the sack for cst_min and ddG calculations:
            update_sack.append(new_prot_folder)
        # Then write the folder of the files that needs updating:
        write_update_folder(args.db_home_dir, update_sack)

    # Read the above writen, or pregenerated, list of proteins to update:
    folders_for_update = args.db_home_dir + '/' + 'folders_for_update.txt'
    cst_filelist_name = 'filenames_for_cst_min.txt'
    with open(folders_for_update) as fh:
        folder_list = fh.read().splitlines()
    # Then loop through all protein folders, submitting them each as a job:
    for idx, folder in enumerate(folder_list):
        if folder == '':
            continue
        cst_filelist_path = folder + '/' + cst_filelist_name
        cst_cmd = args.rosetta_min_cst_app + ' ' + const_flags_min_cst + ' -database ' + args.rosetta_db + ' -l ' + cst_filelist_path
        print('Submitting for:', cst_filelist_path)
        try:
            pbs_submit_cmd(np, cst_cmd, folder, idx)
        except Exception as e:
            print(e)
