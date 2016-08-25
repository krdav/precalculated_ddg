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

db_home_dir = '/home/projects/cu_10020/data/precalculated_ddg'
db_split_dir = db_home_dir + '/split'
# prot_list_file = '/home/projects/cu_10020/data/precalculated_ddg/top1003_constraint_pdbdir.txt'
prot_list_file = '/home/projects/cu_10020/data/precalculated_ddg/prot_list.txt'

rosetta_db = '/services/tools/rosetta/2016.10/main/database'
rosetta_min_cst_app = '/services/tools/rosetta/2016.10/main/source/bin/minimize_with_cst.default.linuxgccrelease'

const_flags_min_cst = '-in:file:fullatom -ignore_unrecognized_res -fa_max_dis 9.0 -ddg::harmonic_ca_tether 0.5 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst_0.5 -ddg::sc_min_only false'
# -l file_list.txt

run_dir = os.getcwd()

os.chdir(db_home_dir)




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
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "UNK": 'X',
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


def QC_prot_check(prot_path):
    res = 10  # Minimum resolution
    resolution = 0
    xray = 0
    canonical = 1
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
                    if res_float < res:
                        resolution = 1
                except:
                    pass
            # Exlude amino acids on the blacklist:
            elif (line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM') and line[17:20] in blacklist:
                canonical = 0
                break
        # A PDB ID is kept if both xray and resolution flags checks out:
        keep = xray * resolution * canonical  # All switches must be on to keep
    return(keep)


class ReadWriteProtein:
    def _fetch_chain_lines(self, prot_path):
        chain_lines = dict()
        with gzip.open(prot_path) as fh:
            lines = fh.readlines()
            for line in lines:
                # Decode from binary to unicode:
                line = str(line, 'utf-8')
                line = line.rstrip('\n')
                # Find residues and ligands:
                if (line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM') and line[17:20] not in greylist:
                    chain_name = line[21]
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
            pass  # Do something here later
        return(split_key)

    def _make_prot_name(self, prot_path):
        '''gggg'''
        path_list = prot_path.split('/')
        if 'pdb' in path_list[-1]:
            prot_name = path_list[-1][3:7]
        else:
            pass  # Do something here later
        return(prot_name)

    def _write_protein(self, new_prot_folder, prot_name, chain_lines):
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
        filenames = list()
        for chain_name, lines in self.chain_lines.items():
            outname = self.new_prot_folder + '/' + chain_name
            filenames.append(outname)
            self.mapping_dict[chain_name] = dict()
            with open(outname, 'w') as fh_out:
                atom_idx = 0
                res_idx = 0
                cur_res = 'XXX'
                for line in lines:
                    lline = list(line)
                    atom_idx += 1
                    if line[17:20] != cur_res:
                        cur_res = line[17:20]
                        res_idx += 1
                    lline[6:11] = list('{:>5}'.format(atom_idx))
                    res_idx_string = '{:>4} '.format(res_idx)  # Notice the blank insertion code
                    lline[22:27] = list(res_idx_string)
                    # Save the mapping:
                    # self.mapping_dict[chain_name][res_idx_string] = line[21:27]  # Notice that the list slice is including both chain name, residue number and insertion code
                    new_line = ''.join(lline)
                    print(new_line, file=fh_out)
                    # Save the mapping:
                    self.mapping_dict[chain_name][new_line[17:27]] = line[17:27]  # Notice that the list slice is including both residue type, chain name, residue number and insertion code

        # Dump the mapping:
        mapping_path = self.new_prot_folder + '/' + 'residue_mapping_dict.p'
        pickle.dump(self.mapping_dict, open(mapping_path, "wb"))
        return(filenames)

    def _find_dimers(self, new_prot_folder, prot_name, chain_lines):
        parser = PDBParser()
        new_prot_path = new_prot_folder + '/' + prot_name
        pdb_obj = parser.get_structure(prot_name, new_prot_path)[0]  # Notice the first model in the structures object is fetch right away
        dist_cut = 20  # Interactions cannot be further away than 20A backbone atom to backbone atom
        interaction_cut = 3  # Maximum distance that defines an interaction
        # R -> E = 7.1A + 4.9A = 12A + bond max 15A to give a bit of freedom go to 20A
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
        self.dimers = self._find_dimers(self.new_prot_folder, self.prot_name, self.chain_lines)
        filenames = list()
        # The pair must be two characters chain_name1chain_name2
        # and the order of the characters must be alphabetic!
        for pair in self.dimers:
            self.mapping_dict[pair] = dict()
            outname = self.new_prot_folder + '/' + pair
            filenames.append(outname)
            lines1 = self.chain_lines[pair[0]]
            lines2 = self.chain_lines[pair[1]]
            all_lines = lines1 + lines2

            with open(outname, 'w') as fh_out:
                atom_idx = 0
                res_idx = 0
                cur_res = 'XXX'
                new_chain_name = 'A'
                for line in all_lines:
                    lline = list(line)
                    atom_idx += 1
                    if line[17:20] != cur_res:
                        cur_res = line[17:20]
                        res_idx += 1
                    lline[6:11] = list('{:>5}'.format(atom_idx))
                    lline[21] = new_chain_name
                    res_idx_string = '{:>4} '.format(res_idx)  # Notice the blank insertion code
                    lline[22:27] = list(res_idx_string)
                    # Save the mapping:
                    # self.mapping_dict[pair][res_idx_string] = line[21:27]  # Notice that the list slice is including both chain name, residue number and insertion code
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
    os.chdir(run_dir)
    log_err = run_dir + '/sub' + str(idx) + '_log.err'
    log_out = run_dir + '/sub' + str(idx) + '_log.out'

    qsub_string = '#!/bin/sh\n\
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
#PBS -l nodes=1:ppn=' + str(np) + ':thinnode\n\
#PBS -l walltime=24:00:00\n\
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
    time.sleep(1)


def make_split_key(prot_path):
    '''Function to make the two character split key from either a normal PDB filename (easy case) or a homology model filename.'''
    path_list = prot_path.split('/')
    if 'pdb' in path_list[-1]:
        split_key = path_list[-1][4:6]
    else:
        pass  # Do something here later
    return(split_key)


def make_prot_name(prot_path):
    '''gggg'''
    path_list = prot_path.split('/')
    if 'pdb' in path_list[-1]:
        prot_name = path_list[-1][3:7]
    else:
        pass  # Do something here later
    return(prot_name)


def entry_exists(prot_path, db_split_dir):
    split_key = make_split_key(prot_path)
    prot_name = make_prot_name(prot_path)
    new_prot_folder = db_split_dir + '/' + split_key + '/' + prot_name
    return(os.path.exists(new_prot_folder))


def write_resfile(new_prot_folder):
    s = 'ALLAA'
    outname = new_prot_folder + '/resfile.txt'
    with open(outname, 'w') as fh_out:
        print(s, file=fh_out)


def write_update_folder(db_home_dir, update_sack):
    outname = db_home_dir + '/' + 'folders_for_update.txt'
    with open(outname, 'w') as fh_out:
        print('\n'.join(update_sack), file=fh_out)


if __name__ == "__main__":
    with open(prot_list_file) as fh:
        prot_list = fh.read().splitlines()
    prot_list = clean_list(prot_list)
    update_sack = list()
    # os.chdir(db_split_dir)
    for prot_path in prot_list:
        if entry_exists(prot_path, db_split_dir):
            continue
        # Skip bad proteins:
        if not QC_prot_check(prot_path):
            print('Following protein did not meet the QC requiemnents:\n{}'.format(prot_path))
            continue
        # Read protein into an object:
        RWprotein_obj = ReadWriteProtein(prot_path, db_split_dir)

        split_key = RWprotein_obj.split_key
        prot_name = RWprotein_obj.prot_name
        new_prot_folder = RWprotein_obj.new_prot_folder

        monomer_filenanes = RWprotein_obj.write_monomers()
        dimer_filenanes = RWprotein_obj.write_dimers()
        all_filenames = monomer_filenanes + dimer_filenanes
        filenames_outname = new_prot_folder + '/' + 'filenames_for_cst_min.txt'
        with open(filenames_outname, 'w') as fh_out:
            print('\n'.join(all_filenames), file=fh_out)

        # Write a resfile to the folder with the newly created PDB files:
        write_resfile(new_prot_folder)
        # Put the folder in the sack for cst_min and ddG calculations:
        update_sack.append(new_prot_folder)
    # Then write the folder of the files that needs updating:
    write_update_folder(db_home_dir, update_sack)

    folders_for_update = db_home_dir + '/' + 'folders_for_update.txt'
    cst_filelist_name = 'filenames_for_cst_min.txt'
    np = 1

    with open(folders_for_update) as fh:
        folder_list = fh.read().splitlines()

    for idx, folder in enumerate(folder_list):
        cst_filelist_path = folder + '/' + cst_filelist_name
        cst_cmd = rosetta_min_cst_app + ' ' + const_flags_min_cst + ' -database ' + rosetta_db + ' -l ' + cst_filelist_path
        print('Submitting for:', cst_filelist_path)
        pbs_submit_cmd(np, cst_cmd, folder, idx)







