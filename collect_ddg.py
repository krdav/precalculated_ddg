
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

db_home_dir = '/home/projects/cu_10020/data/precalculated_ddg'
db_split_dir = db_home_dir + '/split'
prot_list_file = '/home/projects/cu_10020/data/precalculated_ddg/prot_list.txt'

rosetta_db = '/services/tools/rosetta/2016.10/main/database'
rosetta_ddg_app = '/services/tools/rosetta/2016.10/main/source/bin/ddg_monomer.default.linuxgccrelease'

const_flags_ddg = '-resfile resfile.txt -ddg:weight_file soft_rep_design -ddg::iterations 1 -ddg::dump_pdbs false -ignore_unrecognized_res -ddg::local_opt_only true -ddg::suppress_checkpointing true -in::file::fullatom -ddg::mean true -ddg::min false -mute all -ddg::output_silent false'
# -in:file:s min_cst_0.5.XXXX_0001.pdb



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





def collect_ddg(pdb_file_path, mapping_dict):
    # Hacky way of getting the chain/pair name:
    name = pdb_file_path.split('/')[-1][12:-9]
    ddg_rundir = '/'.join(pdb_file_path.split('/')[0:-1]) + '/' + name + '_ddg_rundir'
    os.mkdir(ddg_rundir)
    dst_pdb_file_path = ddg_rundir + '/' + name
    # Copy the relevant pdb file:
    shutil.copy(pdb_file_path, dst_pdb_file_path)
    # Copy the resfile to the ddg_rundir:
    resfile_path = '/'.join(pdb_file_path.split('/')[0:-1]) + '/' + 'resfile.txt'
    dst_resfile_path = ddg_rundir + '/' + 'resfile.txt'
    shutil.copy(resfile_path, dst_resfile_path)


def dgg_success(folder):
    out_log_glob = folder + '/' + 'sub*_log.out'
    if len(out_log_glob) > 1:
        print('This is weird... There are more than submission STDOUT log. Maybe clean folder first? Folder:\n', folder)
        return(False)
    else:
        out_log = out_log_glob[0]
    with open(out_log) as fh:
        lines = fh.readlines()
        # Apparently the cst_min app ends with this on success:
    if lines[-1].startswith('running another iteration of minimization'):
        return(True)
    else:
        print('The cst_min application did not end succesfully in folder:\n', folder)



folders_for_update = db_home_dir + '/' + 'folders_for_update.txt'

np = 1

with open(folders_for_update) as fh:
    folder_list = fh.read().splitlines()

for idx, folder in enumerate(folder_list):
    mapping_path = folder + '/' + 'residue_mapping_dict.p'
    mapping_dict = pickle.load(open(mapping_path, "rb"))
#    if not dgg_success(folder):
#        continue
    ddg_file_glob = folder + '/' + '*.pdb'
    files_for_ddg = glob.glob(ddg_file_glob)
    for ddg_file in files_for_ddg:
        collect_ddg(ddg_file, mapping_dict)












