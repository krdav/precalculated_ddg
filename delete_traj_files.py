import glob
import multiprocessing
import os
import sys
import argparse


# Build commandline parser
parser = argparse.ArgumentParser(description="Remove the traj files created by the Rosetta ddg_monomer app.")

# Arguments
parser.add_argument(
    "-db_split_dir",
    "--db_split_dir",
    type=str,
    dest="db_split_dir",
    metavar="PATH",
    help="Absolute path to the split directory containing all the two character folders, pdb style.",
    required=True)

# Set default arguments
parser.set_defaults()

# Put arguments into args object:
args = parser.parse_args()


def delete_traj(split):
    '''Find all protein directories in the split directory. Walk thorugh them and delete all _traj files.'''
    # Find the protein folders:
    prot_glob = split + '*/'
    prots = glob.glob(prot_glob)
    for prot in prots:
        # Find all the chain folders:
        chain_glob = prot + '*/'
        chains = glob.glob(chain_glob)
        for chain in chains:
            # Then find the all _traj files
            traj_glob = chain + '*_traj*'
            trajs = glob.glob(traj_glob)
            for traj in trajs:
                # Delete the _traj files one by one.
                # Use try to avoid termination if moving/removing while the script is running.
                try:
                    os.remove(traj)
                except:
                    pass


if __name__ == "__main__":
    # Find all the sub directories under split:
    splits_glob = args.db_split_dir + '/' + '*/'
    splits = glob.glob(splits_glob)
    # Paralellise the process:
    pool = multiprocessing.Pool(32)
    output = pool.map(delete_traj, splits)
