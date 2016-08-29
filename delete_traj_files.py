import glob
import multiprocessing
import os


def delete_traj(split):
    prot_glob = split + '/' + '*/'
    prots = glob.glob(prot_glob)
    for prot in prots:
        chain_glob = prot + '/' + '*/'
        chains = glob.glob(chain_glob)
        for chain in chains:
            traj_glob = chain + '/' + 'mutant_traj*'
            trajs = glob.glob(traj_glob)
            for traj in trajs:
                print(traj)
                # os.remove(traj)


split_dir = '/home/projects/cu_10020/data/precalculated_ddg/prot_list.txt/split'
splits_glob = split_dir + '/' + '*/'
splits = glob.glob(splits_glob)

pool = multiprocessing.Pool(28)
output = pool.map(delete_traj, splits)
