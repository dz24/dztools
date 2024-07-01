import numpy as np
import matplotlib.pyplot as plt
import argparse

import warnings
warnings.filterwarnings('ignore')


import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from dztools.ruru_tools.orderp import calc_orderp
from dztools.ruru_tools.locate import locate

BOX = [12.4138, 12.4138, 12.4138, 90.0, 90.0, 90.0]
AT, O, RU = 193, 64, 2
NAMES = ['Ru']*RU + ['O']*O + ['H']*(AT-O-RU)

def calc_op_ruru(arguments):
    parser = argparse.ArgumentParser(
        description="Calculate ruru op from xyz files."
    )

    parser.add_argument(
        "-i",
        help="xyz and homos",
        nargs="+",
        required=True
    )

    args = parser.parse_args(arguments)
    unis = []
    for inp in args.i:
        unis.append(mda.Universe(inp, box=BOX))

    assert len(unis) == 3

    for i in range(len(unis[0].trajectory)):
        for uni in unis:
            uni.trajectory[i]
        merged = mda.Merge(unis[0].atoms,
                           unis[1].select_atoms('name X'),
                           unis[1].select_atoms('name X'))
        op, _ = calc_orderp(merged.atoms.positions)
        print(f'{i}\t{op:.06f}')


def calc_ruru_dist(arguments):
    parser = argparse.ArgumentParser(
        description="Calculate ruru op from xyz files."
    )

    parser.add_argument(
        "-i",
        help="Input files",
        required=True
    )

    args = parser.parse_args(arguments)
    uni = mda.Universe(args.i, box=BOX)
    for i, xyz in enumerate(uni.trajectory):
        dist = distance_array(xyz[0], xyz[1])[0][0]
        print(f'{i}\t{dist:.06f}')


def calc_ruru_x(arguments):
    parser = argparse.ArgumentParser(
        description="Calculate ruru op from xyz files."
    )

    parser.add_argument(
        "-i",
        help="Input files",
        required=True
    )

    args = parser.parse_args(arguments)
    uni = mda.Universe(args.i, box=BOX)
    x_locs = []
    for i, xyz in enumerate(uni.trajectory):
        # dist = distance_array(xyz[0], xyz[1])[0][0]
        # print(f'{i}\t{dist:.06f}')
        op, x_loc = calc_orderp(xyz)
        x_locs.append(x_loc)
        # op = calc_orderp(merged.atoms.positions)

    x_locs_set = list(set(x_locs))
    print('#\t\t', '\t\t'.join([str(i) for i in x_locs_set]))
    # exit('baka')
    for i, xyz in enumerate(uni.trajectory):
        ops = []
        for xloc in x_locs_set:
            dist1 = distance_array(xyz[0], xyz[xloc], box=BOX)
            dist2 = distance_array(xyz[1], xyz[xloc], box=BOX)
            dist3 = distance_array(xyz[0], xyz[1], box=BOX)
            op = (dist1[0][0]-dist2[0][0])/dist3[0][0]
            ops.append(f"{op: 12.10f}")
        print(f"{i}\t{x_locs[i]}\t", '\t'.join(ops))
        # print(f'{i}\t{op:.06f}')


def ruru_cat(arguments):
    import numpy as np
    import MDAnalysis as mda
    import os
    import argparse

    parser = argparse.ArgumentParser(
            description="Reverse and concatenate .xyz trajectories from an infretis simulati  on.")

    parser.add_argument("-out", help = "the outfile trajectory name")
    parser.add_argument("-traj", help="the traj.txt file. Trajectories should be in same fol  der")
    parser.add_argument("--selection", help="The selection, e.g. 'protein' or 'resname UNL'   (default 'all')",
            default="all")

    args = parser.parse_args(arguments)

    traj_file_arr, index_arr = np.loadtxt(args.traj,
                                          usecols=[1,2],
                                          comments="#",
                                          dtype=str,
                                          unpack=True)
    index_arr = index_arr.astype(int)

    U = {}
    U1 = {}
    U2 = {}
    for traj_file in np.unique(traj_file_arr):
        if not os.path.exists(traj_file):
            raise FileNotFoundError(
                    f"\n No such file {traj_file}, maybe you forgot to 'cd accepted/'?")
        U[traj_file] = mda.Universe(traj_file)
        h1 = traj_file[:-4] + '-HOMO_centers_s1-1.xyz'
        h2 = traj_file[:-4] + '-HOMO_centers_s2-1.xyz'
        U1[h1] = mda.Universe(h1)
        U2[h2] = mda.Universe(h2)
    with open('merged.xyz', 'w') as write:
        for traj_file, index in zip(traj_file_arr,index_arr):
            print('file', traj_file, 'index', index)
            h1 = traj_file[:-4] + '-HOMO_centers_s1-1.xyz'
            h2 = traj_file[:-4] + '-HOMO_centers_s2-1.xyz'
            uni1, uni2, uni3 = U[traj_file], U1[h1], U2[h2]
            uni1.trajectory[index]
            uni2.trajectory[index]
            uni3.trajectory[index]
            # merged = mda.Merge(uni1.atoms,
            #                    uni2.select_atoms('name X'),
            #                    uni3.select_atoms('name X'))
            merged = mda.Merge(uni2.atoms,
                               uni3.select_atoms('name X'))
            names = [i.name for i in merged.atoms]
            poss = merged.atoms.positions
            write.write(f"{len(names)}\n")
            write.write(f"# file: {traj_file}, index: {index}\n")
            for name, pos in zip(names, poss):
                write.write(f"{name}\t{pos[0]: 12.10f}\t{pos[1]: 12.10f}\t{pos[2]: 12.10f}\n")

def find_rel_idxes(arguments):
    parser = argparse.ArgumentParser(
        description="Calculate ruru op from xyz files."
    )

    parser.add_argument(
        "-i",
        help="Input files",
        required=True
    )

    args = parser.parse_args(arguments)
    uni = mda.Universe(args.i, box=BOX)
    for i, xyz in enumerate(uni.trajectory):
        opx_loc = locate(xyz)
        exit('a')
        print(opx_loc)
        # dist = distance_array(xyz[0], xyz[1])[0][0]
        # print(f'{i}\t{dist:.06f}')


def calc_ruru_x(arguments):
    parser = argparse.ArgumentParser(
        description="Calculate ruru op from xyz files."
    )

    parser.add_argument(
        "-i",
        help="Input files",
        required=True
    )

    args = parser.parse_args(arguments)
    uni = mda.Universe(args.i, box=BOX)
    x_locs = []
    for i, xyz in enumerate(uni.trajectory):

        op, x_loc = calc_orderp(xyz)
        x_locs.append(x_loc)
        # op = calc_orderp(merged.atoms.positions)

    x_locs_set = list(set(x_locs))
    print('#\t\t', '\t\t'.join([str(i) for i in x_locs_set]))
    # exit('baka')
    for i, xyz in enumerate(uni.trajectory):
        ops = []
        for xloc in x_locs_set:
            dist1 = distance_array(xyz[0], xyz[xloc], box=BOX)
            dist2 = distance_array(xyz[1], xyz[xloc], box=BOX)
            dist3 = distance_array(xyz[0], xyz[1], box=BOX)
            op = (dist1[0][0]-dist2[0][0])/dist3[0][0]
            ops.append(f"{op: 12.10f}")
        print(f"{i}\t{x_locs[i]}\t", '\t'.join(ops))
