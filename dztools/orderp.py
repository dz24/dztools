import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.lib.distances import calc_angles
import argparse
import numpy as np


def find_neighs_o(c_idx, o_group, xyz, box, rcut):
    dists = distance_array(xyz[c_idx],
                           xyz[o_group],
                           box=box)[0]
    minidxs = np.argpartition(dists, 3)
    neighs = []
    for idx in minidxs[:3]:
        if dists[idx] < rcut:
            neighs.append(o_group[idx])
    return neighs


def count_w0(o_group, h_group, xyz, box):
    dinfo = {i: {'hidxs': [], 'dists': []} for i in o_group}
    for hidx in h_group:
        dists = distance_array(xyz[hidx],
                               xyz[o_group],
                               box=box)[0]

        o_close = o_group[np.argmin(dists)]
        oh_dist = dists[np.argmin(dists)]
        dinfo[o_close]['hidxs'].append(hidx)
        dinfo[o_close]['dists'].append(oh_dist)
    olist = [len(i['hidxs']) for i in dinfo.values()]
    # if 3 in olist:
    #     for key, item in dinfo.items():
    #         print('tiger', key, item)
    dict0 = {}
    for i in range(4):
        if i in olist:
            dict0[i] = sum([j==i for j in olist])

    return dict0


def count_w(arguments):
    parser = argparse.ArgumentParser(
        description="Calculate the theta and phi angle for an \
                .sdf file given a set of indices"
    )

    parser.add_argument(
        "-i", help="Input file"
    )

    parser.add_argument(
        "-c",
        help="The length of a cubic cell",
        type=float,
    )

    args = parser.parse_args(arguments)

    box = [args.c]*3 + [90]*3
    uni = mda.Universe(args.i)
    o_group = []
    h_group = []
    for idx, atom in enumerate(uni.atoms):
        if atom.name == 'O':
            o_group.append(idx)
        if atom.name == 'H':
            h_group.append(idx)
    for frame in uni.trajectory:
        print('frame:', box)
        dic = count_w0(o_group, h_group, frame.positions, box)
        for key, item in dic.items():
            print(key, item)


def co2_op(arguments):
    parser = argparse.ArgumentParser(
        description="Calculate the theta and phi angle for an \
                .sdf file given a set of indices"
    )

    parser.add_argument(
        "-i", help="Input file"
    )

    parser.add_argument(
        "-c",
        help="The length of a cubic cell",
        type=float,
    )

    parser.add_argument(
        "-a",
        help="State A",
        type=float,
    )

    parser.add_argument(
        "-b",
        help="State B",
        type=float,
    )

    args = parser.parse_args(arguments)

    box = [args.c]*3 + [90]*3
    uni = mda.Universe(args.i)
    o_group = []
    h_group = []
    c_idx = 0
    for idx, atom in enumerate(uni.atoms):
        if atom.name == 'O':
            o_group.append(idx)
        if atom.name == 'H':
            h_group.append(idx)
        if atom.name == 'C':
            c_idx = idx
    for idx, frame in enumerate(uni.trajectory):
        pos = frame.positions
        dic = count_w0(o_group, h_group, pos, box)
        # trimer = 1 if 3 in dic else 0
        trimer = 0 if 3 not in dic else dic[3]
        monomers = 0 if 1 not in dic else dic[1]

        o_neighs = find_neighs_o(c_idx,
                                 o_group,
                                 pos,
                                 box,
                                 rcut=1.41)
        angles = []
        if len(o_neighs) == 2:
            angle = calc_angles(pos[o_neighs[0]],
                                pos[c_idx],
                                pos[o_neighs[1]],
                                box=box)
        else:
            angle1 = calc_angles(pos[o_neighs[0]],
                                 pos[c_idx],
                                 pos[o_neighs[1]],
                                 box=box)
            angle2 = calc_angles(pos[o_neighs[0]],
                                 pos[c_idx],
                                 pos[o_neighs[2]],
                                 box=box)
            angle3 = calc_angles(pos[o_neighs[1]],
                                 pos[c_idx],
                                 pos[o_neighs[2]],
                                 box=box)
            angle = min(angle1, angle2, angle3)
        # a and b: pure water, and double OH
        # print(idx, angle, trimer, monomers, len(o_neighs), o_neighs[0], o_neighs[1])

        # idx, min_angle, trimers, monomers, len(o_neighs)
        
        # Check that nothing serious have gone wrong:
        assert len(o_neighs) in (2, 3)
        
        composite_op =  np.pi  - angle + (len(o_neighs) - 2)/4
        fake_op = np.pi - angle + (len(o_neighs) - 2)/4

        # check if stable
        # if trimer == 0:
        #     return [fake_op, composite_op, angle, trimer, monomers, len(o_neighs)]
        if trimer > 0:
            if composite_op < args.a:
                fake_op = args.a+ 0.05
            elif composite_op > args.b:
                fake_op = args.b - 0.05

        print(idx, fake_op, composite_op, angle, trimer, monomers, len(o_neighs))
        # return [fake_op, composite_op, angle, trimer, monomers, len(o_neighs)]
            
            


        # print(idx, angle, trimer, monomers, len(o_neighs))
        # exit('ape')


def calc_dist(arguments):
    parser = argparse.ArgumentParser(
        description="Calculate the theta and phi angle for an \
                .sdf file given a set of indices"
    )

    parser.add_argument(
        "-i", help="Input file"
    )

    parser.add_argument(
        "-c",
        help="The length of a cubic cell",
        type=float,
    )

    parser.add_argument(
        "-idx",
        help="The angle idxes",
        type=int,
        nargs="+",
    )

    args = parser.parse_args(arguments)

    box = [args.c]*3 + [90]*3
    uni = mda.Universe(args.i)
    for idx, frame in enumerate(uni.trajectory):
        xyz = frame.positions
        dists1 = distance_array(xyz[args.idx[0] - 1],
                                xyz[args.idx[1] - 1],
                                box=box)[0][0]
        dists2 = distance_array(xyz[args.idx[0] - 1],
                                xyz[args.idx[2] - 1],
                                box=box)[0][0]
        dists3 = distance_array(xyz[args.idx[0] - 1],
                                xyz[args.idx[3] - 1],
                                box=box)[0][0]
        print(idx, dists1, dists2, dists3)
