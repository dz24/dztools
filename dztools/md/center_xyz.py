import numpy as np
import argparse
from dztools.subparsers import add_io


def calc_center(idx, xyzs, cell):
    xyzs_p = []
    center = xyzs[idx-1][1:]
    for xyz in xyzs:
        xyz_p = [xyz[0]]
        for dim, center_i in zip(xyz[1:], center):
            dim_c = dim - center_i
            if dim_c < -0.5 * cell:
                mult = 1
            elif dim_c > 0.5*cell:
                mult = -1
            else:
                mult = 0
            xyz_p.append(f"{dim_c+cell*mult:10.8f}")
        xyzs_p.append('\t'.join(xyz_p) + '\n')
    return xyzs_p


def calc_center2(center, xyzs, cell):
    xyzs_p = []
    # center = xyzs[idx-1][1:]
    for xyz in xyzs:
        xyz_p = [xyz[0]]
        for dim, center_i in zip(xyz[1:], center):
            dim_c = dim - center_i
            if dim_c < -0.5 * cell:
                mult = 1
            elif dim_c > 0.5*cell:
                mult = -1
            else:
                mult = 0
            xyz_p.append(f"{dim_c+cell*mult:10.8f}")
        xyzs_p.append('\t'.join(xyz_p) + '\n')
    return xyzs_p


def center_periodic(arguments):
    parser = argparse.ArgumentParser(
        description="Make xyz file periodic for Jmol viewing."
    )

    parser.add_argument(
        "-c",
        help="The length of a cubic cell",
        type=float,
    )

    parser.add_argument(
        "-idx",
        help="The particle idxes at center",
        type=int,
    )
    add_io(parser)

    args = parser.parse_args(arguments)

    atoms = None
    cnt = 0
    with open(args.i, 'r') as read:
        with open(args.o, 'w') as write:
            while True:
                # read frame header
                header = [read.readline(), read.readline()]
                if all(i == '' for i in header):
                    break
                atoms = int(header[0].rstrip())

                # read xyz
                xyzs = []
                for _ in range(atoms):
                    line = read.readline().split()
                    xyzs.append([line[0]] + [float(i) for i in line[1:]])

                # center xyzs
                # xyzs_c = calc_center(args.idx, xyzs, args.c)
                if cnt == 0:
                    center = xyzs[args.idx-1][1:]
                xyzs_c = calc_center2(center, xyzs, args.c)

                # write frame
                for line_w in header + xyzs_c:
                    write.write(line_w)
                cnt += 1
