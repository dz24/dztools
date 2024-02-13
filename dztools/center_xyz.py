import numpy as np
import argparse


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
    

def center_periodic(arguments):
    parser = argparse.ArgumentParser(
        description="Make xyz file periodic for Jmol viewing."
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
        help="The particle idxes at center",
        type=int,
        # nargs="+",
    )

    parser.add_argument(
        "-o",
        help="Output file",
        nargs='?',
        const='periodic.xyz',
    )


    args = parser.parse_args(arguments)
    print(args)

    atoms = None
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
                xyzs_c = calc_center(args.idx, xyzs, args.c)

                print(header)
                # write frame
                for line_w in header + xyzs_c:
                    write.write(line_w)



            # for idx, line in enumerate(read):
            #     # read frame
            #     rip = line.rstrip().split()
            #     mod = idx%(atoms+2)
            #     if mod in (0, 1):
            #         header.append(line)
            #     else:
            #         xyzs.append(rip)

#                 else:
#                     if mod == args.idx - 2:
#                         center = [float(i) for i in rip[1:]]
#                     xyz_p = [rip[0]]
#                     for dim, center_i in zip([float(i) for i in rip[1:]], center):
#                         dim_c = dim - center_i
#                         if mod == 2:
#                             print(idx, center, dim_c)
#                         if mod == 2:
#                             print(dim, dim_c)
#                         if dim_c < -0.5*boxl:
#                             mult = 1
#                         elif dim_c > 0.5*boxl:
#                             mult = -1
#                         else:
#                             mult = 0
#                         xyz_p.append(f"{dim_c+boxl*mult:10.8f}")

#                    write.write('\t'.join(xyz_p) + '\n')
