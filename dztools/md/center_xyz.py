# from dztools.misc.subparsers import add_io
from dztools.misc.xyz_help import calc_center2
import typer
from typing_extensions import Annotated


def center_periodic(
    i: Annotated[str, typer.Option("-i", help="input")],
    o: Annotated[str, typer.Option("-o", help="output")],
    c: Annotated[float, typer.Option("-c", help="The length of a cubic cell")],
    idx: Annotated[int, typer.Option("-idx", help="The particle idxes at center")],
):
    # parser = argparse.ArgumentParser(
    #     description="Make xyz file periodic for Jmol viewing."
    # )

    # parser.add_argument(
    #     "-c",
    #     help="The length of a cubic cell",
    #     type=float,
    # )

    # parser.add_argument(
    #     "-idx",
    #     help="The particle idxes at center",
    #     type=int,
    # )
    # add_io(parser)

    # args = parser.parse_args(arguments)

    atoms = None
    cnt = 0
    center = [0, 0, 0]
    with open(i, 'r') as read:
        with open(o, 'w') as write:
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
                if cnt == 0:
                    center = xyzs[idx-1][1:]
                xyzs_c = calc_center2(center, xyzs, c)

                # write frame
                for line_w in header + xyzs_c:
                    write.write(line_w)
                cnt += 1
