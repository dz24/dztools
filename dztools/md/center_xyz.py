from typing import Annotated as And
from typer import Option as Opt
import typer


def center_periodic(
    i: And[str, Opt("-i", help="input")],
    o: And[str, Opt("-o", help="output")],
    c: And[float, Opt("-c", help="The length of a cubic cell")],
    idx: And[int, Opt("-idx", help="The particle idxes at center")
    ],
):
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis import transformations as trans
    box = np.array([c, c, c, 90, 90, 90])
    u = mda.Universe(i, box=box)
    sel = u.atoms[[0, 1]]
    ag = u.select_atoms("all")

    # Save the reference center from the first frame
    u.trajectory[0]
    ref_center = sel.positions.mean(axis=0)
    diff = box[:3]/2 - ref_center
    with mda.Writer(o, ag.n_atoms) as wfile:
        for idx, ts in enumerate(u.trajectory):
            u.atoms.translate(diff)
            ag.wrap(box=box)
            wfile.write(ag.atoms)
            # sel.wrap(box=box)
    exit()

    # Define a transformation that recenters each frame to the reference
    def recenter_to_ref(ts):
        u.atoms.translate(diff)
        sel.wrap(box=box)

        # # trans.wrap(ts)
        # current_center = sel.positions.mean(axis=0)
        # shift = ref_center - current_center
        # ts.positions += shift
        return ts

    # Apply transformation
    u.trajectory.add_transformations(recenter_to_ref)

    # Write out the recentered trajectory
    with mda.Writer(o, u.atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(u.atoms)

    # from dztools.misc.xyz_help import calc_center2
    # atoms = None
    # cnt = 0
    # center = [0, 0, 0]
    # with open(i) as read:
    #     with open(o, "w") as write:
    #         while True:
    #             # read frame header
    #             header = [read.readline(), read.readline()]
    #             if all(i == "" for i in header):
    #                 break
    #             atoms = int(header[0].rstrip())

    #             # read xyz
    #             xyzs = []
    #             for _ in range(atoms):
    #                 line = read.readline().split()
    #                 xyzs.append([line[0]] + [float(i) for i in line[1:]])

    #             # center xyzs
    #             if cnt == 0:
    #                 center = xyzs[idx - 1][1:]
    #             xyzs_c = calc_center2(center, xyzs, c)

    #             # write frame
    #             for line_w in header + xyzs_c:
    #                 write.write(line_w)
    #             cnt += 1
