from typing import Annotated as And

# import typer
from typer import Option as Opt

BOX = [12.4138, 12.4138, 12.4138, 90.0, 90.0, 90.0]

def ru_op(
    xyz: And[str, Opt("-xyz", help="xyz file")] = None,
    plot: And[bool, Opt("-plot", help="plot")] = False,
    out: And[str, Opt("-out", help="string")] = "",
):
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.ru_help import get_frame, calc_ruop, calc_hop
    import MDAnalysis as mda
    import glob

    orig = mda.Universe(xyz)
    list_of_files = [i for i in glob.glob(xyz[:-9] + '*') if "HOMO" in i]
    homo1 = mda.Universe(list_of_files[0])
    homo2 = mda.Universe(list_of_files[1])
    homo11 = homo1.select_atoms("name X")
    homo22 = homo2.select_atoms("name X")
    results = []
    for idx, ts in enumerate(orig.trajectory):
        print("frame", idx, len(orig.trajectory))

        # calculate oh, oho and solvation shell
        hop, h_idx, o_idx, oho, solv1, solv2 = calc_hop(orig)

        homo1.trajectory[idx]
        homo2.trajectory[idx]

        homos = np.concatenate((
            # get_frame(orig.atoms.positions, homo1.atoms.positions),
            # get_frame(orig.atoms.positions, homo2.atoms.positions),
            homo11.atoms.positions,
            homo22.atoms.positions,
        ))

        # calculate ru_op
        ru_op, x_idx, status = calc_ruop(orig, homos)
        results.append([
            idx, ru_op, x_idx, status, hop, h_idx, o_idx, oho, solv1, solv2
        ])

    results = np.array(results)
    if len(out) > 0:
        np.savetxt(
            out,
            results,
            fmt="%s",
            header="idx, ru_op, x_idx, status, hop, h_idx, o_idx, oho, solv1, solv2"
        )

    plt.plot(results[:, 0], results[:, 1])
    plt.scatter(results[:, 0], results[:, 1])
    plt.show()
    return




