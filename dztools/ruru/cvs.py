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
    from dztools.misc.ru_help import get_frame, calc_ruop
    import MDAnalysis as mda
    import glob

    orig = mda.Universe(xyz)
    # list_of_files = [i for i in glob.glob(xyz[:-9] + '*') if "HOMO" in i]
    list_of_files = [i for i in glob.glob(xyz.split("-")[0] + '*') if "HOMO" in i]
    ru_ops = []
    for idx, ts in enumerate(orig.trajectory):
        print("frame", idx, len(orig.trajectory))
        homos = np.concatenate((
            get_frame(orig.atoms.positions, mda.Universe(list_of_files[0])),
            get_frame(orig.atoms.positions, mda.Universe(list_of_files[1])),
        ))

        # calculate ru_op
        ru_op, x_idx, status = calc_ruop(orig, homos)
        ru_ops.append(ru_op)

        # calculate_h_op
        # calculate first solvation_shell
        # calculate h-o-h wholeness, diff between 1 and 2 should b large

    if len(out) > 0:
        np.savetxt(out, np.array([np.arange(len(ru_ops)), ru_ops]).T)
    plt.plot(np.arange(len(ru_ops)), ru_ops)
    plt.scatter(np.arange(len(ru_ops)), ru_ops)
    plt.show()
    return




