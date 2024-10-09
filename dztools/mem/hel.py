from typing import Annotated

import typer


def helicity(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = None,
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
    plot: Annotated[str, typer.Option("-plot", help="plot")] = False,
):
    """DSSP"""
    import subprocess

    import MDAnalysis as mda
    import numpy as np
    import matplotlib.pyplot as plt

    from dztools.misc.mem_help import calc_helicity

    # load gro and xtc into MDA
    u = mda.Universe(top, xtc)

    # get helixes
    hels, hels_avg = calc_helicity(u, num_resi=25)
    idxs = list(range(len(u.trajectory)))

    # plot
    if plot:
        for hel in hels:
            plt.plot(idxs, hel, alpha=0.2)
        plt.plot(idxs, hels_avg, color='r')
        plt.show()
