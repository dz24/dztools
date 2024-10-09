from typing import Annotated

import typer


def com_met(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="lip file")] = "POPC",
    plot: Annotated[str, typer.Option("-plot", help="plot")] = True,
):
    """Currently for MEL system only."""
    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    import numpy as np

    from dztools.misc.mem_help import calc_met_com

    # load gro and xtc into MDA
    u = mda.Universe(top, xtc)
    idxs = list(range(len(u.trajectory)))

    # get coms
    pcoms_z, pcoms_z_avg, lcoms_z = calc_met_com(u, lip=lip, num_resi=25)

    # plot
    if plot:
        for i, pcom in enumerate(pcoms_z):

            plt.plot(idxs, abs(pcom - lcoms_z), c=f"C{i+1}",
                     label=f"p{i+1}", alpha=0.2)
        plt.plot(idxs, abs(pcoms_z_avg), color='r', lw=2.0)
        plt.show()
