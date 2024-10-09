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
    idxs = list(range(len(u.trajectory)))

    # get helixes
    hels, hels_avg = calc_helicity(u, num_resi=25)

    # plot
    if plot:
        for hel in hels:
            print(hel)
            plt.plot(idxs, hel, alpha=0.2)
        plt.plot(idxs, hels_avg, color='r')
        plt.show()


def com_met(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="lip file")] = "POPC",
    plot: Annotated[str, typer.Option("-plot", help="plot")] = True,
):
    """Calculates COM. Currently for MEL system only."""
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
        plt.plot(idxs, pcoms_z_avg, color='r', lw=2.0)
        plt.show()

def hel_com(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = False,
    lip: Annotated[str, typer.Option("-lip", help="lip file")] = "POPC",
    plot: Annotated[str, typer.Option("-plot", help="plot")] = True,
):
    """Calculate com vs hel"""
    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    import numpy as np

    from dztools.misc.mem_help import calc_met_com, calc_helicity

    # load gro and xtc into MDA
    if xtc == "False":
        u = mda.Universe(top)
    else:
        u = mda.Universe(top, xtc)
    idxs = list(range(len(u.trajectory)))

    # get helixes
    hels, hels_avg = calc_helicity(u, num_resi=26)

    # get coms
    pcoms_z, pcoms_z_avg, lcoms_z = calc_met_com(u, lip=lip, num_resi=26)

    # plot
    if plot:
        # for i, pcom in enumerate(pcoms_z):
        for i, (hel, pcom) in enumerate(zip(hels, pcoms_z)):
            plt.plot(np.array(hel)*100, pcom - lcoms_z, c=f"C{i+1}",
                     label=f"p{i+1}", alpha=0.2)
        plt.plot(hels_avg*100, pcoms_z_avg, color='r', ls="--", lw=2.0)
        plt.xlim([0, 100])
        plt.ylim([0, 40])
        plt.xlabel("Helicity [%]")
        plt.ylabel("Center Of Mass [z]")
        plt.show()
