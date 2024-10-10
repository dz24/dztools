from typing import Annotated

import typer


def mem_helicity(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = None,
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
    plot: Annotated[str, typer.Option("-plot", help="plot")] = False,
):
    """DSSP"""

    import matplotlib.pyplot as plt
    import MDAnalysis as mda

    from dztools.misc.mem_help import calc_helicity

    # load gro and xtc into MDA
    u = mda.Universe(top, xtc, dt=100)
    idxs = list(range(len(u.trajectory)))

    # get helixes
    hels, hels_avg = calc_helicity(u, num_resi=25)

    # plot
    if plot:
        for hel in hels:
            plt.plot(idxs, hel, alpha=0.2)
        plt.plot(idxs, hels_avg, color="r")
        plt.ylabel("Helicity [%]")
        plt.xlabel("Time")
        plt.show()

    return hels_avg


def mem_com(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="lip file")] = "POPC",
    plot: Annotated[str, typer.Option("-plot", help="plot")] = "False",
):
    """Calculates COM. Currently for MEL system only."""
    import matplotlib.pyplot as plt
    import MDAnalysis as mda

    from dztools.misc.mem_help import calc_met_com

    # load gro and xtc into MDA
    u = mda.Universe(top, xtc)
    idxs = list(range(len(u.trajectory)))

    # get coms
    pcoms_z, pcoms_z_avg = calc_met_com(u, lip=lip, num_resi=26)

    # plot
    if "rue" in plot:
        for i, pcom in enumerate(pcoms_z):

            plt.plot(
                idxs,
                pcom,
                c=f"C{i+1}",
                label=f"p{i+1}",
                alpha=0.2,
            )
        plt.plot(idxs, pcoms_z_avg, color="r", lw=2.0)
        plt.ylabel(r"Delta COM_z [$\AA$]")
        plt.xlabel("Time")
        plt.show()

    print("bu", pcoms_z_avg)
    return pcoms_z_avg


def mem_chain(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    hoxy: Annotated[
        str, typer.Option("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
    coord_n: Annotated[int, typer.Option("-coord_n", help="n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d", help="xtc file")] = 1,
    coord_r: Annotated[float, typer.Option("-coord_r", help="xtc file")] = 9,
    coord_z: Annotated[float, typer.Option("-coord_z", help="xtc file")] = 0.75,
    plot: Annotated[int, typer.Option("-plot", help="plot")] = 0,
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""
    import matplotlib.pyplot as plt
    import MDAnalysis as mda
    import numpy as np

    from dztools.misc.mem_help import f_axial, f_radial, psi_switch

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)
    # hoxy "name OH2 O11 O12 O13 O14"
    # hoxy "name OW OA OB OC OD"
    hoxys = u.select_atoms(f"name {hoxy}")
    lipid = u.select_atoms(f"resname {lip}")
    epsilons = []

    for idx, ts in enumerate(u.trajectory):
        z_mem = lipid.atoms.center_of_mass()[-1]
        z_s = [
            z_mem + (s + 1 / 2 - coord_n / 2) * coord_d for s in range(coord_n)
        ]
        box = ts.dimensions

        nsp = [0 for i in range(coord_n)]
        f_norm = [0 for i in range(coord_n)]
        fx_acc_s = [0 for i in range(coord_n)]
        fx_acc_c = [0 for i in range(coord_n)]
        fy_acc_s = [0 for i in range(coord_n)]
        fy_acc_c = [0 for i in range(coord_n)]
        ws_cyl = [0 for i in range(coord_n)]
        in_axis = [0 for i in range(coord_n)]
        in_radi = [0 for i in range(coord_n)]
        x_sincyl, x_coscyl = 0, 0
        y_sincyl, y_coscyl = 0, 0

        for s in range(coord_n):
            in_axis[s] = f_axial(hoxys.atoms.positions[:, 2], z_s[s], coord_d)

            f_norm[s] = np.sum(in_axis[s])
            ws_cyl[s] = np.tanh(f_norm[s])
            if f_norm[s] == 0:
                continue

            fx_acc_s[s] = np.sum( in_axis[s] * np.sin(2 * np.pi * hoxys.atoms.positions[:, 0] / box[0]))
            fx_acc_c[s] = np.sum( in_axis[s] * np.cos(2 * np.pi * hoxys.atoms.positions[:, 0] / box[0]))
            fy_acc_s[s] = np.sum( in_axis[s] * np.sin(2 * np.pi * hoxys.atoms.positions[:, 1] / box[1]))
            fy_acc_c[s] = np.sum( in_axis[s] * np.cos(2 * np.pi * hoxys.atoms.positions[:, 1] / box[1]))

            x_sincyl += ws_cyl[s] * fx_acc_s[s] / f_norm[s]
            x_coscyl += ws_cyl[s] * fx_acc_c[s] / f_norm[s]
            y_sincyl += ws_cyl[s] * fy_acc_s[s] / f_norm[s]
            y_coscyl += ws_cyl[s] * fy_acc_c[s] / f_norm[s]
        x_sincyl /= np.sum(ws_cyl)
        x_coscyl /= np.sum(ws_cyl)
        y_sincyl /= np.sum(ws_cyl)
        y_coscyl /= np.sum(ws_cyl)
        x_cyl = (
            (np.arctan2(-x_sincyl, -x_coscyl) + np.pi) * box[0] / (2 * np.pi)
        )
        y_cyl = (
            (np.arctan2(-y_sincyl, -y_coscyl) + np.pi) * box[1] / (2 * np.pi)
        )

        in_radi = f_radial(
            hoxys.atoms.positions[:, 0],
            hoxys.atoms.positions[:, 1],
            x_cyl,
            y_cyl,
            coord_r,
        )

        epsilon = 0
        for s in range(coord_n):
            nsp[s] = np.sum(in_axis[s] * in_radi)
            epsilon += psi_switch(nsp[s], coord_z)
        epsilon /= coord_n
        epsilons.append(epsilon)

    # plot
    if plot:
        plt.plot(np.arange(len(epsilons)), epsilons)
        plt.show()

    return epsilons


def mem_hel_com(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = False,
    lip: Annotated[str, typer.Option("-lip", help="lip file")] = "POPC",
    plot: Annotated[str, typer.Option("-plot", help="plot")] = True,
):
    """Calculate com vs hel"""
    import matplotlib.pyplot as plt
    import MDAnalysis as mda
    import numpy as np

    from dztools.misc.mem_help import calc_helicity, calc_met_com

    # load gro and xtc into MDA
    if xtc == "False":
        u = mda.Universe(top)
    else:
        u = mda.Universe(top, xtc)

    # get helixes
    hels, hels_avg = calc_helicity(u, num_resi=26)

    # get coms
    pcoms_z, pcoms_z_avg, lcoms_z = calc_met_com(u, lip=lip, num_resi=26)

    # plot
    if plot:
        # for i, pcom in enumerate(pcoms_z):
        for i, (hel, pcom) in enumerate(zip(hels, pcoms_z)):
            plt.plot(
                np.array(hel) * 100,
                pcom - lcoms_z,
                c=f"C{i+1}",
                label=f"p{i+1}",
                alpha=0.2,
            )
        plt.plot(hels_avg * 100, pcoms_z_avg, color="r", ls="--", lw=2.0)
        plt.xlim([0, 100])
        plt.ylim([0, 40])
        plt.xlabel("Helicity [%]")
        plt.ylabel("Center Of Mass [z]")
        plt.show()
