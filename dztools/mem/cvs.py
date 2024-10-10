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
    u = mda.Universe(top, xtc)
    idxs = list(range(len(u.trajectory)))

    # get helixes
    hels, hels_avg = calc_helicity(u, num_resi=25)

    # plot
    if plot:
        for hel in hels:
            print(hel)
            plt.plot(idxs, hel, alpha=0.2)
        plt.plot(idxs, hels_avg, color="r")
        plt.show()


def mem_com(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="lip file")] = "POPC",
    plot: Annotated[str, typer.Option("-plot", help="plot")] = True,
):
    """Calculates COM. Currently for MEL system only."""
    import matplotlib.pyplot as plt
    import MDAnalysis as mda

    from dztools.misc.mem_help import calc_met_com

    # load gro and xtc into MDA
    u = mda.Universe(top, xtc)
    idxs = list(range(len(u.trajectory)))

    # get coms
    pcoms_z, pcoms_z_avg, lcoms_z = calc_met_com(u, lip=lip, num_resi=25)

    # plot
    if plot:
        for i, pcom in enumerate(pcoms_z):

            plt.plot(
                idxs,
                abs(pcom - lcoms_z),
                c=f"C{i+1}",
                label=f"p{i+1}",
                alpha=0.2,
            )
        plt.plot(idxs, pcoms_z_avg, color="r", lw=2.0)
        plt.show()


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


def mem_chain(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
    coord_n: Annotated[int, typer.Option("-lip", help="n")] = 26, # ns
    coord_d: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.1*10,
    coord_r: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.9*10,
    coord_z: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.75,
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""
    import matplotlib.pyplot as plt
    import numpy as np

    import MDAnalysis as mda
    from MDAnalysis.analysis import helix_analysis as hel

    from dztools.misc.mem_help import psi_switch, theta, f_axial, f_radial

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)
    # hoxys = u.select_atoms("name OH2 O11 O12 O13 O14")
    # hoxys = u.select_atoms("name OW O11 O12 O13 O14")
    hoxys = u.select_atoms("name OW OA OB OC OD")
    # lipid = u.select_atoms(f"resname POPC")
    lipid = u.select_atoms(f"resname DMPC")

    for idx, ts in enumerate(u.trajectory):
        z_mem = lipid.atoms.center_of_mass()[-1]
        z_s = [z_mem + (s + 1/2 - coord_n/2)*coord_d for s in range(coord_n)]
        box = ts.dimensions

        nsp = [0 for i in range(coord_n)]
        f_norm = [0 for i in range(coord_n)]
        fx_acc_s = [0 for i in range(coord_n)]
        fx_acc_c = [0 for i in range(coord_n)]
        fy_acc_s = [0 for i in range(coord_n)]
        fy_acc_c = [0 for i in range(coord_n)]
        ws_cyl = [0 for i in range(coord_n)]
        theta_ix = [0 for i in range(coord_n)]
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

            fx_acc_s[s] = np.sum(in_axis[s]*np.sin(2*np.pi*hoxys.atoms.positions[:, 0]/box[0]))
            fx_acc_c[s] = np.sum(in_axis[s]*np.cos(2*np.pi*hoxys.atoms.positions[:, 0]/box[0]))
            fy_acc_s[s] = np.sum(in_axis[s]*np.sin(2*np.pi*hoxys.atoms.positions[:, 1]/box[1]))
            fy_acc_c[s] = np.sum(in_axis[s]*np.cos(2*np.pi*hoxys.atoms.positions[:, 1]/box[1]))

            x_sincyl += ws_cyl[s]*fx_acc_s[s]/f_norm[s]
            x_coscyl += ws_cyl[s]*fx_acc_c[s]/f_norm[s]
            y_sincyl += ws_cyl[s]*fy_acc_s[s]/f_norm[s]
            y_coscyl += ws_cyl[s]*fy_acc_c[s]/f_norm[s]
        x_sincyl /= np.sum(ws_cyl)
        x_coscyl /= np.sum(ws_cyl)
        y_sincyl /= np.sum(ws_cyl)
        y_coscyl /= np.sum(ws_cyl)
        x_cyl = (np.arctan2(-x_sincyl,-x_coscyl) + np.pi)*box[0]/(2*np.pi)
        y_cyl = (np.arctan2(-y_sincyl,-y_coscyl) + np.pi)*box[1]/(2*np.pi)

        in_radi = f_radial(hoxys.atoms.positions[:, 0],
                           hoxys.atoms.positions[:, 1],
                           x_cyl, y_cyl, coord_r)

        epsilon = 0
        # print("Slice\tN_p\tpsi(N_p,z)\tComPBCx\tComPBCy")
        for s in range(coord_n):# [:4]:
            nsp[s] = np.sum(in_axis[s]*in_radi)
            epsilon += psi_switch(nsp[s], coord_z)
            # print(f"{s}\t{nsp[s]:.4f}\t{psi_switch(nsp[s], coord_z):.4f}\t\t{x_cyl:0.2f}\t{y_cyl:0.2f}")
        epsilon /= coord_n
        # print('XYZ', f"{x_cyl:0.2f}", f"{y_cyl:0.2f}", z_mem)
        # print("min max cylinder", f"{z_s[0]:0.2f}", f"{z_s[-1]:0.2f}")

        print(idx, epsilon)
