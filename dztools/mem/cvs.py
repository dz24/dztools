from typing import Annotated

import typer



def mem_helicity(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = None,
    lip: Annotated[str, typer.Option("-lip", help="lipid type")] = "POPC",
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """DSSP"""

    from dztools.misc.mem_help import calc_helicity
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda

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

    if out:
        with open(out, 'w') as write:
            for idx, cv in zip(np.arange(len(hels_avg)), hels_avg):
                write.write(f"{idx}\t{cv*100:.08f}\t{cv:.08f}\n")

    return hels_avg


def mem_com(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="lip type")] = "POPC",
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Calculates COM. Currently for MEL system only."""

    from dztools.misc.mem_help import calc_met_com
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda

    # load gro and xtc into MDA
    u = mda.Universe(top, xtc)
    idxs = list(range(len(u.trajectory)))

    # get coms
    pcoms_z, pcoms_z_avg = calc_met_com(u, lip=lip, num_resi=26)

    # plot
    if plot:
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

    if out:
        with open(out, 'w') as write:
            for idx, cv in zip(np.arange(len(pcoms_z_avg)), pcoms_z_avg):
                write.write(f"{idx}\t{cv:.08f}\n")

    return pcoms_z_avg


def mem_chain(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    hoxy: Annotated[
        str, typer.Option("-hoxy", help="xtc file")
    ] = "O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
    coord_n: Annotated[int, typer.Option("-coord_n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d")] = 1.0,
    coord_r: Annotated[float, typer.Option("-coord_r")] = 8.0,
    coord_z: Annotated[float, typer.Option("-coord_z")] = 0.75,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    from dztools.misc.mem_help import f_axial, f_radial, psi_switch
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    lipid = u.select_atoms(f"resname {lip}")
    epsilons = []
    tpi = 2 * np.pi

    totlen = len(u.trajectory)
    for idx, ts in enumerate(u.trajectory):
        # Frame properties
        z_mem = lipid.atoms.center_of_mass()[-1]
        z_s = z_mem + (np.arange(coord_n) + 1 / 2 - coord_n / 2) * coord_d
        hoxys = u.select_atoms(f"name {hoxy} and prop z < {z_s[-1]+0.5} and prop z > {z_s[0]-0.5}")
        atoms_x = hoxys.atoms.positions[:, 0]
        atoms_y = hoxys.atoms.positions[:, 1]
        box = ts.dimensions

        ws_cyl = [0 for i in range(coord_n)]
        in_axis = [0 for i in range(coord_n)]
        in_radi = [0 for i in range(coord_n)]
        x_sincyl, x_coscyl = 0, 0
        y_sincyl, y_coscyl = 0, 0
        ang_xs, ang_xc = np.sin(tpi * atoms_x / box[0]), np.cos(tpi * atoms_x / box[0])
        ang_ys, ang_yc = np.sin(tpi * atoms_y / box[1]), np.cos(tpi * atoms_y / box[1])

        for s in range(coord_n):
            in_axis[s] = f_axial(hoxys.atoms.positions[:, 2], z_s[s], coord_d)
            f_norm = np.sum(in_axis[s])
            ws_cyl[s] = np.tanh(f_norm)
            if f_norm == 0:
                continue
            x_sincyl += np.sum(in_axis[s] * ang_xs) * ws_cyl[s] / f_norm
            x_coscyl += np.sum(in_axis[s] * ang_xc) * ws_cyl[s] / f_norm
            y_sincyl += np.sum(in_axis[s] * ang_ys) * ws_cyl[s] / f_norm
            y_coscyl += np.sum(in_axis[s] * ang_yc) * ws_cyl[s] / f_norm

        x_sincyl /= np.sum(ws_cyl)
        x_coscyl /= np.sum(ws_cyl)
        y_sincyl /= np.sum(ws_cyl)
        y_coscyl /= np.sum(ws_cyl)
        x_cyl = (np.arctan2(-x_sincyl, -x_coscyl) + np.pi) * box[0] / tpi
        y_cyl = (np.arctan2(-y_sincyl, -y_coscyl) + np.pi) * box[1] / tpi
        print('hunger')
        in_radi = f_radial(hoxys.atoms.positions, x_cyl, y_cyl, coord_r, box)
        epsilon = 0
        nsp = [0 for i in range(coord_n)]
        for s in range(coord_n):
            nsp[s] = np.sum(in_axis[s] * in_radi)
            epsilon += psi_switch(nsp[s], coord_z)
        epsilon /= coord_n
        epsilons.append(epsilon)
        print(idx, epsilon)

    # plot
    if plot:
        plt.plot(np.arange(len(epsilons)), epsilons)
        plt.show()

    if out:
        with open(out, 'w') as write:
            for idx, cv in zip(np.arange(len(epsilons)), epsilons):
                write.write(f"{idx}\t{cv:.08f}\n")

    return np.array(epsilons)


def mem_hel_com(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = False,
    lip: Annotated[str, typer.Option("-lip", help="lip file")] = "POPC",
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Calculate com vs hel"""
    from dztools.misc.mem_help import calc_helicity, calc_met_com
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda

    # load gro and xtc into MDA
    if xtc == "False":
        u = mda.Universe(top)
    else:
        u = mda.Universe(top, xtc)

    print("get helixes")
    hels, hels_avg = calc_helicity(u, num_resi=26)

    # get coms
    print("get coms")
    pcoms_z, pcoms_z_avg = calc_met_com(u, lip=lip, num_resi=26)

    # plot
    print('cheesze', plot)
    if plot:
        # for i, pcom in enumerate(pcoms_z):
        for i, (hel, pcom) in enumerate(zip(hels, pcoms_z)):
            print('peptide', i)
            plt.plot(
                np.array(hel) * 100,
                pcom,
                c=f"C{i+1}",
                label=f"p{i+1}",
                alpha=0.2,
            )
        plt.plot(hels_avg * 100, pcoms_z_avg, color="r", ls="--", lw=2.0)
        plt.xlim([0, 100])
        plt.ylim([0, max(pcoms_z_avg)])
        plt.xlabel("Helicity [%]")
        plt.ylabel("Center Of Mass [z]")
        plt.show()

    if out:
        with open(out, 'w') as write:
            for idx, cv in zip(np.arange(len(hels_avg)), hels_avg):
                write.write(f"{idx}\t{cv:.08f}\n")


def met_rsmd(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")]='',
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Calculate com vs hel"""
    from dztools.misc.mem_help import calc_helicity, calc_met_com
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.analysis.rms import rmsd
    from dztools.misc.mem_help import calc_helicity
    import os
    from MDAnalysis.lib.mdamath import make_whole

    # load gro and xtc into MDA
    if xtc == "False":
        u = mda.Universe(top)
    else:
        u = mda.Universe(top, xtc)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    met = u.select_atoms("protein").select_atoms("backbone")
    # print('wii', len(met.atoms))
    # print('wii', len(met.select_atoms("backbone").atoms))
    # exit()
    make_whole(met)
    base = os.path.dirname(os.path.abspath(__file__))
    cry_u = mda.Universe(base + "/pdbs/2mlt_proa_h.pdb")
    # cry_u.dimensions = u.dimensions
    cry = cry_u.select_atoms("protein").select_atoms("backbone")

    print(met.positions[:, 0])
    print(cry.positions[:, 0])
    mmax = max(met.positions[:, 0])
    cmax = max(cry.positions[:, 0])
    # plt.plot(np.arange(len(met.positions[:, 0])), met.positions[:, 0]/mmax)
    # plt.plot(np.arange(len(met.positions[:, 0])), cry.positions[:, 0]/cmax)
    # plt.show()
    # exit()

    print("get helixes")

    R = mda.analysis.rms.RMSD(u,cry_u,select="protein")
    R.run()
    rmsd = R.results.rmsd.T
    hels, hels_avg = calc_helicity(u, num_resi=26)
    print(len(rmsd[1]), len(set(rmsd[1])), len(hels_avg))
    # plt.plot(np.arange(len(rmsd[1])), rmsd[1])
    # plt.show()
    # exit()
    # plt.plot(rmsd[1], rmsd[2])
    ax1.plot(np.arange(len(rmsd[1])), rmsd[2], color='C0')
    ax2.plot(np.arange(len(rmsd[1])), np.array(hels_avg), color='C1')
    ax1.set_ylabel('rsmd', color='C0')
    ax2.set_ylabel('heli', color='C1')
    # plt.scatter(rmsd[1], rmsd[2])
    # plt.axhline(rsmd0, color='k')
    plt.show()
    print(rmsd[1])
    print(rmsd[2])
    # idx, rsmd = [], []
    #for i, ts in enumerate(u.trajectory):

    print("met", len(met.atoms))
    print("cry", len(cry.atoms))

    if out:
        with open(out, 'w') as write:
            for idx, cv in zip(np.arange(len(hels_avg)), hels_avg):
                write.write(f"{idx}\t{cv:.08f}\n")
