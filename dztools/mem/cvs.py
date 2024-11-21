from typing import Annotated

import typer



def mem_helicity(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = None,
    lip: Annotated[str, typer.Option("-lip", help="lipid type")] = "POPC",
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
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

    return hels_avg


def mem_com(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="lip type")] = "POPC",
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
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

    return pcoms_z_avg


def mem_chain(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    hoxy: Annotated[
        str, typer.Option("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
    coord_n: Annotated[int, typer.Option("-coord_n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d")] = 1.0,
    coord_r: Annotated[float, typer.Option("-coord_r")] = 9.0,
    coord_z: Annotated[float, typer.Option("-coord_z")] = 0.75,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    from dztools.misc.mem_help import f_axial, f_radial, psi_switch
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda

    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    # hoxy examples:
    # "name OH2 O11 O12 O13 O14"
    # "name OW OA OB OC OD"
    lipid = u.select_atoms(f"resname {lip}")
    # lzmax = max(lipid.positions[:, 2])
    # lzmin = min(lipid.positions[:, 2])
    # hoxys = u.select_atoms(f"name {hoxy}")

    # lipid = u.select_atoms(f"r")
    epsilons = []
    tpi = 2 * np.pi

    totlen = len(u.trajectory)
    for idx, ts in enumerate(u.trajectory):
        # Frame properties
        z_mem = lipid.atoms.center_of_mass()[-1]
        z_s = z_mem + (np.arange(coord_n) + 1 / 2 - coord_n / 2) * coord_d
        hoxys = u.select_atoms(f"name {hoxy} and prop z < {z_s[-1]+0.5} and prop z > {z_s[0]-0.5}")
        # hoxys = u.select_atoms(f"resname TIP3 and prop z < {z_s[-1]+0.5} and prop z > {z_s[0]-0.5}")
        atoms_x = hoxys.atoms.positions[:, 0]
        atoms_y = hoxys.atoms.positions[:, 1]
        box = ts.dimensions

        ws_cyl = [0 for i in range(coord_n)]
        in_axis = [0 for i in range(coord_n)]
        in_radi = [0 for i in range(coord_n)]
        x_sincyl, x_coscyl = 0, 0
        y_sincyl, y_coscyl = 0, 0

        plt.axhline(20, ls='--', color='k')
        plt.axhline(box[2]-20, ls='--', color='k')
        for s in range(coord_n):
            in_axis[s] = f_axial(hoxys.atoms.positions[:, 2], z_s[s], coord_d)
            f_norm = np.sum(in_axis[s])
            ws_cyl[s] = np.tanh(f_norm)
            print('ummm', hoxys.atoms.positions[:, 2][0], ws_cyl[s],  f_norm, coord_d)
            if f_norm == 0:
                # print(f"{s:5.0f}\t0\t0\tnan\tnan")
                continue
            # plt.axhline(z_s[s], color='k', alpha=0.5, ls='--')
            # z_in = [i for i, j in zip(hoxys.atoms.positions[:, 2], in_axis[s]) if j > 0]
            # alpha = [j for i, j in zip(hoxys.atoms.positions[:, 2], in_axis[s]) if j > 0]
            # plt.scatter(np.arange(len(z_in)), z_in, alpha=alpha)

            ang_x = tpi * atoms_x / box[0]
            ang_y = tpi * atoms_y / box[1]
            x_sincyl += np.sum(in_axis[s] * np.sin(ang_x)) * ws_cyl[s] / f_norm
            x_coscyl += np.sum(in_axis[s] * np.cos(ang_x)) * ws_cyl[s] / f_norm
            y_sincyl += np.sum(in_axis[s] * np.sin(ang_y)) * ws_cyl[s] / f_norm
            y_coscyl += np.sum(in_axis[s] * np.cos(ang_y)) * ws_cyl[s] / f_norm
            # print(f"{s:5.0f}\t{sum(in_axis[s]):.3f}\t{f_norm:.03f}\t\t{ws_cyl[s]:.3f}")

        x_sincyl /= np.sum(ws_cyl)
        x_coscyl /= np.sum(ws_cyl)
        y_sincyl /= np.sum(ws_cyl)
        y_coscyl /= np.sum(ws_cyl)
        x_cyl = (np.arctan2(-x_sincyl, -x_coscyl) + np.pi) * box[0] / tpi
        y_cyl = (np.arctan2(-y_sincyl, -y_coscyl) + np.pi) * box[1] / tpi
        in_radi = f_radial(atoms_x, atoms_y, x_cyl, y_cyl, coord_r)
        print('bean', atoms_x[0], atoms_y[0], x_cyl, y_cyl, coord_r)

        epsilon = 0
        nsp = [0 for i in range(coord_n)]
        print('')
        print("Slice\tN_p\tpsi(N_p,z)")
        for s in range(coord_n):
            nsp[s] = np.sum(in_axis[s] * in_radi)
            epsilon += psi_switch(nsp[s], coord_z)
            print(f"{s:5.0f}\t{nsp[s]:.5f}\t{psi_switch(nsp[s], coord_z):.5f}")
        print(f"xy: {x_cyl/10:.4f} {y_cyl/10:.4f}, {z_mem/10:.4f}")
        print(f"{z_s[0]/10:.4f}\t{z_s[-1]/10:.4f}")
        epsilon /= coord_n
        epsilons.append(epsilon)
        print(f"frame {idx:5.0f}/{totlen} processeed", epsilon)
        print('flag', box)
        # ax = plt.gca()
        # ax.yaxis.set_inverted(True)

        # ax.set_xlim(0, box[0])
        # ax.set_ylim(0, box[1])
        # ax.set_zlim(0, box[2])
        # Xc,Yc,Zc = data_for_cylinder_along_z(x_cyl,y_cyl,coord_r,z_s[-1]-z_s[0])
        # ax.scatter(hoxys.atoms.positions[:, 0], hoxys.atoms.positions[:, 1], hoxys.atoms.positions[:, 2])
        # ax.plot_surface(Xc, Yc, Zc+z_mem, alpha=0.5)
        # plt.show()
        # exit()
        # exit()

    # plot
    if plot:
        plt.plot(np.arange(len(epsilons)), epsilons)
        plt.show()

    if out:
        with open(out, 'w') as write:
            print('hpew')
            for idx, cv in zip(np.arange(len(epsilons)), epsilons):
                # write.write(f"{idx}\t{2*cv-0.25:.08f}\n")
                write.write(f"{idx}\t{cv:.08f}\n")

    return np.array(epsilons)

def data_for_cylinder_along_z(center_x,center_y,radius,height_z):
    import numpy as np
    z = np.linspace(0, height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid


def mem_hel_com(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = False,
    lip: Annotated[str, typer.Option("-lip", help="lip file")] = "POPC",
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    save: Annotated[bool, typer.Option("-out", help="pdf")] = False,
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
