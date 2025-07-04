from typing import Annotated

import typer


def calc_plateou(rmax, r_cap, plat=5, power=3):
    if rmax < plat:
        lindec = 1
    elif rmax > r_cap:
        lindec = 0.0
    else:
        lindec = (1/(r_cap)**power)*(r_cap - (rmax-plat))**power
    return lindec



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
    ] = "OH2 O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "resname POPC",
    coord_n: Annotated[int, typer.Option("-coord_n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d")] = 1.0,
    coord_r: Annotated[float, typer.Option("-coord_r")] = 8.0,
    coord_z: Annotated[float, typer.Option("-coord_z")] = 0.75,
    padding: Annotated[float, typer.Option("-padding")] = 0.5,
    coord_h: Annotated[float, typer.Option("-coord_h")] = 0.25,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    from dztools.misc.mem_help import calc_chain
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    lipid = u.select_atoms(f"{lip}")
    epsilons = []
    epsilons_p = []
    tpi = 2 * np.pi

    totlen = len(u.trajectory)
    x, y, box = [], [], []
    for idx, ts in enumerate(u.trajectory):
        # Frame properties
        box.append(ts.dimensions[:3].copy())
        epsilon, epsilon_p, x0, y0, _ = calc_chain(
            u,
            lip=lip,
            hoxy=hoxy,
            coord_n=coord_n,
            coord_d=coord_d,
            coord_z=coord_z,
            coord_h=coord_h)
        epsilons.append(epsilon)
        epsilons_p.append(epsilon_p)
        x.append(x0)
        y.append(y0)

    # plot
    if plot:
        plt.plot(np.arange(len(epsilons)), epsilons)
        plt.show()

    if out:
        with open(out, 'w') as write:
            for idx, cv in zip(np.arange(len(epsilons)), epsilons):
                # write.write(f"{idx}\t{cv:.08f}\n")
                towrite = f"{idx}\t{cv:.08f}\t{x[idx]:.08f}\t{y[idx]:.08f}\t{box[idx][0]:.08f}"
                towrite += f"\t{box[idx][2]:.08f}\t{epsilons_p[idx]:.08f}"
                towrite += "\n"
                write.write(towrite)

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

def mem_gyr(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num_prot: Annotated[int, typer.Option("-num_prot", help="residue number per protein")] = 20,
    radius: Annotated[int, typer.Option("-r", help="radius around center")] = 2,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    from dztools.misc.mem_help import f_axial, f_radial, psi_switch, gyr_com
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda
    # import scienceplots
    # plt.style.use('science')
    # plt.axis('off')
    # plt.box(True)
    plt.tick_params(top='off', bottom='off', left='off', right='off',
                    labelleft='off', labelbottom='off')

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)
    protein = u.select_atoms("protein")
    num_atoms = len(protein.atoms)//num_prot

    # Get the atom indexes of protein starts and ends
    idxes = sorted([num_atoms*i for i in range(num_prot)])
    idxes_e = sorted([num_atoms*(i+1)-1 for i in range(num_prot)])
    # idx_sort = " ".join(sorted([str(i) for i in idxes1 + idxes2]))
    idx_sort = " ".join([str(i) for i in idxes])
    idx_sort_e = " ".join([str(i) for i in idxes_e])
    print("bum ", " ".join([str(i) for i in idxes]))
    print("dum ", " ".join([str(i) for i in idxes_e]))
    print(idx_sort_e)
    atoms = u.select_atoms("index " + idx_sort)
    atoms_e = u.select_atoms("index " + idx_sort_e)

    ax = plt.gca()
    ax.set_frame_on(True)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_aspect('equal', adjustable='box')

    for ts in u.trajectory[-2:]:
    # for ts in u.trajectory:
        box = ts.dimensions
        x_com = gyr_com(atoms, box, 0)
        y_com = gyr_com(atoms, box, 1)

        Rcyc = 50
        radii = np.linspace(0, 2*Rcyc, 100)
        alphaa = []
        for radi in radii:
            pos = np.array([[radi, 0, 0]])
            alphaa.append(f_radial(pos, 0., 0., Rcyc, box, h=.75)[0])

        print("sweat", alphaa)
        for alpha, radi in zip(alphaa, radii):
            circle1 = plt.Circle((x_com, y_com), radi, color='r', fill=False, alpha=alpha)
            print(radi, alpha, box[0])
            ax.add_patch(circle1)

        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for idx in range(len(idxes)):
                    # x0, x1 = atoms.atoms[idx*2].position[0], atoms.atoms[idx*2+1].position[0]
                    # y0, y1 = atoms.atoms[idx*2].position[1], atoms.atoms[idx*2+1].position[1]
                    # print(x0, x1, y0, y1)
                    x0, y0 = atoms.atoms[idx].position[0], atoms.atoms[idx].position[1]
                    x1, y1 = atoms_e.atoms[idx].position[0], atoms_e.atoms[idx].position[1]

                    plt.plot([x0 + box[0]*i, x1 + box[0]*i], [y0 + box[0]*j, y1 + box[0]*j])
                    plt.scatter([x0 + box[0]*i], [y0 + box[0]*j])

        plt.plot([-box[0], box[0]*2], [0]*2, ls='--', color='k', alpha=0.5)
        plt.plot([-box[0], box[0]*2], [box[1]]*2, ls='--', color='k', alpha=0.5)
        plt.plot([0]*2, [-box[0], box[0]*2], ls='--', color='k', alpha=0.5)
        plt.plot([box[0]]*2, [-box[0], box[0]*2], ls='--', color='k', alpha=0.5)

		# calculate circles
        # plt.scatter([x_com], [y_com], color='k', s=5000, facecolors="none")

        plt.scatter([x_com], [y_com], color='k', marker='x')

        plt.plot([x_com-box[0]/2, x_com+box[0]/2], [y_com+box[1]/2]*2, color='r', ls='--', alpha=0.5)
        plt.plot([x_com-box[0]/2, x_com+box[0]/2], [y_com-box[1]/2]*2, color='r', ls='--', alpha=0.5)
        plt.plot([x_com-box[0]/2]*2, [y_com-box[1]/2, y_com+box[1]/2], color='r', ls='--', alpha=0.5)
        plt.plot([x_com+box[0]/2]*2, [y_com-box[1]/2, y_com+box[1]/2], color='r', ls='--', alpha=0.5)

        plt.xlim([-box[0], 2*box[0]])
        plt.ylim([-box[1], 2*box[1]])
        # plt.savefig("gyr.pdf", bbox_inches='tight')
        plt.show()
        exit()


    # for i in range(num_resi):
    #     mel = u.select_atoms(f"resid {1+num_prot*i}:{num_prot*(i+1)}")
    #     atoms.append(mel.atoms[0])
    #     atoms.append(mel.atoms[-1])


        # print(i, len(mel), f"resid {1+num_resi*i}:{num_resi*(i+1)}", num_prot*(i+1)-(1+num_prot*i))
#
#     # lipid = u.select_atoms(f"resname {lip}")
#     epsilons = []
#     tpi = 2 * np.pi
#
#     totlen = len(u.trajectory)
#     for idx, ts in enumerate(u.trajectory):
#         # Frame properties
#         z_mem = lipid.atoms.center_of_mass()[-1]
#         z_s = z_mem + (np.arange(coord_n) + 1 / 2 - coord_n / 2) * coord_d
#         hoxys = u.select_atoms(f"name {hoxy} and prop z < {z_s[-1]+0.5} and prop z > {z_s[0]-0.5}")
#         atoms_x = hoxys.atoms.positions[:, 0]
#         atoms_y = hoxys.atoms.positions[:, 1]
#         box = ts.dimensions
# 
#         ws_cyl = [0 for i in range(coord_n)]
#         in_axis = [0 for i in range(coord_n)]
#         in_radi = [0 for i in range(coord_n)]
#         x_sincyl, x_coscyl = 0, 0
#         y_sincyl, y_coscyl = 0, 0
#         ang_xs, ang_xc = np.sin(tpi * atoms_x / box[0]), np.cos(tpi * atoms_x / box[0])
#         ang_ys, ang_yc = np.sin(tpi * atoms_y / box[1]), np.cos(tpi * atoms_y / box[1])
# 
#         for s in range(coord_n):
#             in_axis[s] = f_axial(hoxys.atoms.positions[:, 2], z_s[s], coord_d)
#             f_norm = np.sum(in_axis[s])
#             ws_cyl[s] = np.tanh(f_norm)
#             if f_norm == 0:
#                 continue
#             x_sincyl += np.sum(in_axis[s] * ang_xs) * ws_cyl[s] / f_norm
#             x_coscyl += np.sum(in_axis[s] * ang_xc) * ws_cyl[s] / f_norm
#             y_sincyl += np.sum(in_axis[s] * ang_ys) * ws_cyl[s] / f_norm
#             y_coscyl += np.sum(in_axis[s] * ang_yc) * ws_cyl[s] / f_norm
# 
#         x_sincyl /= np.sum(ws_cyl)
#         x_coscyl /= np.sum(ws_cyl)
#         y_sincyl /= np.sum(ws_cyl)
#         y_coscyl /= np.sum(ws_cyl)
#         x_cyl = (np.arctan2(-x_sincyl, -x_coscyl) + np.pi) * box[0] / tpi
#         y_cyl = (np.arctan2(-y_sincyl, -y_coscyl) + np.pi) * box[1] / tpi
#         print('hunger')
#         in_radi = f_radial(hoxys.atoms.positions, x_cyl, y_cyl, coord_r, box)
#         epsilon = 0
#         nsp = [0 for i in range(coord_n)]
#         for s in range(coord_n):
#             nsp[s] = np.sum(in_axis[s] * in_radi)
#             epsilon += psi_switch(nsp[s], coord_z)
#         epsilon /= coord_n
#         epsilons.append(epsilon)
#         print(idx, epsilon)
# 
#     # plot
#     if plot:
#         plt.plot(np.arange(len(epsilons)), epsilons)
#         plt.show()
# 
#     if out:
#         with open(out, 'w') as write:
#             for idx, cv in zip(np.arange(len(epsilons)), epsilons):
#                 write.write(f"{idx}\t{cv:.08f}\n")
# 
#     return np.array(epsilons)
# 
# 

def mem_closest(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt

    u = mda.Universe(top, xtc)
    popc = u.select_atoms("name P")
    plen = len(popc)//2

    x = []
    zdeltas = []
    zmems = []
    zavgs = {i: [] for i in range(num)}
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions

        z_mem = popc.atoms.center_of_mass()[2]
        posz = sorted(np.abs(popc.atoms.positions[:, 2] - z_mem))
        maxz = max(posz)

        for i in zavgs.keys():
            zavgs[i].append(np.average(posz[0:i+1])/maxz)

    for i in zavgs.keys():
        plt.plot(x, zavgs[i])
        # zavgs[i].append(np.average(posz[0:i+1])/maxz)

    if out:
        with open(out, 'w') as write:
            # for idx, cv in zip(x, epsilons):
            for idx in x:
                line = f"{idx}"
                for i in zavgs.keys():
                    line += f"\t{zavgs[i][idx]:.4f}"
					# for zavgs[i]
                print(line)
                write.write(line + "\n")
    # zdeltas = np.array(zdeltas)
    # for i in range(len(popc)):
    #     plt.plot(x, zdeltas[:, i])
    # # plt.plot(x, zmems, color="k")
    plt.show()


def mem_perm(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt

    u = mda.Universe(top, xtc)
    popc = u.select_atoms("name P")
    plen = len(popc)//2

    x = []
    y = []
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions

        z_mem = popc.atoms.center_of_mass()[2]
        farthest_lipid_z = sorted(np.abs(popc.atoms.positions[:, 2] - z_mem))[-1]
        zlim_up = z_mem+farthest_lipid_z
        zlim_dw = z_mem-farthest_lipid_z

        hos = u.select_atoms(f"name OH* and prop z < {zlim_up} and prop z > {zlim_dw}")
        posz = sorted(np.abs(hos.atoms.positions[:, 2] - z_mem))

        y.append(-posz[0]/farthest_lipid_z)

    # plt.plot(x, y)

    if out:
        with open(out, 'w') as write:
            # for idx, cv in zip(x, epsilons):
            for idx in x:
                line = f"{idx}\t{y[idx]}"
                print(line)
                write.write(line + "\n")
    plt.show()


def mem_lippair(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.mem_help import gyr_com2

    u = mda.Universe(top, xtc)
    popc = u.select_atoms("name P")
    plen = len(popc)//2

    x  = []
    y1 = []
    y2 = []
    y3 = []
    y4 = []
    y5 = []
    y6 = []
    y7 = []
    y8 = []
    y9 = []
    wavea = np.loadtxt("/home/daniel/Documents/md/ps/chain/128dmpc/gromacs_input/chaingromacs/1unbias/300b/wave.txt")
    wave = np.zeros(len(popc))
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions
        pos = popc.atoms.positions
        y7.append(len(pos)/(box[0]*box[1]))
        z_mem = popc.atoms.center_of_mass()[2]

        z_dist = np.abs(pos[:, 2] - z_mem)
        z_max = max(z_dist)
        z_norm = z_dist/z_max
        y5.append(np.average(sorted(z_norm)[:10]))
        wave += sorted(z_norm)

        z_norm_avg = np.average(z_norm)
        less_than_avg = z_norm[z_norm<z_norm_avg]
        y4.append(np.var(z_norm_avg -less_than_avg))

        closest = np.argmin(z_norm)
        closest_p = distances.distance_array(pos[closest], pos, box=box)[0]
        closest_num = sorted(closest_p)[20]
        closest_numz = z_norm[closest_p < closest_num]
        y6.append(sum(closest_numz))

        z_op = min(z_norm)
        y8.append(sorted(z_norm)[0] + sorted(z_norm)[1])
        summmm = sorted(z_norm)[0] + sorted(z_norm)[1]
        mult = sorted(z_norm)[0]/sorted(z_norm)[1]
        # print(mult)
        # if mult < 1:
        #     print("whada")
        #     exit()
        y9.append(summmm * mult ** 3)
        x_min, y_min = pos[closest][0], pos[closest][1]

        print('bom', sorted(z_norm))
        if sum(z_norm<0.2) > 0:
            y3.append(sum(z_norm<0.2))
            print('cheeze', len(z_norm<0.2))
            clips = pos[z_norm<0.2]
            com_x = gyr_com2(clips, box, 0)
            com_y = gyr_com2(clips, box, 1)
        else:
            y3.append(0)
            com_x = pos[closest][0]
            com_y = pos[closest][1]
        com_d = distances.distance_array(np.array([x_min, y_min, 0]), np.array([com_x, com_y, 0]), box=box)[0][0]
        y1.append(z_op)
        y2.append(com_d)
        print(len(y1), len(y2), len(y3), len(y4))
        if len(set((len(y1), len(y2), len(y3), len(y4)))) > 1:
            print("panda")
            break

    if out:
        with open(out, 'w') as write:
            # for idx, cv in zip(x, epsilons):
            for idx in x:
                line = f"{idx}\t{y1[idx]}\t{y2[idx]}\t{y3[idx]}\t{y4[idx]}\t{y5[idx]}\t{y6[idx]}\t{y7[idx]}\t{y8[idx]}\t{y9[idx]}"
                print(line)
                write.write(line + "\n")
        with open("wave.txt", 'w') as write:
            # for idx, cv in zip(x, epsilons):
            for idx, w in enumerate(wave/len(x)):
                line = f"{idx}\t{w}"
                print(line)
                write.write(line + "\n")
    plt.show()


def mem_rdf(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.mem_help import gyr_com2

    u = mda.Universe(top, xtc)
    popc = u.select_atoms("name P")

    x = []
    # count from 0 - 10
    r_cap = 50
    h_len = 120
    r_lin = np.linspace(0, 50, h_len)
    hist = np.zeros(h_len)
    dens_avg = np.average([0.00041632855753776275, 0.0004151167125516087])
    traj_len2 = len(u.trajectory)
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions
        pos = popc.atoms.positions
        # dens.append(len(popc)/(box[0]*box[1]*box[2]))
        z_mem = popc.atoms.center_of_mass()[2]

        z_dist = np.abs(pos[:, 2] - z_mem)
        z_max = max(z_dist)
        z_norm = z_dist/z_max
        closest = np.argmin(z_norm)
        r_dists = distances.distance_array(popc.atoms.positions[np.argmax(z_norm)], popc.atoms.positions, box=box)[0]
        for r_dist in sorted(r_dists)[1:]:
            # floor
            r_idx = int((r_dist/r_cap)*100)
            if r_idx >= 100:
                print('hehe', r_dist, r_idx, r_lin[-1])
                continue
            hist[r_idx] += 1

    hist_norm = np.zeros(h_len)
    dr = r_lin[1] - r_lin[0]
    traj_len = len(x)
    for idx, hi in enumerate(hist):
        ri0 = r_lin[idx]
        ri1 = ri0 + dr
        div = (4/3)*np.pi*(ri1**3 - ri0**3)*dens_avg
        hist_norm[idx] = hi/(div*traj_len)

    if out:
        with open(out, 'w') as write:
            # for idx, cv in zip(x, epsilons):
            for idx in range(h_len):
                line = f"{r_lin[idx]}\t{hist_norm[idx]}"
                print(line)
                write.write(line + "\n")
    plt.show()

def mem_rdf2(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
    hoxy: Annotated[
        str, typer.Option("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "resname POPC",
    coord_n: Annotated[int, typer.Option("-coord_n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d")] = 1.0,
    coord_r: Annotated[float, typer.Option("-coord_r")] = 8.0,
    coord_z: Annotated[float, typer.Option("-coord_z")] = 0.75,
    padding: Annotated[float, typer.Option("-padding")] = 0.5,
    coord_h: Annotated[float, typer.Option("-coord_h")] = 0.25,
    color: Annotated[str, typer.Option("-color")] = "C0",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.mem_help import f_axial, f_radial, psi_switch

    u = mda.Universe(top, xtc)
    print(u.select_atoms("resid 1"))
    print(u.select_atoms("resid 2"))
    lipid = u.select_atoms("name P")
    tpi = 2 * np.pi

    x = []
    # count from 0 - 10
    r_cap = 50
    h_len = 120
    r_lin = np.linspace(0, r_cap, h_len)
    dr = r_lin[1] - r_lin[0]
    hist = np.zeros(h_len)
    dens_avg = np.average([0.00041632855753776275, 0.0004151167125516087])
    traj_len2 = len(u.trajectory)

    # ax = plt.figure().add_subplot(projection="3d")
    maxh = []
    r_dists_list = []
    blasts = []
    blasts2 = []
    hist_sorted_acc = np.zeros(h_len)

    npsave = []
    minlp = []

    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions
        z_mem = lipid.atoms.center_of_mass()[-1]
        z_s = z_mem + (np.arange(coord_n) + 1 / 2 - coord_n / 2) * coord_d
        hoxys = u.select_atoms(f"name {hoxy} and prop z < {z_s[-1]+padding} and prop z > {z_s[0]-padding}")
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

        print("box", box)
        print('xyc center', x_cyl, y_cyl, z_mem)

        r_dists = distances.distance_array(np.array([x_cyl, y_cyl, z_mem]), lipid.atoms.positions, box=box)[0]
        r_dists_list.append(min(r_dists))
        hist_frame = np.zeros(h_len)
        for r_dist in sorted(r_dists):
            # floor
            r_idx = int((r_dist/r_cap)*100)
            if r_idx >= 100:
                print('hehe', r_dist, r_idx, r_lin[-1])
                continue
            hist_frame[r_idx] += 1

        hist_frame_norm = np.zeros(h_len)
        for idx2, hi in enumerate(hist_frame):
            ri0 = r_lin[idx2]
            ri1 = ri0 + dr
            div = (4/3)*np.pi*(ri1**3 - ri0**3)*dens_avg
            hist_frame_norm[idx2] = hi/div

        # r_dists_list.append(np.average(hist_frame_norm[:10]**2))
        # if r_dists_list[-1] == 0:
        # f_max = hist_frame_norm[hist_frame_norm!=0][0]
        # f_max_idx = np.where(hist_frame_norm==f_max)[0][0]
        # f_max_r = r_lin[f_max_idx]

        npsave.append(hist_frame_norm)
        sorted_hist_frame_norm = sorted(hist_frame_norm)
        hist_sorted_acc += sorted_hist_frame_norm[::-1]

        withinr = 10
        pairs = []
        for r000, h000 in zip(r_lin, hist_frame_norm):
            if r000 > withinr:
                break
            if h000 > 0:
                pairs.append([r000, h000])
        hinsort = sorted([i[1] for i in pairs])
        # blasts.append(0 if not hinsort else max(hinsort)/min(hinsort))
        blasts.append(sum(sorted_hist_frame_norm[::-1][:5]))
        blasts2.append(sum(sorted_hist_frame_norm))
                # if len(pairs) == 2:
                #     twodist = pairs[1][0] - pairs[0][0]
                #     twodistf = 1 if twodist < 2 else -(1/4)*twodist+1.5
        # if len(pairs) < 2:
        #     twodistf = 0.1

        op = 0

        # loner = 0.5 if len(pairs) == 1 else 1

        for pidx, pair in enumerate(pairs):

            lindec = calc_plateou(pair[0], r_cap=10, plat=2.5, power=1)
            if pair[0] != pairs[-1][0]:
                twodist = pairs[pidx+1][0] - pair[0]
                if pidx == 0:
                    twodist = min([twodist, pair[0]])
                    print("sel", [twodist, pair[0]], twodist)
                lindec2 = calc_plateou(twodist, r_cap=5, plat=2.0, power=3)
                # print("blast", max(hinsort)/min(hinsort))
            else:
                twodist = None
                lindec2 = 1

            op1 = pair[1]*lindec*lindec2
            op += op1
            print('hhh', pair[0], pair[1], lindec, lindec2, twodist, op1, op)
            # print('aya', twodist, lindec2)

        # f_max1 = sorted_hist_frame_norm[-1]
        # f_max_idx1 = np.where(hist_frame_norm==f_max1)[0][0]
        # f_max_r1 = r_lin[f_max_idx1]

        # f_max2 = sorted_hist_frame_norm[-2]
        # f_max_idx2 = np.where(hist_frame_norm==f_max2)[0][0]
        # f_max_r2 = r_lin[f_max_idx2]

        # print('bea', f_max1, f_max_r1, f_max2, f_max_r2)
        # exit()

        # argmax = np.argmax(hist_frame_norm)
        # maxdiff = (f_max/max(hist_frame_norm))**2
        # plateou = 5
        # if f_max_r < plateou:
        #     lindec = 1
        # else:
        #     lindec = (1/(r_cap)**3)*(r_cap - (f_max_r-plateou))**3
        # maxh.append(lindec * maxdiff * f_max)

        # lindec_1 = calc_plateou(f_max_r1, 15)
        # lindec_2 = calc_plateou(f_max_r2, 15)
        # op = lindec_1 * f_max1 + lindec_2 * f_max2
        # print('tiger', lindec_1, f_max1, lindec_2, f_max2)
        # exit()
        maxh.append(200 if op > 200 else op)
        # if op > 200:
            # maxh.append(200)

        # if idx%100 == 0:
        # if 10 < maxh[-1] < 15:
        # if  35 > maxh[-1] > 20:
        # if  30 > maxh[-1] > 22:
        # if  125 > blasts[-1] > 115:
        # if  180 > blasts2[-1] > 160:
        # if True:
        if False:
            print("tf")
            # print(hist_frame_norm)
            # print(idx, f_max1, lindec, f_max, maxh[-1])
            print("ope", maxh[-1])
            # print("twodist", twodist, twodistf)
            # plt.plot(r_lin, np.log(hist_frame_norm+0.5))

            # plt.plot(r_lin, sorted(hist_frame_norm)[::-1], color=color)
            # plt.scatter(r_lin, sorted(hist_frame_norm)[::-1], color=color)

            # temp00 = np.array(sorted(hist_frame_norm)[::-1])[:20]
            # plt.plot(r_lin[:19], temp00[:-1]/temp00[1:], color=color, ls='--')
            # plt.scatter(r_lin[:19], temp00[:-1]/temp00[1:], color=color, marker='x')

            # lip_com calc
            p


            plt.show()
            # exit()
        # if idx%100 == 0:
        # plt.plot(r_lin, hist_frame_norm, idx, alpha=0.5, color='k')
        # for x00, y00 in zip(r_lin, hist_frame_norm):
        #     if y00 == 0:
        #         continue
        #     write.write(f"Ar\t{x00}\t{y00}\t{idx}\n")

        hist += hist_frame

    # np.save(out[:4], npsave)
    # exit()
    hist_norm = np.zeros(h_len)
    # dens_avg = np.average(dens)
    # print("dens_avg", dens_avg)
    traj_len = len(x)
    for idx, hi in enumerate(hist):
        ri0 = r_lin[idx]
        ri1 = ri0 + dr
        div = (4/3)*np.pi*(ri1**3 - ri0**3)*dens_avg
        hist_norm[idx] = hi/(div*traj_len)

    plt.plot(r_lin, hist_norm, color='k')
    if out:
        with open(out, 'w') as write:
            for idx, cv in zip(x, maxh):
            # for x00, x11 in zip(r_lin, hist_norm):
            # for idx in range(h_len):
                line = f"{idx}\t{cv}\t{r_dists_list[idx]}\t{blasts[idx]}\t{blasts2[idx]}"
                # line = f"{x00}\t{x11}"
                print(line)
                write.write(line + "\n")
    plt.show()
    print("gorilla")
    for x, y in zip(r_lin, hist_sorted_acc/traj_len):
        print(x, y)


def mem_rdf3(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
    hoxy: Annotated[
        str, typer.Option("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "resname POPC",
    coord_n: Annotated[int, typer.Option("-coord_n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d")] = 1.0,
    coord_r: Annotated[float, typer.Option("-coord_r")] = 8.0,
    coord_z: Annotated[float, typer.Option("-coord_z")] = 0.75,
    padding: Annotated[float, typer.Option("-padding")] = 0.5,
    coord_h: Annotated[float, typer.Option("-coord_h")] = 0.25,
    color: Annotated[str, typer.Option("-color")] = "C0",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.mem_help import f_axial, f_radial, psi_switch

    u = mda.Universe(top, xtc)
    lipid = u.select_atoms("name P")
    tpi = 2 * np.pi

    x = []
    # count from 0 - 10
    r_cap = 30
    h_len = 120
    r_lin = np.linspace(0, r_cap, h_len)
    dr = r_lin[1] - r_lin[0]
    hist = np.zeros(h_len)
    dens_avg = np.average([0.00041632855753776275, 0.0004151167125516087])
    traj_len2 = len(u.trajectory)

    maxh = []
    r_dists_list = []
    blasts = []
    blasts2 = []
    hist_sorted_acc = np.zeros(h_len)
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions
        z_mem = lipid.atoms.center_of_mass()[-1]
        z_s = z_mem + (np.arange(coord_n) + 1 / 2 - coord_n / 2) * coord_d
        hoxys = u.select_atoms(f"name {hoxy} and prop z < {z_s[-1]+padding} and prop z > {z_s[0]-padding}")
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


        print("box", box)
        print('xyc center', x_cyl, y_cyl, z_mem)

        pcoms = []
        for i in range(len(lipid)):
            resid = u.select_atoms(f"resid {i+1}")
            pcoms.append(resid.center_of_mass(unwrap=True))
        r_dists = distances.distance_array(np.array([x_cyl, y_cyl, z_mem]),
                                           np.array(pcoms), box=box)[0]
        # plt.plot(np.arange(len(r_dists)), r_dists)
        # plt.plot(np.arange(len(r_dists)), sorted(r_dists))
        # plt.show()
        # print("sorted", sorted(r_dists))
        # exit()

        hist_frame = np.zeros(h_len)
        for r_dist in sorted(r_dists):
            r_idx = int((r_dist/r_cap)*100)
            if r_idx >= 100:
                print('hehe', r_dist, r_idx, r_lin[-1])
                continue
            hist_frame[r_idx] += 1

        hist_frame_norm = np.zeros(h_len)
        for idx2, hi in enumerate(hist_frame):
            ri0 = r_lin[idx2]
            ri1 = ri0 + dr
            div = (4/3)*np.pi*(ri1**3 - ri0**3)*dens_avg
            hist_frame_norm[idx2] = hi/div

        r_dists_list.append(np.average(hist_frame_norm[:10]**2))

        # if r_dists_list[-1] == 0:
        # f_max = hist_frame_norm[hist_frame_norm!=0][0]
        # f_max_idx = np.where(hist_frame_norm==f_max)[0][0]
        # f_max_r = r_lin[f_max_idx]

        sorted_hist_frame_norm = sorted(hist_frame_norm)
        hist_sorted_acc += sorted_hist_frame_norm[::-1]

        withinr = 10
        pairs = []
        for r000, h000 in zip(r_lin, hist_frame_norm):
            if r000 > withinr:
                break
            if h000 > 0:
                pairs.append([r000, h000])
        hinsort = sorted([i[1] for i in pairs])
        # blasts.append(0 if not hinsort else max(hinsort)/min(hinsort))
        blasts.append(sum(sorted_hist_frame_norm[::-1][:5]))
        blasts2.append(sum(sorted_hist_frame_norm))
                # if len(pairs) == 2:
                #     twodist = pairs[1][0] - pairs[0][0]
                #     twodistf = 1 if twodist < 2 else -(1/4)*twodist+1.5
        # if len(pairs) < 2:
        #     twodistf = 0.1

        op = 0

        # loner = 0.5 if len(pairs) == 1 else 1

        for pidx, pair in enumerate(pairs):

            lindec = calc_plateou(pair[0], r_cap=10, plat=2.5, power=1)
            if pair[0] != pairs[-1][0]:
                twodist = pairs[pidx+1][0] - pair[0]
                if pidx == 0:
                    twodist = min([twodist, pair[0]])
                    print("sel", [twodist, pair[0]], twodist)
                lindec2 = calc_plateou(twodist, r_cap=5, plat=2.0, power=3)
                # print("blast", max(hinsort)/min(hinsort))
            else:
                twodist = None
                lindec2 = 1

            op1 = pair[1]*lindec*lindec2
            op += op1
            print('hhh', pair[0], pair[1], lindec, lindec2, twodist, op1, op)
            # print('aya', twodist, lindec2)

        # f_max1 = sorted_hist_frame_norm[-1]
        # f_max_idx1 = np.where(hist_frame_norm==f_max1)[0][0]
        # f_max_r1 = r_lin[f_max_idx1]

        # f_max2 = sorted_hist_frame_norm[-2]
        # f_max_idx2 = np.where(hist_frame_norm==f_max2)[0][0]
        # f_max_r2 = r_lin[f_max_idx2]

        # print('bea', f_max1, f_max_r1, f_max2, f_max_r2)
        # exit()

        # argmax = np.argmax(hist_frame_norm)
        # maxdiff = (f_max/max(hist_frame_norm))**2
        # plateou = 5
        # if f_max_r < plateou:
        #     lindec = 1
        # else:
        #     lindec = (1/(r_cap)**3)*(r_cap - (f_max_r-plateou))**3
        # maxh.append(lindec * maxdiff * f_max)

        # lindec_1 = calc_plateou(f_max_r1, 15)
        # lindec_2 = calc_plateou(f_max_r2, 15)
        # op = lindec_1 * f_max1 + lindec_2 * f_max2
        # print('tiger', lindec_1, f_max1, lindec_2, f_max2)
        # exit()
        maxh.append(200 if op > 200 else op)
        # if op > 200:
            # maxh.append(200)

        # if idx%100 == 0:
        # if 10 < maxh[-1] < 15:
        # if  35 > maxh[-1] > 20:
        # if  30 > maxh[-1] > 22:
        # if  125 > blasts[-1] > 115:
        # if  180 > blasts2[-1] > 160:
        # if True:
        if False:
            print("tf")
            # print(hist_frame_norm)
            # print(idx, f_max1, lindec, f_max, maxh[-1])
            print("ope", maxh[-1])
            # print("twodist", twodist, twodistf)
            # plt.plot(r_lin, np.log(hist_frame_norm+0.5))

            # plt.plot(r_lin, sorted(hist_frame_norm)[::-1], color=color)
            # plt.scatter(r_lin, sorted(hist_frame_norm)[::-1], color=color)

            # temp00 = np.array(sorted(hist_frame_norm)[::-1])[:20]
            # plt.plot(r_lin[:19], temp00[:-1]/temp00[1:], color=color, ls='--')
            # plt.scatter(r_lin[:19], temp00[:-1]/temp00[1:], color=color, marker='x')

            # lip_com calc
            p


            plt.show()
            # exit()
        # if idx%100 == 0:
        # plt.plot(r_lin, hist_frame_norm, idx, alpha=0.5, color='k')
        # for x00, y00 in zip(r_lin, hist_frame_norm):
        #     if y00 == 0:
        #         continue
        #     write.write(f"Ar\t{x00}\t{y00}\t{idx}\n")

        hist += hist_frame

    # exit()
    hist_norm = np.zeros(h_len)
    # dens_avg = np.average(dens)
    # print("dens_avg", dens_avg)
    traj_len = len(x)
    for idx, hi in enumerate(hist):
        ri0 = r_lin[idx]
        ri1 = ri0 + dr
        div = (4/3)*np.pi*(ri1**3 - ri0**3)*dens_avg
        hist_norm[idx] = hi/(div*traj_len)

    plt.plot(r_lin, hist_norm, color='k')
    if out:
        with open(out, 'w') as write:
            for idx, cv in zip(x, maxh):
            # for x00, x11 in zip(r_lin, hist_norm):
            # for idx in range(h_len):
                line = f"{idx}\t{cv}\t{r_dists_list[idx]}\t{blasts[idx]}\t{blasts2[idx]}"
                # line = f"{x00}\t{x11}"
                print(line)
                write.write(line + "\n")
    plt.show()
    print("gorilla")
    for x, y in zip(r_lin, hist_sorted_acc/traj_len):
        print(x, y)


def mem_lin(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
    color: Annotated[str, typer.Option("-color")] = "C0",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.mem_help import f_axial, f_radial, psi_switch

    u = mda.Universe(top, xtc)
    lipid = u.select_atoms("name P")

    x = []
    # up_max, up_min, dw_max, dw_min = [], [], [], []
    zmax, zmin = [], []
    # smal1, smal2 = [], []
    z_mems, z_0, z_1 = [], [], []
    z_2, z_3 = [], []
    z_max, z_min = [], []
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions
        z_mem = lipid.atoms.center_of_mass()[-1]
        z_mems.append(z_mem)

        # l_up = lipid.select_atoms(f"prop z > {z_mem}")
        # l_dw = lipid.select_atoms(f"prop z < {z_mem}")

        pos = lipid.atoms.positions[:, 2]
        z_max.append(max(pos))
        z_min.append(min(pos))
        z_dist = np.abs(pos - z_mem)
        z_idxes = np.argsort(z_dist)
        z_idxes2 = np.argsort(pos)

        z_0.append(pos[z_idxes[0]])
        z_1.append(pos[z_idxes[1]])

        z_hmax = z_idxes2[len(z_idxes2)//2]
        z_hmin = z_idxes2[len(z_idxes2)//2-1]
        z_2.append(pos[z_hmax])
        z_3.append(pos[z_hmin])

        # dist = distances.distance_array(pos[z_hmin], pos[z_hmax], box=box)[0][0]
        # dz_list.append(np.abs(pos[z_hmin][2] - pos[z_hmax][2]))

    if out:
        with open(out, 'w') as write:
            for idx in x:
                # line = f"{idx}\t{up_max[idx]}\t{up_min[idx]}\t{dw_max[idx]}\t{dw_min[idx]}\t{smal1[idx]}\t{smal2[idx]}"
                line = f"{idx}\t{z_mems[idx]}\t{z_0[idx]}\t{z_1[idx]}"
                line += f"\t{z_2[idx]}\t{z_3[idx]}\t{z_max[idx]}\t{z_min[idx]}"
                # line = f"{x00}\t{x11}"
                print(line)
                write.write(line + "\n")

def mem_sphere(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
    color: Annotated[str, typer.Option("-color")] = "C0",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.mem_help import f_axial, f_radial, psi_switch

    u = mda.Universe(top, xtc)
    lipid = u.select_atoms("name PH")

    x = []
    dist = []
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions
        z_mem = lipid.atoms.center_of_mass()
        up_dists = distances.distance_array(z_mem, lipid.atoms.positions, box=box)[0]
        dist.append(max(up_dists))
    # plt.plot(x, dist)
    # plt.show()

    if out:
        with open(out, 'w') as write:
            for idx in x:
                line = f"{idx}\t{dist[idx]}"
                # line = f"{x00}\t{x11}"
                print(line)
                write.write(line + "\n")


def mem_flat(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
    color: Annotated[str, typer.Option("-color")] = "C0",
):
    """Plot how flat lipids are"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.mem_help import f_axial, f_radial, psi_switch

    u = mda.Universe(top, xtc)
    lipid = u.select_atoms("name P")
    arange = np.arange(1,1+len(lipid))
    lipids = [u.select_atoms(f"resid {i}") for i in arange]
    allpid = u.select_atoms(f"resid {arange[0]}-{arange[-1]}")

    x = []
    dist1 = []
    dist2 = []
    dist3 = []
    dist4 = []
    dist5 = []
    dist6 = []
    zcom1 = []
    zcom2 = []
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions
        zlist = []
        zcom1.append(lipid.atoms.center_of_mass()[-1])
        z_mem = lipid.atoms.center_of_mass()[-1]

        pos = lipid.atoms.positions
        z_idxes = np.argsort(pos[:, 2])

        dz_list = []
        for ii in range(5):
            z_hmax = z_idxes[len(z_idxes)//2 +ii]
            z_hmin = z_idxes[len(z_idxes)//2-(1+ii)]
            dist = distances.distance_array(pos[z_hmin], pos[z_hmax], box=box)[0][0]
            dz_list.append(np.abs(pos[z_hmin][2] - pos[z_hmax][2]))

        if False:
            ax = plt.figure().add_subplot(projection='3d')
            lower = z_idxes[:len(z_idxes)//2]
            highe = z_idxes[len(z_idxes)//2:]

            ax.scatter(pos[lower, 0], pos[lower, 1], pos[lower, 2], color='C0')
            ax.scatter(pos[highe, 0], pos[highe, 1], pos[highe, 2], color='C1')

            print("no cap", len(pos[len(z_idxes)//2:, 2]), len(pos[:len(z_idxes)//2, 0]))
            plt.show()
        dist1.append(dist)
        dist2.append(dz_list[0])
        dist3.append(dz_list[1])
        dist4.append(dz_list[2])
        dist5.append(dz_list[3])
        dist6.append(dz_list[4])

    if out:
        with open(out, 'w') as write:
            for idx in x:
                line = f"{idx}\t{dist1[idx]}\t{dist2[idx]}\t{dist3[idx]}\t{dist4[idx]}\t{dist5[idx]}\t{dist6[idx]}\t{zcom1[idx]}"

                write.write(line + "\n")


def mem_ff3(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
    color: Annotated[str, typer.Option("-color")] = "C0",
):
    """Plot how flat lipids are"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.mem_help import f_axial, f_radial, psi_switch

    u = mda.Universe(top, xtc)
    lipid = u.select_atoms("name P")
    arange = np.arange(1,1+len(lipid))
    lipids = [u.select_atoms(f"resid {i}") for i in arange]
    allpid = u.select_atoms(f"resid {arange[0]}-{arange[-1]}")

    x = []
    ops1 = []
    ops2 = []
    ops3 = []
    ops4 = []
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions
        # zlist = []
        # zcom1.append(lipid.atoms.center_of_mass()[-1])
        # z_mem = lipid.atoms.center_of_mass()[-1]

        pos = lipid.atoms.positions
        z_idxes = np.argsort(pos[:, 2])
        lower = pos[:, 2][z_idxes[:len(z_idxes)//2]]
        upper = pos[:, 2][z_idxes[len(z_idxes)//2:]]
        up_avg = np.average(upper)
        dw_avg = np.average(lower)
        mlen = up_avg - dw_avg
        dup_avg = abs(min(upper)-up_avg)/mlen
        ddw_avg = abs(max(lower)-dw_avg)/mlen
        op1 = max([dup_avg, ddw_avg])
        op2 = sum([dup_avg, ddw_avg])
        op3 = abs(dup_avg - ddw_avg)

        ops1.append(op1)
        ops2.append(op2)
        ops3.append(op3 + op1)
        ops4.append(op3 + op2)

    if out:
        with open(out, 'w') as write:
            for idx in x:
                line = f"{idx}\t{ops1[idx]}\t{ops2[idx]}\t{ops3[idx]}\t{ops4[idx]}"

                write.write(line + "\n")

def mem_leafneig(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt

    u = mda.Universe(top, xtc)
    popc = u.select_atoms("name P")
    plen = len(popc)//2

    x = []
    y = []
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions

        pos = lipid.atoms.positions
        z_idxes = np.argsort(pos[:, 2])
        lower_idx = z_idxes[:len(z_idxes)//2]
        upper_idx = z_idxes[len(z_idxes)//2:]
        lower = pos[:, 2][lower_idx]
        upper = pos[:, 2][upper_idx]

        dists = distances.distance_array(lower, upper, box=box)[0]

        z_mem = popc.atoms.center_of_mass()[2]
        d_pot = sorted(np.abs(pot.atoms.positions[:, 2] - z_mem))[0]
        y.append(-d_pot)

    if out:
        with open(out, 'w') as write:
            # for idx, cv in zip(x, epsilons):
            for idx in x:
                line = f"{idx}\t{y[idx]}"
                print(line)
                write.write(line + "\n")
    plt.show()


def mem_ff4(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Plot how flat lipids are"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    from dztools.misc.mem_help import f_axial, f_radial, psi_switch

    print(top, xtc)
    u = mda.Universe(top, xtc)
    lipid = u.select_atoms("name P")
    p_upper = u.select_atoms("name P and resid 1 to 64")
    p_lower = u.select_atoms("name P and resid 65 to 128")

    x = []
    ops1 = []
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions

        com_u = p_upper.atoms.center_of_mass()[2]
        com_l = p_lower.atoms.center_of_mass()[2]

        op_1 = np.abs(min(p_upper.atoms.positions[:, 2]) - com_u)
        op_2 = np.abs(max(p_lower.atoms.positions[:, 2]) - com_l)

        op1 = op_1 + op_2
        ops1.append(op1)

    if out:
        with open(out, 'w') as write:
            for idx in x:
                line = f"{idx}\t{ops1[idx]}"
                write.write(line + "\n")

def mem_cff(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    hoxy: Annotated[
        str, typer.Option("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "resname DMPC",
    coord_n: Annotated[int, typer.Option("-coord_n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d")] = 1.0,
    coord_r: Annotated[float, typer.Option("-coord_r")] = 8.0,
    coord_z: Annotated[float, typer.Option("-coord_z")] = 0.75,
    padding: Annotated[float, typer.Option("-padding")] = 0.5,
    coord_h: Annotated[float, typer.Option("-coord_h")] = 0.25,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    from dztools.misc.mem_help import calc_chain
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    lipid = u.select_atoms(f"{lip}")
    lipid_p = u.select_atoms(f"name P")
    epsilons = []
    tpi = 2 * np.pi

    totlen = len(u.trajectory)
    ff = []
    ff2 = []
    cff = []
    cffe = []
    z_mins = []
    for idx, ts in enumerate(u.trajectory):
        # Frame properties
        epsilon, epsilon_e, _, _, _ = calc_chain(u, lip=lip, coord_r = coord_r, coord_n=coord_n)
        epsilons.append(epsilon)
        cffe.append(epsilon_e)

        # flipflop 2
        lip_z = lipid_p.atoms.positions[:, 2]
        z_mem = lipid.atoms.center_of_mass()[-1]
        z_mins.append(np.min(np.abs(lip_z - z_mem)))
        z_idxes = np.argsort(lip_z)
        dz_list = []
        for ii in range(2):
            z_hmax = z_idxes[len(z_idxes)//2 +ii]
            z_hmin = z_idxes[len(z_idxes)//2-(1+ii)]
            dz_list.append(np.abs(lip_z[z_hmin] - lip_z[z_hmax]))

        # OP:
        # lambda = epsilon - b*ff
        b = 0.01
        op = epsilon_e - b*(dz_list[0] + dz_list[1])
        ff.append(-b*dz_list[0])
        ff2.append(-b*dz_list[1])

        # IF IN UNSTABLE FF STATE, GO TO B...
        if epsilon < 0.6 and -dz_list[0] > -5:
            op1 = 0.99
        else:
            op1 = op

        cff.append(op1)

    # plot
    if plot:
        plt.plot(np.arange(len(epsilons)), epsilons)
        plt.show()

    if out:
        with open(out, 'w') as write:
            for idx, cv in zip(np.arange(len(epsilons)), epsilons):
                write.write(f"{idx}\t{cff[idx]:.08f}\t{ff[idx]:.08f}\t{ff2[idx]:.08f}\t{cv:.08f}\t{cffe[idx]:.08f}\t{z_mins[idx]:.08f}\n")

    return np.array(epsilons)



def mem_rdf4(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
    hoxy: Annotated[
        str, typer.Option("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "resname DMPC",
    coord_n: Annotated[int, typer.Option("-coord_n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d")] = 1.0,
    coord_r: Annotated[float, typer.Option("-coord_r")] = 8.0,
    coord_z: Annotated[float, typer.Option("-coord_z")] = 0.75,
    padding: Annotated[float, typer.Option("-padding")] = 0.5,
    coord_h: Annotated[float, typer.Option("-coord_h")] = 0.25,
    color: Annotated[str, typer.Option("-color")] = "C0",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.mem_help import calc_chain

    u = mda.Universe(top, xtc)
    lipid = u.select_atoms("name P")
    tpi = 2 * np.pi

    x = []
    # count from 0 - 10
    r_cap = 30
    h_len = 120
    r_lin = np.linspace(0, r_cap, h_len)
    dr = r_lin[1] - r_lin[0]
    hist = np.zeros(h_len)
    dens_avg = np.average([0.00041632855753776275, 0.0004151167125516087])
    traj_len2 = len(u.trajectory)

    maxh = []
    r_dists_list = []
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions
        z_mem = lipid.atoms.center_of_mass()[-1]

        epsilon, _, x_cyl, y_cyl = calc_chain(u, lip=lip)

        # flipflop 2
        lip_z = lipid.atoms.positions[:, 2]
        z_idxes = np.argsort(lip_z)
        dz_list = []
        for ii in range(1):
            z_hmax = z_idxes[len(z_idxes)//2 +ii]
            z_hmin = z_idxes[len(z_idxes)//2-(1+ii)]
            dz_list.append(np.abs(lip_z[z_hmin] - lip_z[z_hmax]))

        b = 0.01
        op = epsilon - b*dz_list[0]
        if op > 0.3:
            print(epsilon, b*dz_list[0])
            continue

        z_s = z_mem + (np.arange(coord_n) + 1 / 2 - coord_n / 2) * coord_d
        r_dists = distances.distance_array(np.array([x_cyl, y_cyl, z_mem]),
                                           lipid.atoms.positions, box=box)[0]
        for r_dist in sorted(r_dists):
            r_idx = int((r_dist/r_cap)*100)
            if r_idx >= 100:
                continue
            hist[r_idx] += 1

    hist_norm = np.zeros(h_len)
    dr = r_lin[1] - r_lin[0]
    traj_len = len(x)
    for idx, hi in enumerate(hist):
        ri0 = r_lin[idx]
        ri1 = ri0 + dr
        div = (4/3)*np.pi*(ri1**3 - ri0**3)
        hist_norm[idx] = hi/(div*traj_len)
    print("cheeze")

    plt.plot(r_lin, hist_norm, color='k')
    if out:
        with open(out, 'w') as write:
            for idx in range(h_len):
                line = f"{r_lin[idx]}\t{hist_norm[idx]}"
                write.write(line + "\n")
    plt.show()


def mem_cyl(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    num: Annotated[int, typer.Option("-num", help="lipid number")] = 10,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
    hoxy: Annotated[
        str, typer.Option("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "resname DMPC",
    coord_n: Annotated[int, typer.Option("-coord_n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d")] = 1.0,
    coord_r: Annotated[float, typer.Option("-coord_r")] = 8.0,
    coord_z: Annotated[float, typer.Option("-coord_z")] = 0.75,
    padding: Annotated[float, typer.Option("-padding")] = 0.5,
    coord_h: Annotated[float, typer.Option("-coord_h")] = 0.25,
    color: Annotated[str, typer.Option("-color")] = "C0",
):
    """Plot closest lipids to MEMCOM"""
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    import matplotlib.pyplot as plt
    from dztools.misc.mem_help import calc_chain

    u = mda.Universe(top, xtc)
    lipid = u.select_atoms("name P")
    tpi = 2 * np.pi

    x = []
    nsp_h = np.zeros(26)
    nsp_h2 = np.zeros(26)
    nsp_h3 = np.zeros(26)
    xcyls, ycyls = [], []
    xcyls2, ycyls2 = [], []
    xcyls3, ycyls3 = [], []
    eps, ds, ops = [], [], []
    for idx, ts in enumerate(u.trajectory):
        box = ts.dimensions
        z_mem = lipid.atoms.center_of_mass()[-1]

        epsilon, _, x_cyl, y_cyl, nsp = calc_chain(u, lip=lip)
        _, _, x_cyl2, y_cyl2, nsp2 = calc_chain(u, lip=lip, hoxy="O11 O12 O13 O14")
        _, _, x_cyl3, y_cyl3, nsp3 = calc_chain(u, lip=lip, hoxy="OH2")
        xcyls.append(x_cyl)
        ycyls.append(y_cyl)
        xcyls2.append(x_cyl2)
        ycyls2.append(y_cyl2)
        xcyls3.append(x_cyl3)
        ycyls3.append(y_cyl3)

        # flipflop 2
        lip_z = lipid.atoms.positions[:, 2]
        z_idxes = np.argsort(lip_z)
        dz_list = []
        for ii in range(1):
            z_hmax = z_idxes[len(z_idxes)//2 +ii]
            z_hmin = z_idxes[len(z_idxes)//2-(1+ii)]
            dz_list.append(np.abs(lip_z[z_hmin] - lip_z[z_hmax]))

        b = 0.01
        op = epsilon - b*dz_list[0]
        if op > 0.2550:
            continue

        eps.append(epsilon)
        ds.append(dz_list[0])
        ops.append(op)

        x.append(idx)
        nsp_h += nsp
        nsp_h2 += nsp2
        nsp_h3 += nsp3


    if out:
        if 0 in [np.sum(nsp_h), np.sum(nsp_h2), np.sum(nsp_h3)]:
            return False
        with open(out, 'w') as write:
            line = f"# recalc, samples: {len(x)}\n"
            line += f"# avg_eps {np.average(eps)}\n"
            line += f"# avg_d1 {np.average(ds)}\n"
            line += f"# avg_ops {np.average(ops)}"
            write.write(line + "\n")
            for idx in range(len(nsp_h)):
                line = f"{idx}\t{nsp_h[idx]/len(x):.04f}\t{nsp_h2[idx]/len(x):.04f}\t{nsp_h3[idx]/len(x):.04f}"
                #line += f"\t{box[0]:.04f}\t{xcyl[idx]:.04f}\t{y_cyl:.04f}\t{x_cyl2:.04f}\t{y_cyl2:.04f}\t{x_cyl3:.04f}\t{y_cyl3:.04f}"
                write.write(line + "\n")
    # plt.show()
    return True


def mem_mic(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Plot how flat lipids are"""
    import MDAnalysis as mda
    from MDAnalysis.analysis.distances import distance_array
    from MDAnalysis.lib.pkdtree import PeriodicKDTree
    import numpy as np
    from dztools.misc.mem_help import f_axial, f_radial, psi_switch
    from dztools.misc.mem_help import pcom_axis
    import matplotlib.pyplot as plt

    print(top, xtc)
    u = mda.Universe(top, xtc)
    plip = u.select_atoms("name P and resid 1 to 178")
    ball = u.select_atoms("resid 1 to 50 and name P")
    memb = u.select_atoms("resid 51 to 306 and name P")

    x = []
    ops1 = []
    ops2 = []
    ops3 = []
    ops4 = []
    for idx, ts in enumerate(u.trajectory):
        x.append(idx)
        box = ts.dimensions[:3]
        tree = PeriodicKDTree(box=ts.dimensions)
        mcom = memb.center_of_mass()
        bcom = ball.center_of_mass()

        zcom = np.array([0, 0, mcom[2]])
        zball = ball.atoms.positions[:, 2]

        ballp = ball.atoms.positions.copy()
        ballp[:, 0] *= 0
        ballp[:, 1] *= 0

        boxx = list(box[:3]) + [90.]*3
        da = distance_array(zcom, ballp, box=boxx)[0]
        da_avg = np.average(da)

        bcom_z = pcom_axis(ball.positions[:, 2], box[2])
        lz, hz = min((mcom[2], bcom_z)), max((mcom[2], bcom_z))
        half = u.select_atoms(f"name P and prop z >= {lz} and prop z <= {hz}")
        hpos = half.atoms.positions.copy()
        hpos[:, 2] *= 0
        tree.set_coords(hpos, cutoff=20) # pbc cutoff or something

        samples = 50
        xlim = np.linspace(0, box[0], samples)
        X, Y = np.meshgrid(xlim, xlim)  # Create 2D grids
        w = np.array([X.flatten(), Y.flatten(), np.zeros(samples**2)]).T
        hehe = tree.search_tree(w, radius=8.)
        all_indices = np.arange(len(X.flatten()))
        idxes = np.setdiff1d(all_indices, hehe[:, 0])
        print(idx, len(idxes))
        # print("nani", idxes, len(idxes))
        # if len(idxes) > 0:
        #     bcom_x = pcom_axis(w[:, 0][idxes], box[0])
        #     bcom_y = pcom_axis(w[:, 1][idxes], box[1])
        #     pcom = np.array([bcom_x, bcom_y, 0])

        #     daz = distance_array(pcom, hpos, box=ts.dimensions)[0]
        #     asort = np.argsort(daz)
        #     plt.scatter(hpos[:, 0][asort[:5]], hpos[:, 1][asort[:5]], alpha=0.8, s=100)
        #     op = np.sum(sorted(daz)[:5])

        #     print("uhuh", daz)
        #     # exit()

        #     plt.scatter(X.flatten()[idxes], Y.flatten()[idxes], marker="x")
        #     plt.scatter(bcom_x, bcom_y, s=50)
        # else:
        #     op = 0
        # plt.scatter(hpos[:, 0], hpos[:, 1])
        # plt.title(f"{len(idxes)}")
        # plt.show()
        # exit()

        # bcom2_x = pcom_axis(half.positions[:, 0], box[0])
        # bcom2_y = pcom_axis(half.positions[:, 1], box[1])
        # bcom2_z = pcom_axis(half.positions[:, 2], box[2])
        # # half = u.select_atoms(f"name P and resid 1 to 50")
        # print("len", len(half.atoms))
        # half_posz = half.atoms.positions.copy()
        # half_posz[:, 0] *= 0
        # half_posz[:, 1] *= 0
        # # daz = distance_array(half_posz, half, box=boxx)[0]
        # daz = distance_array(zcom, half_posz, box=boxx)[0]
        # asort = np.argsort(daz)
        # zsum = 0
        # bcom2 = np.array([bcom2_x, bcom2_y, 0])

        # # print(asort)
        # for at_idx in asort[:10]:
        #     bpos = half.atoms[at_idx].position.copy()
        #     ax.scatter(
        #         bpos[ 0],
        #         bpos[ 1],
        #         bpos[ 2],
        #     color="C2", s=100)
        #     bpos[2] = 0
        #     atd = distance_array(bpos, bcom2, box=boxx)[0][0]
        #     print('sn', atd)
        #     zsum += atd
        # ops4.append(zsum)
        # print("snow", zsum)

        # ax.scatter(
        #     half.atoms.positions[:, 0],
        #     half.atoms.positions[:, 1],
        #     half.atoms.positions[:, 2],
        # color="C3", s=50)

        # ax.scatter(
        #     bcom2_x,
        #     bcom2_y,
        #     bcom2_z,
        #     color="C1", s=50)
        # print("oprd", zsum)
        # plt.title(f"{zsum}")
        # plt.show()
        # exit()

        # half = u.select_atoms(f"name P and resid 1 to 50 and prop z >= {lz} and prop z <= {hz}")
        # print("heiH", len(half.atoms))
        # exit()
        # half_pos = half.atoms.positions
        # # half_pos = ball.atoms.positions
        # half_pos[:, 2] *= 0
        # zda = sorted(distance_array(np.array([bcom_x, bcom_y, 0]), half_pos, box=boxx)[0])[:5]
        # ops4.append(sum(zda))
        # # print(" ".join([f"{i:.04f}" for i in zda]))
        # print("ha?", sum(zda))


        # half = u.select_atoms(f"resid 1 to 50 and name P and prop z <= {hz}")
        # exit()

        # cyl
        daball = sorted(distance_array(bcom, ballp, box=boxx)[0])
        # bdadel = max(daball) - min(daball)
        bdadel = np.average(daball[-5:]) - np.average(daball[:5])

        ops1.append(da_avg)
        ops2.append(bdadel)
        ops3.append(sum(daball[-5:]))
        ops4.append(len(idxes))

    if out:
        with open(out, 'w') as write:
            for idx in x:
                line = f"{idx}\t{ops1[idx]}\t{ops2[idx]}\t{ops3[idx]}\t{ops4[idx]}"
                write.write(line + "\n")


def mem_pfcvs(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    hoxy: Annotated[
        str, typer.Option("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "resname DMPC",
    coord_n: Annotated[int, typer.Option("-coord_n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d")] = 1.0,
    coord_r: Annotated[float, typer.Option("-coord_r")] = 8.0,
    coord_z: Annotated[float, typer.Option("-coord_z")] = 0.75,
    padding: Annotated[float, typer.Option("-padding")] = 0.5,
    coord_h: Annotated[float, typer.Option("-coord_h")] = 0.25,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    from dztools.misc.mem_help import calc_chain
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.analysis.distances import distance_array

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    lipid = u.select_atoms(f"{lip}")
    lipid_p = u.select_atoms(f"name P")
    tpi = 2 * np.pi

    totlen = len(u.trajectory)
    eps_ch, eps_p = [], []
    cylx, cyly = [], []
    d1, d2, d3 = [], [], []
    d12, d123 = [], []
    dph, symff = [], []
    ff3_1 = []
    ff3_2 = []
    ff3_3 = []
    ff3_4 = []
    ff3_5 = []
    ff3_12 = []
    ff3_123 = []
    ff3_1234 = []
    ff3_12345 = []

    lip_z = lipid_p.atoms.positions[:, 2]
    z_idxes0 = np.argsort(lip_z)[:64]

    coms1 = []
    coms2 = []
    coms3 = []
    coms12 = []
    coms123 = []

    for idx, ts in enumerate(u.trajectory):
        # Frame properties
        epsilon, epsilon_e, x0, y0, _ = calc_chain(u, lip=lip, coord_r = coord_r, coord_n=coord_n)
        eps_ch.append(epsilon)
        eps_p.append(epsilon_e)
        cylx.append(x0)
        cyly.append(y0)

        # flipflop 2
        lip_z = lipid_p.atoms.positions[:, 2]
        z_mem = lipid.atoms.center_of_mass()[-1]
        z_idxes = np.argsort(lip_z)
        dz_list = []
        for ii in range(3):
            z_hmax = z_idxes[len(z_idxes)//2 +ii]
            z_hmin = z_idxes[len(z_idxes)//2-(1+ii)]
            dz_list.append(np.abs(lip_z[z_hmin] - lip_z[z_hmax]))

        d1.append(dz_list[0])
        d2.append(dz_list[1])
        d3.append(dz_list[2])
        d12.append(dz_list[0] + dz_list[1])
        d123.append(dz_list[0] + dz_list[1] + dz_list[2])
        dph.append(np.min(np.abs(lip_z - z_mem)))

        symff.append(64 - len([i for i in z_idxes[:64] if i in z_idxes0]))

        lip_xyz = lipid_p.atoms.positions
        lipnum = len(z_idxes)
        idx_up = z_idxes[:lipnum//2]
        idx_dw = z_idxes[lipnum//2:]
        boxx = np.copy(ts.dimensions)
        boxx[2] = 1000
        dists = distance_array(lip_xyz[idx_up], lip_xyz[idx_dw], box=boxx)
        indices = np.dstack(np.unravel_index(np.arange(dists.size), dists.shape))[0]
        values = dists.flatten()
        sorted_pairs = sorted(zip(values, indices), key=lambda x: x[0])

        fdw, fup, ffds = [], [], []
        N = 5
        for val, (x, y) in sorted_pairs:
            if idx_up[x] not in fup and idx_dw[y] not in fdw:
                fup.append(idx_up[x])
                fdw.append(idx_dw[y])
                ffds.append(val)
            if len(ffds) >= N:
                break
        ff3_1.append(np.abs(lip_z[fdw[0]] - lip_z[fup[0]]))
        ff3_2.append(np.abs(lip_z[fdw[1]] - lip_z[fup[1]]))
        ff3_3.append(np.abs(lip_z[fdw[2]] - lip_z[fup[2]]))
        ff3_4.append(np.abs(lip_z[fdw[3]] - lip_z[fup[3]]))
        ff3_5.append(np.abs(lip_z[fdw[4]] - lip_z[fup[4]]))
        ff3_12.append(ff3_1[-1] + ff3_2[-1])
        ff3_123.append(ff3_1[-1] + ff3_2[-1] + ff3_3[-1])
        ff3_1234.append(ff3_1[-1] + ff3_2[-1] + ff3_3[-1] + ff3_4[-1])
        ff3_12345.append(ff3_1[-1] + ff3_2[-1] + ff3_3[-1] + ff3_4[-1] + ff3_5[-1])

        # print("f", fdw, fup)
        resid1d = u.select_atoms(f"resid {fdw[0]+1}").center_of_mass(unwrap=True)
        resid1u = u.select_atoms(f"resid {fup[0]+1}").center_of_mass(unwrap=True)
        resid2d = u.select_atoms(f"resid {fdw[1]+1}").center_of_mass(unwrap=True)
        resid2u = u.select_atoms(f"resid {fup[1]+1}").center_of_mass(unwrap=True)
        resid3d = u.select_atoms(f"resid {fdw[2]+1}").center_of_mass(unwrap=True)
        resid3u = u.select_atoms(f"resid {fup[2]+1}").center_of_mass(unwrap=True)

        # print("bonk a", resid1d, lip_z[fdw[0]])
        # print("bonk b", resid1u, lip_z[fup[0]])
        coms1.append(np.abs(resid1d[2] - resid1u[2]))
        # print(resid1d[2], resid1u[2])
        # print(coms1[-1], ff3_1[-1])
        # print("mouse")
        # exit()
        coms2.append(np.abs(resid2d[2] - resid2u[2]))
        coms3.append(np.abs(resid3d[2] - resid3u[2]))
        coms12.append(coms1[-1] + coms2[-1])
        coms123.append(coms1[-1] + coms2[-1] + coms3[-1])


    if out:
        with open(out, 'w') as write:
            for idx in range(len(u.trajectory)):
                string = f"{idx}\t{eps_ch[idx]:.08f}\t{eps_p[idx]:.08f}\t{cylx[idx]:.08f}\t{cylx[idx]:.08f}"                        # 0 1 2 3 4
                string += f"\t{d1[idx]:.08f}\t{d2[idx]:.08f}\t{d3[idx]:.08f}"                                                       # 5 6 7
                string += f"\t{d12[idx]:.08f}\t{d123[idx]:.08f}"                                                                    # 8 9
                string += f"\t{dph[idx]:.08f}\t{symff[idx]:.08f}"                                                                   # 10 11
                string += f"\t{ff3_1[idx]:.08f}\t{ff3_2[idx]:.08f}\t{ff3_3[idx]:.08f}\t{ff3_4[idx]:.08f}\t{ff3_5[idx]:.08f}"        # 12 13 14 15 16
                string += f"\t{ff3_12[idx]:.08f}\t{ff3_123[idx]:.08f}\t{ff3_1234[idx]:.08f}\t{ff3_12345[idx]:.08f}"                 # 17 18 19 20
                string += f"\t{coms1[idx]:.08f}\t{coms2[idx]:.08f}\t{coms3[idx]:.08f}\t{coms12[idx]:.08f}\t{coms123[idx]:.08f}\n"   # 21 22 23 24 25
                write.write(string)


def mem_thin(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    hoxy: Annotated[
        str, typer.Option("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "resname DMPC",
    coord_n: Annotated[int, typer.Option("-coord_n")] = 26,
    coord_d: Annotated[float, typer.Option("-coord_d")] = 1.0,
    coord_r: Annotated[float, typer.Option("-coord_r")] = 8.0,
    coord_z: Annotated[float, typer.Option("-coord_z")] = 0.75,
    padding: Annotated[float, typer.Option("-padding")] = 0.5,
    coord_h: Annotated[float, typer.Option("-coord_h")] = 0.25,
    lmt_n: Annotated[int, typer.Option("-lmt_n")] = 12,
    lmt_k: Annotated[int, typer.Option("-lmt_k")] = 30,
    plot: Annotated[bool, typer.Option("-plot", help="plot")] = False,
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    from dztools.misc.mem_help import calc_chain, calc_thin
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.analysis.distances import distance_array
    from MDAnalysis.lib.distances import calc_angles

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    lipid = u.select_atoms(f"{lip}")
    lipid_p = u.select_atoms(f"name P")
    if len(lipid_p.atoms) == 0:
        lipid_p = u.select_atoms(f"name PH")
    tpi = 2 * np.pi

    totlen = len(u.trajectory)
    eps_ch, eps_p = [], []
    eps_ch2, eps_p2 = [], []
    mint = []
    track = []

    lip_z = lipid_p.atoms.positions[:, 2]
    z_idxes0 = np.argsort(lip_z)[:64]

    for idx, ts in enumerate(u.trajectory):
        epsilon, epsilon_e, x0, y0, dlmt = calc_thin(u,
                                                     lip=lip,
                                                     coord_r = coord_r,
                                                     coord_n=coord_n,
                                                     hoxy=hoxy,
                                                     coord_d=coord_d,
                                                     coord_z=coord_z,
                                                     lmt_n=lmt_n,
                                                     lmt_k=lmt_k,)

        epsilon2, epsilon_e2, _, _, _ = calc_chain(u,
                                                   lip=lip,
                                                   coord_r = coord_r,
                                                   coord_n=coord_n,
                                                   hoxy="OH2",
                                                   coord_d=coord_d,
                                                   coord_z=coord_z)
        eps_ch.append(epsilon)
        eps_p.append(epsilon_e)
        eps_ch2.append(epsilon2)
        eps_p2.append(epsilon_e2)
        mint.append(dlmt)

        zcom = lipid.center_of_mass()[2]
        p_z = lipid_p.atoms.positions[:, 2]
        p_zu = p_z > zcom
        track.append(np.abs(64-np.sum(p_zu)))

    if out:
        with open(out, 'w') as write:
            min0 = np.min([len(u.trajectory), len(eps_ch)])
            for idx in range(min0):
                string = f"{idx}\t{eps_ch[idx]:.08f}\t{eps_p[idx]:.08f}"                    # 0 1 2 3 4
                string += f"\t{mint[idx]:.08f}"                                                   # 5 6 7
                string += f"\t{eps_ch2[idx]:.08f}\t{eps_p2[idx]:.08f}"                    # 0 1 2 3 4
                string += f"\t{track[idx]:.08f}\n"                                                   # 5 6 7
                write.write(string)


def mem_void(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    rad: Annotated[float, typer.Option("-coord_z")] = 6.00,
    samples: Annotated[int, typer.Option("-coord_z")] = 100,
    # lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "resname DMPC",
    out: Annotated[str, typer.Option("-out", help="string")] = "",
):
    """stick ball"""

    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.analysis.distances import distance_array
    from MDAnalysis.lib.pkdtree import PeriodicKDTree

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    lipid_p = u.select_atoms(f"name P")
    if len(lipid_p.atoms) == 0:
        lipid_p = u.select_atoms(f"name PH")
    tpi = 2 * np.pi

    totlen = len(u.trajectory)

    cnts_k = []
    cnts_ku = []
    cnts_kl = []
    for idx, ts in enumerate(u.trajectory):
        # get standard info from frame
        pos = lipid_p.atoms.positions
        boxx = ts.dimensions[0]
        zavg = np.average(pos[:, 2])
        upper = lipid_p.select_atoms(f"prop z > {zavg}")
        lower = lipid_p.select_atoms(f"prop z <= {zavg}")

        # project P into x, y
        pos_u = np.copy(upper.atoms.positions * [1, 1, 0])
        pos_l = np.copy(lower.atoms.positions * [1, 1, 0])

        # create 2D xy dotted grid
        xlim = np.linspace(0, boxx, samples)
        X, Y = np.meshgrid(xlim, xlim)  # Create 2D grids
        all_indices = np.arange(len(X.flatten()))
        w = np.array([X.flatten(), Y.flatten(), np.zeros(samples**2)]).T

        cnts = []

        for pos0 in [pos_u, pos_l]:
            tree = PeriodicKDTree(box=ts.dimensions)
            tree.set_coords(pos0, cutoff=40) # pbc cutoff or something
            dots = tree.search_tree(w, radius=rad)
            idxes = np.setdiff1d(all_indices, dots[:, 0])
            cnts.append(len(X.flatten()[idxes]))
        cnts_k.append(np.average(cnts))
        cnts_ku.append(cnts[0])
        cnts_kl.append(cnts[1])
        print(idx, totlen, cnts_k[-1])

    if out:
        with open(out, 'w') as write:
            # min0 = np.min([len(u.trajectory), len(eps_ch)])
            for idx in range(len(cnts_k)):
                string = f"{idx}\t{cnts_k[idx]:.08f}\t{cnts_ku[idx]:.08f}\t{cnts_kl[idx]:.08f}\n"
                write.write(string)
