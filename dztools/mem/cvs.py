from typing import Annotated as And

from typer import Option as Opt

def mem_helicity(
    top: And[str, Opt("-top", help="gro/pdb/tpr file")],
    xtc: And[str, Opt("-xtc", help="xtc file")] = None,
    lip: And[str, Opt("-lip", help="lipid type")] = "POPC",
    plot: And[bool, Opt("-plot", help="plot")] = False,
    out: And[str, Opt("-out", help="string")] = "",
):
    """DSSP"""

    import matplotlib.pyplot as plt
    import MDAnalysis as mda
    import numpy as np

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

    if out:
        with open(out, 'w') as write:
            for idx, cv in zip(np.arange(len(hels_avg)), hels_avg):
                write.write(f"{idx}\t{cv*100:.08f}\t{cv:.08f}\n")

    return hels_avg


def mem_chain(
    top: And[str, Opt("-top", help="gro/pdb/tpr file")],
    xtc: And[str, Opt("-xtc", help="xtc file")],
    hoxy: And[
        str, Opt("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: And[str, Opt("-lip", help="xtc file")] = "resname POPC",
    coord_n: And[int, Opt("-coord_n")] = 26,
    coord_d: And[float, Opt("-coord_d")] = 1.0,
    coord_r: And[float, Opt("-coord_r")] = 8.0,
    coord_z: And[float, Opt("-coord_z")] = 0.75,
    padding: And[float, Opt("-padding")] = 0.5,
    coord_h: And[float, Opt("-coord_h")] = 0.25,
    plot: And[bool, Opt("-plot", help="plot")] = False,
    out: And[str, Opt("-out", help="string")] = "",
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    import matplotlib.pyplot as plt
    import MDAnalysis as mda
    import numpy as np

    from dztools.misc.mem_help import calc_chain

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    u.select_atoms(f"{lip}")
    epsilons = []
    epsilons_p = []
    2 * np.pi

    len(u.trajectory)
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


def mem_pfcvs(
    top: And[str, Opt("-top", help="gro/pdb/tpr file")],
    xtc: And[str, Opt("-xtc", help="xtc file")],
    hoxy: And[
        str, Opt("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: And[str, Opt("-lip", help="xtc file")] = "resname DMPC",
    coord_n: And[int, Opt("-coord_n")] = 26,
    coord_d: And[float, Opt("-coord_d")] = 1.0,
    coord_r: And[float, Opt("-coord_r")] = 8.0,
    coord_z: And[float, Opt("-coord_z")] = 0.75,
    padding: And[float, Opt("-padding")] = 0.5,
    coord_h: And[float, Opt("-coord_h")] = 0.25,
    plot: And[bool, Opt("-plot", help="plot")] = False,
    out: And[str, Opt("-out", help="string")] = "",
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.analysis.distances import distance_array

    from dztools.misc.mem_help import calc_chain

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    lipid = u.select_atoms(f"{lip}")
    lipid_p = u.select_atoms("name P")
    2 * np.pi

    len(u.trajectory)
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

        resid1d = u.select_atoms(f"resid {fdw[0]+1}").center_of_mass(unwrap=True)
        resid1u = u.select_atoms(f"resid {fup[0]+1}").center_of_mass(unwrap=True)
        resid2d = u.select_atoms(f"resid {fdw[1]+1}").center_of_mass(unwrap=True)
        resid2u = u.select_atoms(f"resid {fup[1]+1}").center_of_mass(unwrap=True)
        resid3d = u.select_atoms(f"resid {fdw[2]+1}").center_of_mass(unwrap=True)
        resid3u = u.select_atoms(f"resid {fup[2]+1}").center_of_mass(unwrap=True)

        coms1.append(np.abs(resid1d[2] - resid1u[2]))
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
    top: And[str, Opt("-top", help="gro/pdb/tpr file")],
    xtc: And[str, Opt("-xtc", help="xtc file")],
    hoxy: And[
        str, Opt("-hoxy", help="xtc file")
    ] = "OH2 O11 O12 O13 O14",
    lip: And[str, Opt("-lip", help="xtc file")] = "resname DMPC",
    coord_n: And[int, Opt("-coord_n")] = 26,
    coord_d: And[float, Opt("-coord_d")] = 1.0,
    coord_r: And[float, Opt("-coord_r")] = 8.0,
    coord_z: And[float, Opt("-coord_z")] = 0.75,
    padding: And[float, Opt("-padding")] = 0.5,
    coord_h: And[float, Opt("-coord_h")] = 0.25,
    lmt_n: And[int, Opt("-lmt_n")] = 12,
    lmt_k: And[int, Opt("-lmt_k")] = 30,
    plot: And[bool, Opt("-plot", help="plot")] = False,
    out: And[str, Opt("-out", help="string")] = "",
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    import MDAnalysis as mda
    import numpy as np

    from dztools.misc.mem_help import calc_chain, calc_thin

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    lipid = u.select_atoms(f"{lip}")
    lipid_p = u.select_atoms("name P")
    if len(lipid_p.atoms) == 0:
        lipid_p = u.select_atoms("name PH")
    2 * np.pi

    len(u.trajectory)
    eps_ch, eps_p = [], []
    eps_ch2, eps_p2 = [], []
    mint = []
    track = []
    lavgs = []
    stdavgs = []

    lip_z = lipid_p.atoms.positions[:, 2]
    np.argsort(lip_z)[:64]

    for idx, ts in enumerate(u.trajectory):
        epsilon, epsilon_e, x0, y0, dlmt, lavg, stdavg = calc_thin(u,
                                                                    lip=lip,
                                                                    coord_r = coord_r,
                                                                    coord_n=coord_n,
                                                                    hoxy=hoxy,
                                                                    coord_d=coord_d,
                                                                    coord_z=coord_z,
                                                                    lmt_n=lmt_n,
                                                                    lmt_k=lmt_k,)

        # epsilon2, epsilon_e2, _, _, _ = calc_chain(u,
        #                                            lip=lip,
        #                                            coord_r = coord_r,
        #                                            coord_n=coord_n,
        #                                            hoxy="OH2",
        #                                            coord_d=coord_d,
        #                                            coord_z=coord_z)
        eps_ch.append(epsilon)
        eps_p.append(epsilon_e)
        # eps_ch2.append(epsilon2)
        # eps_p2.append(epsilon_e2)
        mint.append(dlmt)
        lavgs.append(lavg)
        stdavgs.append(stdavg)

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
                # string += f"\t{eps_ch2[idx]:.08f}\t{eps_p2[idx]:.08f}"                    # 0 1 2 3 4
                string += f"\t{track[idx]:.08f}\t{lavgs[idx]:.08f}\t{stdavgs[idx]:.08f}\n"                     # 5 6 7
                write.write(string)


def mem_void(
    top: And[str, Opt("-top", help="gro/pdb/tpr file")],
    xtc: And[str, Opt("-xtc", help="xtc file")],
    rad: And[float, Opt("-coord_z")] = 6.00,
    samples: And[int, Opt("-coord_z")] = 100,
    # lip: And[str, Opt("-lip", help="xtc file")] = "resname DMPC",
    out: And[str, Opt("-out", help="string")] = "",
):
    """stick ball"""

    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.lib.pkdtree import PeriodicKDTree

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)

    lipid_p = u.select_atoms("name P")
    if len(lipid_p.atoms) == 0:
        lipid_p = u.select_atoms("name PH")
    2 * np.pi

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


def mem_lipodiso(
    top: And[str, Opt("-top", help="gro/pdb/tpr file")],
    xtc: And[str, Opt("-xtc", help="xtc file")],
    com: And[str, Opt("-com", help="C6")],
    out: And[str, Opt("-out", help="string")] = "lipodiso.txt",
):
    """farthest clostest neighbor

                           tot time                 time per
    just neighbor: 1.102856159210205 0.0026768353378888473
    neigh and double: 0m1.698s


    """

    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.analysis.distances import distance_array

    import time

    u = mda.Universe(top, xtc)
    resids = set([str(i.resid) for i in u.atoms if i.resname == "SDS"])
    coms = u.select_atoms(f"name {com}")
    assert len(coms) > 0

    results = []
    time0 = time.time()
    for idx, ts in enumerate(u.trajectory):
        da1 = distance_array(coms, coms, box=ts.dimensions)

        # first, second, third minimal neighbour
        sc = np.partition(da1, 1, axis=1)[:, 1]
        sort = np.argsort(sc)[-3:,][::-1]
        ops = sc[sort]

        resid1 = u.select_atoms(f"resid {coms[sort[0]].resid}")
        resid2 = u.select_atoms(f"resid " + " ".join(resids - {str(coms[sort[0]].resid)}))

        da2 = distance_array(resid1, resid2, box=ts.dimensions)
        minda2 = np.min(da2)

        results.append([idx, ops[0], ops[1], ops[2], sort[0], sort[1], sort[2], minda2])
    time1 = time.time()

    print("wut", idx, len(u.trajectory), time1-time0, (time1-time0)/(idx+1))

    np.savetxt(out, results)


def mem_lipodiso2(
    top: And[str, Opt("-top", help="gro/pdb/tpr file")],
    xtc: And[str, Opt("-xtc", help="xtc file")],
    out: And[str, Opt("-out", help="string")] = "lipodiso.txt",
    lim: And[int, Opt("-lim", help="string")] = 10,
):
    """farthest clostest neighbor

                           tot time                 time per
    just neighbor: 1.102856159210205 0.0026768353378888473
    neigh and double: 0m1.698s


    """

    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.analysis.distances import distance_array
    from dztools.misc.mem_help import pcom_axis

    import time

    u = mda.Universe(top, xtc)
    pots = u.select_atoms(f"name POT")
    carbs = u.select_atoms(f"name C*")[::2]
    lips = len(u.select_atoms(f"name S"))
    lip_ats = int(len(carbs)/lips + 1)

    carbs_name = " ".join(["S"] + [i.name for i in carbs])
    group = u.select_atoms(f"name {carbs_name}")
    resids = group.resids
    mask = resids[:, None] == resids[None, :]

    results = []
    for idx, ts in enumerate(u.trajectory):
        result = [idx]
        da = distance_array(group, group, box=ts.dimensions)

        # disregard interlipid distances
        da[mask] += 1000

        # calc percentage
        sc = np.partition(da, 1, axis=1)[:, 0]
        sc2 = sc.copy()
        sc2r = sc2.reshape(lips, lip_ats)
        sc[sc>lim] = lim
        # avg = sc.reshape(lips, lip_ats).mean(axis=1)/lim
        avg = sc.reshape(lips, lip_ats).mean(axis=1)/lim

        sort = np.argsort(avg)[-3:,][::-1]
        ops = avg[sort]

        disso_idx = sort[0]
        disso_avrg = np.average(sc2r[disso_idx])
        disso_head = sc2r[disso_idx][0]
        disso_tail = sc2r[disso_idx][-1]

        # mask = np.zeros(len(group), dtype=bool)
        # start = disso_idx * lip_ats
        # stop  = (disso_idx + 1) * lip_ats
        # mask[start:stop] = True

        # group2 = group[mask]      # the dissociating lipid
        # group3 = group[~mask]     # everything else

        group2 = group[disso_idx*lip_ats:(disso_idx+1)*lip_ats]
        da2 = distance_array(group2, pots, box=ts.dimensions)
        minpot = np.min(da2)
        # sc2 = np.partition(da2, 1, axis=1)[:, 0]
        # minpot = np.min(sc2)

        # resid = resids[disso_idx*lip_ats]
        # print(sc2r[disso_idx])
        # print(np.min(da2, axis=1))
        # print(sc2)
        # exit()

        result += list(ops) + list(sort) + [disso_avrg, disso_head, disso_tail, minpot]
        print(idx, ops[0], sort[0])
        results.append(result)

    np.savetxt(out, results)
