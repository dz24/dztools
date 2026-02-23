def calc_met_com(mdau, lip, num_resi):
    """Calculates individual and average MELITTIN COMS wrt. specified lipid
    membrane."""
    import numpy as np
    import MDAnalysis as mda

    # Select individual proteins and membrane
    lipid = mdau.select_atoms(f"resname {lip}")
    protein = mdau.select_atoms("protein")
    num_prot = len(protein.residues)//num_resi

    # define coms and iterate over all xtc frames
    lcoms_z, pcoms_z_avg = [], []
    pcoms_z = [[] for _ in range(num_prot)]
    for idx, ts in enumerate(mdau.trajectory):
        pcoms_z_avg.append(0)
        lcoms_z.append(lipid.atoms.center_of_mass()[-1])

        for i in range(num_prot):
            mel = mdau.select_atoms(f"resid {1+num_resi*i}:{num_resi*(i+1)}")
            pcoms_z[i].append(mel.center_of_mass(unwrap=True)[-1])
            pcoms_z_avg[-1] += abs(pcoms_z[i][-1] - lcoms_z[-1])
        pcoms_z_avg[-1] /= num_prot

    return abs(np.array(pcoms_z) - lcoms_z), np.abs(np.array(pcoms_z_avg))


def calc_helicity(mdau, num_resi):
    """Calculates individual and average MELITTIN helicities using DSSP."""
    import numpy as np
    from MDAnalysis.analysis.dssp import DSSP

    protein = mdau.select_atoms("protein")
    hels = []
    for i in range(int(len(protein.residues) / num_resi)):
        mel = mdau.select_atoms(f"resid {num_resi*i}:{num_resi*(i+1)}")
        hels.append(
            [
                sum(np.array(f) == "H") / num_resi
                for f in DSSP(mel).run().results.dssp
            ]
        )

    return hels, np.average(np.array(hels), axis=0)


def psi_switch(x, zeta):
    import numpy as np
    if x <= 1:
        return x*zeta
    else:
        b = zeta/(1-zeta)
        c = (1-zeta)*np.exp(b)
        return 1 - c*np.exp(-b*x)

def theta(x, h):
    import numpy as np
    xnew = np.zeros(len(x))

    # if
    xnew += (-1 + h <= x) * (x <= 1 - h) * 1

    # elif
    xnew += (1 - h < x) * (x < 1 + h) * (1/2 - (3/(4*h))*(x-1) +
                                         (1/(4*h**3))*(x-1)**3)
    # elif
    xnew += (-1 - h < x) * (x < -1 + h) * (1/2 + (3/(4*h))*(x+1) -
                                           (1/(4*h**3))*(x+1)**3)

    return xnew

def f_axial(zi, zs, ds, h=1/4):
    return theta((zi-zs)/(ds/2),h)

def f_radial(pos, Xcyl, Ycyl, Rcyl, box, h=1/4):
    import numpy as np
    from MDAnalysis.analysis import distances
    pos[:, 2] = 0
    ri = distances.distance_array(np.array([Xcyl, Ycyl, 0]), pos, box=box)[0]
    return theta(ri/Rcyl, h)

def gyr_com(atoms, box, dim):
    import numpy as np

    angle = [2*np.pi*i.position[dim]/box[dim] for i in atoms.atoms]
    cos, sin =  np.cos(angle), np.sin(angle)
    avg_cos, avg_sin = np.average(cos), np.average(sin)
    theta = np.arctan2(-avg_sin, -avg_cos) + np.pi
    com = box[dim]*theta/(2*np.pi)

    return com

def gyr_com2(pos, box, dim):
    import numpy as np

    angle = [2*np.pi*i[dim]/box[dim] for i in pos]
    cos, sin =  np.cos(angle), np.sin(angle)
    avg_cos, avg_sin = np.average(cos), np.average(sin)
    theta = np.arctan2(-avg_sin, -avg_cos) + np.pi
    com = box[dim]*theta/(2*np.pi)

    return com


def calc_chain(
    u,
    hoxy = "OH2 O11 O12 O13 O14",
    lip= "resname POPC",
    coord_n = 26,
    coord_d = 1.0,
    coord_r = 8.0,
    coord_z = 0.75,
    padding = 0.5,
    coord_h = 0.25,
    r0 = 4.05,
    d0 = 10.0,
    v0 = 29.96,
    ep0 = 0.925,
):
    """Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

    from dztools.misc.mem_help import f_axial, f_radial, psi_switch
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda

    lipid = u.select_atoms(f"{lip}")
    epsilons = []
    tpi = 2 * np.pi

    z_mem = lipid.atoms.center_of_mass()[-1]
    z_s = z_mem + (np.arange(coord_n) + 1 / 2 - coord_n / 2) * coord_d
    hoxys = u.select_atoms(f"name {hoxy} and prop z < {z_s[-1]+padding} and prop z > {z_s[0]-padding}")
    atoms_x = hoxys.atoms.positions[:, 0]
    atoms_y = hoxys.atoms.positions[:, 1]
    box = u.dimensions

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
    in_radi = f_radial(hoxys.atoms.positions, x_cyl, y_cyl, coord_r, box, coord_h)
    epsilon = 0
    nsp = [0 for i in range(coord_n)]
    for s in range(coord_n):
        nsp[s] = np.sum(in_axis[s] * in_radi)
        epsilon += psi_switch(nsp[s], coord_z)
    epsilon /= coord_n

    # expansion
    exp_num = u.select_atoms(f"name {hoxy} and prop z > {z_mem - d0 - padding} and prop z < {z_mem + d0 + padding}")
    exp_z = exp_num.atoms.positions[:, 2]
    npr = np.sum(theta((exp_z - z_mem)/(d0/2), h=0.1))
    # due to numerical error, npr can become negative when really it is just should become 0.
    if npr < 0:
        npr = 0
    rr = np.sqrt((npr*v0)/(np.pi*d0))

    epsilon_e = epsilon + heavys(np.array([epsilon - ep0]), 0.05)*((rr-r0)/r0)
    # print(epsilon, epsilon_e, heavys(np.array([epsilon - ep0]), 0.05)*((rr-r0)/r0))

    return epsilon, epsilon_e[0], x_cyl, y_cyl, nsp


def heavys(x, eps):
    import numpy as np
    xnew = np.zeros(len(x))

    # elif
    xnew += (x <= eps) * (-eps <= x) * (1/2 + ((3*x)/(4*eps)) - (x**3)/(4*eps**3))

    # else 
    xnew += (x > eps) * 1

    return xnew


def pcom_axis(pos, box, mass=None):
    import numpy as np

    if mass is None:
        mass = np.ones_like(pos)

    tpi = np.pi * 2
    ang_s = np.sin(tpi * pos / box)
    ang_c = np.cos(tpi * pos / box)

    avg_s = np.average(ang_s, weights=mass)
    avg_c = np.average(ang_c, weights=mass)

    com = (np.arctan2(-avg_s, - avg_c) + np.pi) * box / tpi
    return com

def periodic_com(pos, box):
    com_x = pcom_axis(pos[:, 0], box[0])
    com_y = pcom_axis(pos[:, 1], box[1])
    return com_x, com_y


def calc_thin(
    u,
    hoxy = "OH2 O11 O12 O13 O14",
    lip= "resname POPC",
    coord_n = 26,
    coord_d = 1.0,
    coord_r = 8.0,
    coord_z = 0.75,
    padding = 0.5,
    coord_h = 0.25,
    r0 = 4.05,
    d0 = 10.0,
    v0 = 29.96,
    ep0 = 0.925,
    lmt_n = 12,
    lmt_k = 30,
):
    from MDAnalysis.analysis.distances import distance_array
    from MDAnalysis.lib.distances import calc_angles
    import numpy as np

    epsilon, epsilon_e, x0, y0, _ = calc_chain(u, lip=lip, coord_r = coord_r, coord_n=coord_n, hoxy=hoxy, coord_d=coord_d, coord_z=coord_z)

    lipid = u.select_atoms(f"{lip}")
    lipid_p = u.select_atoms(f"name P")
    if len(lipid_p.atoms) == 0:
        lipid_p = u.select_atoms(f"name PH")


    xyz = lipid_p.atoms.positions
    zavg = np.average(xyz[:, 2])
    radial = distance_array(np.array([x0, y0, 0]),
                            xyz*np.array([1, 1, 0]),
                            box = u.dimensions)[0]

    ups, dws = [], []           # group idxes
    ups_c, dws_c = [], []       # angles
    for i in np.argsort(radial):
        at = lipid_p.atoms[i]
        resid = u.select_atoms(f"resid {at.resid}")
        com = resid.center_of_mass(unwrap=True)
        coords3 = np.array([xyz[i, 0], xyz[i, 1], com[2]])
        ang = calc_angles(xyz[i], com, coords3, box=u.dimensions)
        if xyz[i, 2] > zavg:
            ups.append(i)
            ups_c.append(ang)
        else:
            dws.append(i)
            dws_c.append(ang)
        if len(ups) >= 1 and len(dws) >= 1 and len(ups) + len(dws) >= lmt_n:
            break
    radial1 = radial[i]

    # get their positions, and get their pair stuff
    ups_pos = xyz[ups]
    dws_pos = xyz[dws]

    # standard average
    stdavg = np.average(ups_pos[:, 2]) - np.average(dws_pos[:, 2])

    ip = np.array([(i, j) for i in range(len(ups)) for j in range(len(dws))]).T
    pdist = distance_array(ups_pos, dws_pos, box = u.dimensions * [1,1,10,1,1,1])
    dzs = []
    coms = []

    # sorted pdist indxes
    ind = np.unravel_index(np.argsort(pdist, axis=None), pdist.shape)
    for cnt, (i, j) in enumerate(zip(ind[0], ind[1])):
        dzs.append(ups_pos[i, 2] - dws_pos[j, 2])
        coms.append((ups_c[i] + dws_c[j])/2)
        if cnt > lmt_k:
            break

    lavg = np.average(dzs)
    # print(f"{np.abs(stdavg- lavg):.04f}", stdavg, lavg)
    d_lmt = lavg*np.average(coms)*(2/np.pi)
    return epsilon, epsilon_e, x0, y0, d_lmt, lavg, stdavg


def calc_thin2(
    u,
    hoxy = "OH2 O11 O12 O13 O14",
    lip= "resname POPC",
    coord_n = 26,
    coord_d = 1.0,
    coord_r = 8.0,
    coord_z = 0.75,
    padding = 0.5,
    coord_h = 0.25,
    r0 = 4.05,
    d0 = 10.0,
    v0 = 29.96,
    ep0 = 0.925,
    lmt_n = 12,
    lmt_k = 30,
):
    from MDAnalysis.analysis.distances import distance_array
    from MDAnalysis.lib.distances import calc_angles
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy
    from scipy.interpolate import make_smoothing_spline

    epsilon, epsilon_e, x0, y0, _ = calc_chain(u, lip=lip, coord_r = coord_r, coord_n=coord_n, hoxy=hoxy, coord_d=coord_d, coord_z=coord_z)

    lipid = u.select_atoms(f"{lip}")
    lipid_p = u.select_atoms(f"name P")
    if len(lipid_p.atoms) == 0:
        lipid_p = u.select_atoms(f"name PH")


    xyz = lipid_p.atoms.positions
    zavg = np.average(xyz[:, 2])
    radial = distance_array(np.array([x0, y0, 0]),
                            xyz*np.array([1, 1, 0]),
                            box = u.dimensions)[0]

    ups, dws = [], []           # group idxes
    ups_c, dws_c = [], []       # angles
    for i in np.argsort(radial):
        resid = u.select_atoms(f"resid {i+1}")
        com = resid.center_of_mass(unwrap=True)
        coords3 = np.array([xyz[i, 0], xyz[i, 1], com[2]])
        ang = calc_angles(xyz[i], com, coords3, box=u.dimensions)
        if xyz[i, 2] > zavg:
            ups.append(i)
            ups_c.append(ang)
        else:
            dws.append(i)
            dws_c.append(ang)
        if len(ups) >= 1 and len(dws) >= 1 and len(ups) + len(dws) >= lmt_n:
            break
    radial1 = radial[i]

    # get their positions, and get their pair stuff
    ups_pos = xyz[ups]
    dws_pos = xyz[dws]
    ip = np.array([(i, j) for i in range(len(ups)) for j in range(len(dws))]).T
    pdist = distance_array(ups_pos, dws_pos, box = u.dimensions * [1,1,10,1,1,1])
    dzs = []
    coms = []

    # sorted pdist indxes
    # print(np.arange(len(pdist)))
    # print(pdist.flatten())
    # plt.scatter(np.arange(len(pdist.flatten())), sorted(pdist.flatten()))
    y = np.sort(pdist.flatten())/np.max(pdist)
    x = np.arange(len(y))
    # fn = make_smoothing_spline(x[::3], y[::3], lam=None)
    # x0, y0 = [], []
    # y0 = np.interp(
    coef = np.polyfit(x, y, 20)
    fn = np.poly1d(coef)
    # plt.plot(x, y)
    # plt.plot(x, fn(x))

#     grad = np.gradient(sorted(pdist.flatten())[::10])
    grad = np.gradient(fn(x))
    # plt.plot(np.arange(len(pdist.flatten())), sorted(pdist.flatten())/np.max(pdist))
    ig = np.trapz(grad[:1000])
    print("ig", ig)
    # plt.plot(x, grad)
    return ig

    # print(len(pdist[0]))
    # plt.show()
    ind = np.unravel_index(np.argsort(pdist, axis=None), pdist.shape)
    for cnt, (i, j) in enumerate(zip(ind[0], ind[1])):
        dzs.append(ups_pos[i, 2] - dws_pos[j, 2])
        coms.append((ups_c[i] + dws_c[j])/2)
        if cnt > lmt_k:
            break
    # print(dzs)
    # exit()
    min0 = np.min(np.array(dzs)*np.array(coms))
    y = np.sort(np.array(dzs)*np.array(coms))
    # plt.plot(np.arange(len(dzs)), y)
    # plt.show()
    print(len(pdist[0]))

    lavg = np.average(dzs)
    d_lmt = lavg*np.average(coms)*(2/np.pi)
    return epsilon, epsilon_e, x0, y0, d_lmt


def unwrapper(posw0, posw1, posu0, box1) :
    "https://doi.org/10.1021/acs.jctc.3c00308"
    # posw0 :   wrapped pos at t0
    # posw1 :   wrapped pos at t1
    # posu0 : unwrapped pos at t0
    #  box1 : box at t1

    import numpy as np

    deltaw = (posw1 - posw0)
    posu1 = posu0 + deltaw - np.floor(deltaw/box1 + 0.5)*box1

    return posu1




