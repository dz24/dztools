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

    return epsilon
