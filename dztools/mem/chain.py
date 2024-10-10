from typing import Annotated

import typer

import numpy as np
import time

""" Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""


def mem_chain(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
    coord_n: Annotated[int, typer.Option("-lip", help="n")] = 26, # ns
    coord_d: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.1*10,
    coord_r: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.9*10,
    coord_z: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.75,
):
    """Currently for MEL system only."""
    import matplotlib.pyplot as plt
    import numpy as np

    import MDAnalysis as mda
    from MDAnalysis.analysis import helix_analysis as hel

    # load top and xtc into MDA
    time0 = time.perf_counter()
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

        time0 = time.perf_counter()
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
        # exit('brom')


def psi_switch(x, zeta):
    if x <= 1:
        return x*zeta
    else:
        b = zeta/(1-zeta)
        c = (1-zeta)*np.exp(b)
        return 1 - c*np.exp(-b*x)

def theta(x, h):
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
    # print('pringle', zi, zs, (zi-zs))
    return theta((zi-zs)/(ds/2),h)

def f_radial(xi, yi, Xcyl, Ycyl, Rcyl, h=1/4):
    ri = np.sqrt((xi-Xcyl)**2 + (yi-Ycyl)**2)
    return theta(ri/Rcyl, h)
