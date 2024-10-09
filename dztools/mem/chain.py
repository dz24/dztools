from typing import Annotated

import typer

import numpy as np

""" Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""


def mem_chain(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
    coord_n: Annotated[int, typer.Option("-lip", help="n")] = 26, # ns
    coord_d: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.1*10,
    coord_r: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.9,
    coord_z: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.75,
):
    """Currently for MEL system only."""
    import matplotlib.pyplot as plt
    import numpy as np

    import MDAnalysis as mda
    from MDAnalysis.analysis import helix_analysis as hel

    # load top and xtc into MDA
    u = mda.Universe(top, xtc)
    hoxys = u.select_atoms("name OH2 O11 O12 O13 O14")
    lipid = u.select_atoms(f"resname POPC")

    for idx, ts in enumerate(u.trajectory):
        z_mem = lipid.atoms.center_of_mass()[-1]
        z_s = [z_mem + (s + 1/2 - coord_n)*coord_d for s in range(coord_n)]
        box = ts.dimensions

        nsp = [0 for i in range(coord_n)]
        f_norm = [0 for i in range(coord_n)]
        fx_acc_s = [0 for i in range(coord_n)]
        fx_acc_c = [0 for i in range(coord_n)]
        fy_acc_s = [0 for i in range(coord_n)]
        fy_acc_c = [0 for i in range(coord_n)]
        ws_cyl = [0 for i in range(coord_n)]
        theta_ix = [0 for i in range(coord_n)]
        x_sincyl, x_coscyl = 0, 0
        y_sincyl, y_coscyl = 0, 0

        for s in range(coord_n):
            for hoxy in hoxys.atoms.positions:
                fx_acc_s[s] += f_axial(hoxy[-1], z_s[s],
                                       coord_d)*np.sin(2*np.pi*hoxy[0]/box[0])
                fx_acc_c[s] += f_axial(hoxy[-1], z_s[s],
                                       coord_d)*np.cos(2*np.pi*hoxy[0]/box[0])
                fy_acc_s[s] += f_axial(hoxy[-1], z_s[s],
                                       coord_d)*np.sin(2*np.pi*hoxy[1]/box[1])
                fy_acc_c[s] += f_axial(hoxy[-1], z_s[s],
                                       coord_d)*np.cos(2*np.pi*hoxy[1]/box[1])
                f_norm[s] += f_axial(hoxy[-1], z_s[s], coord_d)
            ws_cyl[s] = np.tanh(f_norm[s])
            print('fl', fy_acc_c[s], f_norm[s])
            x_sincyl += ws_cyl[s]*fx_acc_s[s]/f_norm[s]
            x_coscyl += ws_cyl[s]*fx_acc_c[s]/f_norm[s]
            y_sincyl += ws_cyl[s]*fy_acc_s[s]/f_norm[s]
            y_coscyl += ws_cyl[s]*fy_acc_c[s]/f_norm[s]
            # exit('')
        x_sincyl /= np.sum(ws_cyl)
        x_coscyl /= np.sum(ws_cyl)
        y_sincyl /= np.sum(ws_cyl)
        y_coscyl /= np.sum(ws_cyl)
        x_cyl = (np.arctan2(-x_sincyl,-x_coscyl) + np.pi)*box[0]/(2*np.pi)
        y_cyl = (np.arctan2(-x_sincyl,-y_coscyl) + np.pi)*box[0]/(2*np.pi)


        print('umba x', box[0], x_cyl)
        print('umba y', box[1], y_cyl)
        print('zcom', z_mem)

        exit()




        print(hoxys.atoms.positions)
        # w_norm =
        # f_norm =
        print('chez')
        exit()


    # X_cyl = 
    # z_com_mem =
    # z_mem = 
    # epsilon =
    # X_cyl = 
    # Y_cyl = 


    # group:


def s_center():
    print()


def psi_switch(x, zeta):
    if x <= 1:
        return x*zeta
    else:
        b = zeta/(1-zeta)
        c = (1-zeta)*np.exp(b)
        return 1 - c*np.exp(-b*x)

# def nsp():
#     print()

def theta(x, h):
    if -1 + h <= x <= 1 - h:
        return 1
    elif 1 - h < x < 1 + h:
        return 1/2 - (3/(4*h))*(x-1) + (1/(4*h**3))*(x-1)**3
    elif -1 - h < x < -1 + h:
        return 1/2 + (3/(4*h))*(x+1) - (1/(4*h**3))*(x+1)**3
    else:
        # print('pig', x)
        return 0

def f_axial(zi, zs, ds, h=1/4):
    return theta((zi-zs)/(ds/2),h)

def f_radial(xi, yi, Xcyl, Ycyl, Rcyl):
    ri = np.sqrt((xi-Xcyl)**2 + (yi-Ycyl)**2)
    return theta(ri/Rcyl, h)

        
