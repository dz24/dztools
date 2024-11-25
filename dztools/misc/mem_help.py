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
    # xnew2 = np.zeros(len(x))

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

def f_radial(xi, yi, Xcyl, Ycyl, Rcyl, box, pos, h=1/4):
    import numpy as np
    # import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    ri = np.sqrt((xi-Xcyl)**2 + (yi-Ycyl)**2)
    pos[:, 2] = 0
    #print('clean', type(pos))
    ri1 = distances.distance_array(np.array([Xcyl, Ycyl, 0]), pos, box=box)
    for i, j in zip(ri, ri1[0]):
        d = abs(i-j)
        print(d < 0.001, abs(i-j))
    print('clown', ri)
    print('clown', ri1)
    return theta(ri1[0]/Rcyl, h)
