def calc_met_com(mdau, lip, num_resi):
    """Calculates individual and average MELITTIN COMS wrt. specified lipid
    membrane."""
    import numpy as np

    # Select individual proteins and membrane
    lipid = mdau.select_atoms(f"resname {lip}")
    protein = mdau.select_atoms("protein")
    mels = []
    for i in range(int(len(protein.residues) / num_resi)):
        mels.append(mdau.select_atoms(f"resid {num_resi*i}:{num_resi*(i+1)}"))

    # define coms and iterate over all xtc frames
    lcoms_z, pcoms_z_avg = [], []
    pcoms_z = [[] for _ in range(len(mels))]
    for idx, ts in enumerate(mdau.trajectory):
        pcoms_z_avg.append(0)
        lcoms_z.append(lipid.atoms.center_of_mass()[-1])
        for idx0, mel in enumerate(mels):
            pcoms_z[idx0].append(mel.center_of_mass()[-1])
            pcoms_z_avg[-1] += abs(pcoms_z[idx0][-1] - lcoms_z[-1])
        pcoms_z_avg[-1] /= len(mels)

    return np.array(pcoms_z), np.abs(np.array(pcoms_z_avg)), lcoms_z


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

def f_radial(xi, yi, Xcyl, Ycyl, Rcyl, h=1/4):
    import numpy as np
    ri = np.sqrt((xi-Xcyl)**2 + (yi-Ycyl)**2)
    return theta(ri/Rcyl, h)
