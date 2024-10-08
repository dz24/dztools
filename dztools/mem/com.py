from typing import Annotated

import typer


def com_met(
    pdb: Annotated[str, typer.Option("-pdb", help="pdb file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="lip file")] = "POPC",
):
    """Currently for MEL system only."""
    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    import numpy as np

    # load gro and xtc into MDA
    u = mda.Universe(pdb, xtc)
    # protein = u.select_atoms('protein')
    # for ts in u.trajectory:
    #     protein.unwrap(compound='fragments')

    # exit('a')
    # Select individual proteins and membrane
    # print(protein.atoms.center_of_mass(unwrap=True))
    protein = u.select_atoms("protein")
    nores = 26
    noprot = int(len(protein.residues)/nores)
    mels = []
    for i in range(noprot):
        # mels.append(protein.residues[nores*i:nores*(i+1)])
        # mels.append(protein.residues[nores*i:nores*(i+1)])
        mels.append(u.select_atoms(f"resid {nores*i}:{nores*(i+1)}"))
        print("tiger", type(protein), type(mels[-1]))
        protein.atoms.center_of_mass(unwrap=True)
        mels[-1].atoms.center_of_mass(unwrap=True)
        # protein.unwrap(compound="fragments")
        # mels[-1].unrwap(compound="fragments")
    lipid   = u.select_atoms(f"resname {lip}")

    # define idxs and coms and iterate over all xtc frames
    idxs, box_z, lcoms_z, pcoms_z1 = [], [], [], [[] for _ in range(noprot)]
    for idx, ts in enumerate(u.trajectory):
        idxs.append(idx)
        box_z.append(ts.dimensions[2])
        lcoms_z.append(lipid.atoms.center_of_mass()[-1])
        for idx0, mel in enumerate(mels):
            pcoms_z1[idx0].append(mel.center_of_mass(unwrap=True)[-1])

    # change to np array
    lcoms_z, box_z, pcoms_z2 = np.array(lcoms_z), np.array(box_z), []
    for idx0, pcom_z in enumerate(pcoms_z1):
        pcoms_z2.append(np.array(pcom_z))

    # calculate CVs
    CV0 = lcoms_z*0
    for pcom in pcoms_z2:
        CV0 += abs(pcom - lcoms_z)
    CV0 = -CV0/noprot

    # plot
    plt.plot(idxs, lcoms_z , color='k', lw=2.0)
    plt.plot(idxs, lcoms_z + box_z/2, color='k', ls='--', lw=2.0)
    plt.plot(idxs, lcoms_z - box_z/2, color='k', ls='--', lw=2.0)
    plt.plot(idxs, lcoms_z + 20, color='k', ls='--', lw=2.0, alpha=0.2)
    plt.plot(idxs, lcoms_z - 20, color='k', ls='--', lw=2.0, alpha=0.2)
    for i in range(noprot):
        plt.plot(idxs, pcoms_z2[i], c=f"C{i+1}", label=f"p{i+1}",
                 alpha=0.2)
        plt.plot(idxs, -abs(pcoms_z2[i] - lcoms_z)+ lcoms_z, c=f"C{i+1}",
                 label=f"p{i+1}", alpha=0.2)
    plt.plot(idxs, CV0 + lcoms_z, color='r', lw=2.0)
    plt.show()
