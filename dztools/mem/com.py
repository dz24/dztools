from typing import Annotated

import typer

from dztools.misc.xyz_help import calc_center2


def com_6met(
    gro: Annotated[str, typer.Option("-gro", help="gro file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
):
    """Currently for 6MEL system only."""
    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    import numpy as np

    # load gro and xtc into MDA
    u = mda.Universe(gro)
    u.load_new(xtc)

    # Select individual proteins and membrane
    protein = u.select_atoms("protein")
    melres = int(len(protein.residues)/6)
    mels = []
    for i in range(6):
        mels.append(protein.residues[melres*i:melres*(i+1)])
    lipid   = u.select_atoms(f"resname {lip}")

    # define idxs and coms and iterate over all xtc frames
    idxs, box_z, lcoms_z, pcoms_z1 = [], [], [], [[] for _ in range(6)]
    for idx, ts in enumerate(u.trajectory):
        idxs.append(idx)
        box_z.append(ts.dimensions[2])
        lcoms_z.append(lipid.atoms.center_of_mass()[-1])
        for idx0, mel in enumerate(mels):
            pcoms_z1[idx0].append(mel.center_of_mass()[-1])

        if idx == 50:
            break

    # change to np array
    lcoms_z, pcoms_z2 = np.array(lcoms_z), []
    for idx0, pcom_z in enumerate(pcoms_z1):
        pcoms_z2.append(np.array(pcom_z))
        # plt.plot(idxs, pcoms_z2[-1], c=f"C{idx0+1}", label=f"p{idx0+1}")
        pcoms_z2[-1] += (lcoms_z - pcoms_z2[-1] > 0)*box_z
        plt.plot(idxs, pcoms_z2[-1]- lcoms_z[0], c=f"C{idx0+1}", label=f"p{idx0+1}")
        # plt.plot(idxs, np.abs(np.array(pcom_z)), c=f"C{idx0+1}", label=f"p{idx0+1}")
    plt.plot(idxs, lcoms_z - lcoms_z[0], color='C0', lw=2.0)
    plt.plot(idxs, lcoms_z+box_z - lcoms_z[0], color='C0', lw=2.0)
    plt.show()
    exit()

    # change to np array ##, shift and add pbc
    newzzero = lcoms_z[0]
    box_z = np.array(box_z)
    # lcoms_z = np.array(lcoms_z) - newzzero
    lcoms_z = np.array(lcoms_z)
    plt.plot(idxs, lcoms_z, color='C0', lw=2.0)
    # pcoms_z2 = []
    # for idx0, pcom_z in enumerate(pcoms_z1):
    #     plt.plot(idxs, np.abs(np.array(pcom_z)), c=f"C{idx0+1}", label=f"p{idx0+1}")
    #     pcoms_z2.append(np.array(pcom_z) - newzzero)
    #     # add pbc
    #     print('dawg', sum(pcoms_z2[-1]<0))
    #     pcoms_z2[-1] += (pcoms_z2[-1]<0)*box_z
    #     exit('dog')

    # plt.plot(idxs, np.array(lcoms_z)*0, color='C0', lw=2.0)
    # plt.plot(idxs, np.array(box_z)*0+20, ls='--', color='k', lw=2.0)
    # plt.plot(idxs, np.array(box_z)*0-20, ls='--', color='k', lw=2.0)
    # plt.plot(idxs, np.array(box_z), color='C0', lw=2.0)
    # plt.plot(idxs, -np.array(box_z), color='C0', lw=2.0)
    # for i in range(6):
    #     plt.plot(idxs, np.abs(np.array(pcoms_z[0])-lcoms_z), c=f"C{i+1}", label=f"p{i+1}")
    plt.legend(frameon=False)
    plt.show()
    print('whadab', len(u.trajectory))
