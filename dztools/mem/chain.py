from typing import Annotated

import typer

""" Implementation of https://pubs.acs.org/doi/10.1021/acs.jctc.7b00106"""

def mem_chain(
    gro: Annotated[str, typer.Option("-gro", help="gro file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")],
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
    coord_n: Annotated[int, typer.Option("-lip", help="n")] = 26,
    coord_d: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.1,
    coord_r: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.9,
    coord_z: Annotated[float, typer.Option("-lip", help="xtc file")] = 0.75,
):
    """Currently for MEL system only."""
    import matplotlib.pyplot as plt
    import numpy as np

    import MDAnalysis as mda
    from MDAnalysis.analysis import helix_analysis as hel

    # load gro and xtc into MDA
    u = mda.Universe(gro, xtc)

    # group:

    # for idx, ts in enumerate(u.trajectory):

    # # Select individual proteins and membrane
    # protein = u.select_atoms("protein")
    # nores = 26
    # noprot = int(len(protein.residues)/nores)
    # print(len(protein.residues))
    # mels = []
    # for i in range(noprot):
    #     print(nores*i, nores*(i+1), f"name CA and resnum {nores*i}-{nores*(i+1)}")
    #     h = hel.HELANAL(u, select=f"name CA and resid {nores*i}:{nores*(i+1)}").run()
    #     # h = hel.HELANAL(u, select=f"resid {nores*i}:{nores*(i+1)} and name CA", ref_axis=[0, 0, 1])
    #     print('boomer 1', dir(h))
    #     print('boomer 2', h.results)
    #     return
    #     plt.plot(h.results.local_twists.mean(axis=1))
    #     plt.xlabel('Frame')
    #     plt.ylabel('Average twist (degrees)')

        # hel.HELANAL(u, select=protein.residues[nores*i:nores*(i+1)])
        # mels.append(protein.residues[nores*i:nores*(i+1)])

    # helicity
    # lipid   = u.select_atoms(f"resname {lip}")

    # define idxs and coms and iterate over all xtc frames
    # idxs, box_z, lcoms_z, pcoms_z1 = [], [], [], [[] for _ in range(noprot)]
    # for idx, ts in enumerate(u.trajectory):
    #     idxs.append(idx)
    #     box_z.append(ts.dimensions[2])
    #     lcoms_z.append(lipid.atoms.center_of_mass()[-1])
    #     for idx0, mel in enumerate(mels):
    #         pcoms_z1[idx0].append(mel.center_of_mass(unwrap=True)[-1])

    # # change to np array
    # lcoms_z, box_z, pcoms_z2 = np.array(lcoms_z), np.array(box_z), []
    # for idx0, pcom_z in enumerate(pcoms_z1):
    #     pcoms_z2.append(np.array(pcom_z))

    # # calculate CVs
    # CV0 = lcoms_z*0
    # for pcom in pcoms_z2:
    #     CV0 += abs(pcom - lcoms_z)
    # CV0 = -CV0/noprot

    # # plot
    # plt.plot(idxs, lcoms_z , color='k', lw=2.0)
    # plt.plot(idxs, lcoms_z + box_z/2, color='k', ls='--', lw=2.0)
    # plt.plot(idxs, lcoms_z - box_z/2, color='k', ls='--', lw=2.0)
    # plt.plot(idxs, lcoms_z + 20, color='k', ls='--', lw=2.0, alpha=0.2)
    # plt.plot(idxs, lcoms_z - 20, color='k', ls='--', lw=2.0, alpha=0.2)
    # for i in range(noprot):
    #     plt.plot(idxs, pcoms_z2[i], c=f"C{i+1}", label=f"p{i+1}",
    #              alpha=0.2)
    #     plt.plot(idxs, -abs(pcoms_z2[i] - lcoms_z)+ lcoms_z, c=f"C{i+1}",
    #              label=f"p{i+1}", alpha=0.2)
    # plt.plot(idxs, CV0 + lcoms_z, color='r', lw=2.0)
    # plt.show()
