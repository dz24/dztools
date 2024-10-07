from typing import Annotated

import typer


def helicity(
    gro: Annotated[str, typer.Option("-gro", help="gro file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = None,
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
):
    """Currently for MEL system only."""
    import matplotlib.pyplot as plt
    import numpy as np

    import MDAnalysis as mda
    from MDAnalysis.analysis import helix_analysis as hel
    import subprocess 

    stride = "/home/daniel/Documents/programs/stride/stride"
    out = subprocess.run([stride, gro], capture_output=True, text=True)
    resi = []
    head = True
    # print(dir(out))
    # print(out.stdout)
    # print(type(out.stdout))
    # exit('a')
    for line in out.stdout.rstrip().split("\n"):
        if head:
            if "---Residue---" in line:
                head = False
            continue
        cut = line.rstrip().split()
        # print('prime', cut)
        resi.append(1 if cut[6]=="AlphaHelix" else 0)
        # print('cut', cut[0], cut[6])

    return sum(resi)/(len(resi)-2)
    # print(sum(resi)/(len(resi)-2))
    # exit()
    # load gro and xtc into MDA
    # u = mda.Universe(gro, xtc)
    u = mda.Universe(gro)

    # Select individual proteins and membrane
    protein = u.select_atoms("protein")
    nores = 26
    noprot = int(len(protein.residues)/nores)
    mels = []
    # for i in range(noprot):
    #     print(nores*i, nores*(i+1), f"name CA and resnum {nores*i}-{nores*(i+1)}")
    h = hel.HELANAL(u, select=f"name CA and resid {nores*0}:{nores*(0+1)}").run()
    res = h.results.summary
    keys = h.results.summary.keys()
    # print('w', res["local_bends"]["mean"])
    print('w', res["local_twists"]["mean"])
    # for key in keys:
        # dic = h.results.summary[key].items()
        # print(key, dic)
        # print(key)
        # for key0, val in h.results.summary[key].items():
        #     print(key, key0, val.shape)
        #     print(val)
        #     break
        # break
            # print(key, f"{key0}: {val:.3f}")

    # h = hel.HELANAL(u, select=f"resid {nores*i}:{nores*(i+1)} and name CA", ref_axis=[0, 0, 1])
    # print('boomer 1', h.results.all_bends.shape)
    # print('boomer 2', h.results.summary.keys())
    # print('boomer 3', h.results.summary)
        # print('boomer 2', h.results)
        # plt.plot(h.results.local_twists.mean(axis=1))
        # plt.xlabel('Frame')
        # plt.ylabel('Average twist (degrees)')
        # plt.show()
        # return

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
