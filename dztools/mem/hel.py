from typing import Annotated

import typer


def helicity(
    gro: Annotated[str, typer.Option("-gro", help="gro file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = None,
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
):
    """STRIDE wrapper."""
    import subprocess

    import numpy as np
    import matplotlib.pyplot as plt

    import MDAnalysis as mda
    from MDAnalysis.analysis.dssp import DSSP

    u = mda.Universe(gro, xtc)

    protein = u.select_atoms("protein")
    nores = 26
    noprot = int(len(protein.residues)/nores)
    # print(len(protein.residues))

    s = DSSP(u).run().results.dssp
    print(len(s))
    print(s[0])
    print(s[-1])

    exit('a')



    mels = []
    for i in range(noprot):
        mels.append(protein.residues[nores*i:nores*(i+1)])

    for idx, ts in enumerate(u.trajectory):

        # DSSP(u).run().results.dssp[0]
        # s = DSSP(u).run().results.dssp[0]
        s = DSSP(mels[0]).run().results.dssp[0]

    # print('tiger a', sum(resi), len(resi))
        print('tiger b', sum(np.array(s)=='H'), len(s))
    return sum(resi)/(len(resi)-2)
