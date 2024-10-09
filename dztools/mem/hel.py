from typing import Annotated

import typer


def helicity(
    top: Annotated[str, typer.Option("-top", help="gro/pdb/tpr file")],
    xtc: Annotated[str, typer.Option("-xtc", help="xtc file")] = None,
    lip: Annotated[str, typer.Option("-lip", help="xtc file")] = "POPC",
):
    """STRIDE wrapper."""
    import subprocess

    import numpy as np
    import matplotlib.pyplot as plt

    import MDAnalysis as mda
    from MDAnalysis.analysis.dssp import DSSP

    u = mda.Universe(top, xtc)

    protein = u.select_atoms("protein")
    nores = 25
    noprot = int(len(protein.residues)/nores)
    mels = []
    print('len', len(protein.residues))
    for i in range(noprot):
        # mels.append(u.select_atoms(f"resid {nores*i}:{nores*(i+1)}"))
        mel = u.select_atoms(f"resid {nores*i}:{nores*(i+1)}")
        helicity = []
        x = []
        for x0, f in enumerate(DSSP(mel).run().results.dssp):
            heli = sum(np.array(f)=='H')/nores
            if heli > 0:
                helicity.append(sum(np.array(f)=='H')/nores)
                x.append(x0)
        plt.scatter(x, helicity)

    plt.show()
