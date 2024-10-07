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

    u = mda.Universe(gro)
    s = DSSP(u).run().results.dssp[0]

    stride = "/home/daniel/Documents/programs/stride/stride"
    out = subprocess.run([stride, gro], capture_output=True, text=True)
    resi = []
    head = True

    for line in out.stdout.rstrip().split("\n"):
        if head:
            if "---Residue---" in line:
                head = False
            continue
        cut = line.rstrip().split()
        resi.append(1 if cut[6]=="AlphaHelix" else 0)

    print('tiger a', sum(resi), len(resi))
    print('tiger b', sum(np.array(s)=='H'), len(s))
    return sum(resi)/(len(resi)-2)
