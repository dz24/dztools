from typing import Annotated
import typer


def vel_rev(
    i: Annotated[str, typer.Option("-i", help="input")],
    o: Annotated[str, typer.Option("-o", help="output")] = False
):
    import MDAnalysis as mda
    import os

    # load input to mda universe
    u = mda.Universe(i)

    # check if universe has velocities
    try:
        u.atoms.velocities
    except:
        print(f"{i} likely does not contain any velocities.")
        return

    # reverse velocities
    u.atoms.velocities *= -1

    # prep write
    atoms = u.select_atoms("all")
    split = i.split(".")
    outf = ".".join(split[:-1]) + "_v." + split[-1]
    outf = o if o != "False" else outf
    assert not os.path.isfile(outf)

    # write
    atoms.write(outf)
