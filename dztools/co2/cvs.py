from typing import Annotated as And

# import typer
from typer import Option as Opt


def co2_op(
    xyz: And[str, Opt("-xyz", help="xyz file")],
    box: And[float, Opt("-box", help="box")],
    plot: And[bool, Opt("-plot", help="plot")] = False,
    out: And[str, Opt("-out", help="string")] = "co2_op.txt",
):
    import numpy as np
    import matplotlib.pyplot as plt
    import MDAnalysis as mda
    from MDAnalysis.analysis.distances import distance_array
    from MDAnalysis.lib.distances import calc_angles, calc_dihedrals


    box = np.array([box, box, box, 90, 90, 90])
    u = mda.Universe(xyz)

    carbon = u.select_atoms("element C")
    oxygens = u.select_atoms("element O")
    hydrogens = u.select_atoms("element H")

    results = []
    for idx0, ts in enumerate(u.trajectory):

        # stable
        d_HO = distance_array(
        hydrogens.positions,
        oxygens.positions,
        box=box
        )

        # sort O distances for each H
        d_sorted = np.sort(d_HO, axis=1)
        stable_f = np.sort(d_sorted[:,1]/d_sorted[:,0])
        # stable = stable_f[2]

        pairs = [[carbon[0].index] for _ in range(3)]
        oh_dists = []
        stables = []
        d_CO = distance_array(carbon, oxygens, box=box)[0]
        oidxes = np.argsort(d_CO)
        for idx, oidx in enumerate(oidxes[:3]):
            pairs[idx].append(oxygens[oidx].index)
            d_OH = distance_array(oxygens[oidx].position, hydrogens, box=box)[0]
            asort = np.argsort(d_OH)
            stables.append(stable_f[asort[0]])
            oh_dists.append(d_OH[asort[0]])
            pairs[idx].append(hydrogens[asort[0]].index)

        angs = []
        for tup in [(0, 1), (0, 2), (1, 2)]:
            angle = calc_angles(
            u.atoms[pairs[tup[0]][1]].position,
            u.atoms[pairs[0][0]].position,
            u.atoms[pairs[tup[1]][1]].position,
            box=u.dimensions
            )
            angs.append(np.rad2deg(angle))

        cdouble2 = np.array([np.argmax(oh_dists)])
        csingle2 = np.array([i for i in range(3) if i not in cdouble2])
        oxy0 = pairs[cdouble2[0]]
        ang1 = pairs[csingle2[0]]
        ang2 = pairs[csingle2[1]]

        c_xyz  = u.atoms[carbon[0].index].position
        h1_xyz = u.atoms[ang1[2]].position
        h2_xyz = u.atoms[ang2[2]].position
        v = h2_xyz - h1_xyz
        vv = np.dot(v, v)
        t = np.dot(c_xyz - h1_xyz, v) / vv
        p = h1_xyz + t * v
        chh_line = np.linalg.norm(c_xyz - p)

        phi1 = calc_dihedrals(
            u.atoms[oxy0[1]].position,
            u.atoms[oxy0[0]].position,
            u.atoms[ang1[1]].position,
            u.atoms[ang1[2]].position,
            box=box
        )
        phi2 = calc_dihedrals(
            u.atoms[oxy0[1]].position,
            u.atoms[oxy0[0]].position,
            u.atoms[ang2[1]].position,
            u.atoms[ang2[2]].position,
            box=box
        )
        dihed1 = np.rad2deg(phi1)
        dihed2 = np.rad2deg(phi2)

        results.append([
            idx0,              # time
            np.max(angs),     # max o-c-o angle
            chh_line,
            dihed1,           # dihed angle 1
            dihed2,           # dihed angle 2
            oxy0[1],          # o idx 1 
            ang1[1],          # o idx 2 
            ang2[1],          # o idx 3 
            oxy0[2],          # h idx 1
            ang1[2],          # h idx 2
            ang2[2],          # h idx 3
            d_CO[oidxes[0]],  # c-o dist 1 
            d_CO[oidxes[1]],  # c-o dist 2 
            d_CO[oidxes[2]],  # c-o dist 3 
            oh_dists[0],      # o-h dist 1 
            oh_dists[1],      # o-h dist 2 
            oh_dists[2],      # o-h dist 3 
            stable_f[0],      # smallest 1 o-h bond
            stable_f[1],      # smallest 2 o-h bond
            stable_f[2],      # smallest 3 o-h bond
        ])
        print(results[-1][0], results[-1][1])
    results = np.array(results)
    np.savetxt(out, results)
