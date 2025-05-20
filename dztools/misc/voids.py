import numpy as np
import matplotlib.pyplot as plt

import MDAnalysis as mda
from MDAnalysis.lib.pkdtree import PeriodicKDTree

def leaflet_void(pos, w, dim, X, Y, ax, rad):
    tree = PeriodicKDTree(box=dim)
    tree.set_coords(pos, cutoff=40) # pbc cutoff or something
    dots = tree.search_tree(w, radius=rad)
    boxx = dim[0]

    if len(dots) > 0:
        all_indices = np.arange(len(X.flatten()))
        idxes = np.setdiff1d(all_indices, dots[:, 0])
        ax.scatter(X.flatten()[idxes], Y.flatten()[idxes], color='k', rasterized=True)

        ax.scatter(X.flatten()[idxes] + boxx, Y.flatten()[idxes], color='k', alpha=0.2, rasterized=True)
        ax.scatter(X.flatten()[idxes] - boxx, Y.flatten()[idxes], color='k', alpha=0.2, rasterized=True)
        ax.scatter(X.flatten()[idxes], Y.flatten()[idxes] + boxx, color='k', alpha=0.2, rasterized=True)
        ax.scatter(X.flatten()[idxes], Y.flatten()[idxes] - boxx, color='k', alpha=0.2, rasterized=True)
        ax.scatter(X.flatten()[idxes] + boxx, Y.flatten()[idxes] + boxx, color='k', alpha=0.2, rasterized=True)
        ax.scatter(X.flatten()[idxes] - boxx, Y.flatten()[idxes] + boxx, color='k', alpha=0.2, rasterized=True)
        ax.scatter(X.flatten()[idxes] + boxx, Y.flatten()[idxes] - boxx, color='k', alpha=0.2, rasterized=True)
        ax.scatter(X.flatten()[idxes] - boxx, Y.flatten()[idxes] - boxx, color='k', alpha=0.2, rasterized=True)

    ax.scatter(pos[:, 0], pos[:, 1], color="C0")
    ax.scatter(pos[:, 0] + boxx, pos[:, 1], color="C1", alpha=0.2)
    ax.scatter(pos[:, 0] - boxx, pos[:, 1], color="C1", alpha=0.2)
    ax.scatter(pos[:, 0], pos[:, 1] + boxx, color="C1", alpha=0.2)
    ax.scatter(pos[:, 0], pos[:, 1] - boxx, color="C1", alpha=0.2)
    ax.scatter(pos[:, 0] + boxx, pos[:, 1] + boxx, color="C1", alpha=0.2)
    ax.scatter(pos[:, 0] - boxx, pos[:, 1] + boxx, color="C1", alpha=0.2)
    ax.scatter(pos[:, 0] + boxx, pos[:, 1] - boxx, color="C1", alpha=0.2)
    ax.scatter(pos[:, 0] - boxx, pos[:, 1] - boxx, color="C1", alpha=0.2)

    ax.set_xlim([0 - boxx*0.2, boxx + boxx*0.2])
    ax.set_ylim([0 - boxx*0.2, boxx + boxx*0.2])
    ax.axhline(0, color="r")
    ax.axhline(boxx, color="r")
    ax.axvline(0, color="r")
    ax.axvline(boxx, color="r")
    ax.set_xlabel("X [nm]")
    ax.set_ylabel("Y [nm]")

    return len(X.flatten()[idxes])

rad = 4
top = "conf.gro"
xtc = "pbc.xtc"
out = f"void_r{rad}.pdf"

u = mda.Universe(top, xtc)

lipid_p = u.select_atoms(f"name P")
samples = 100

for idx, ts in enumerate(u.trajectory):
    # get standard info from frame
    pos = lipid_p.atoms.positions
    boxx = ts.dimensions[0]
    zavg = np.average(pos[:, 2])
    upper = lipid_p.select_atoms(f"prop z > {zavg}")
    lower = lipid_p.select_atoms(f"prop z <= {zavg}")

    # project P into x, y
    pos_u = np.copy(upper.atoms.positions * [1, 1, 0])
    pos_l = np.copy(lower.atoms.positions * [1, 1, 0])

    # create 2D xy dotted grid
    xlim = np.linspace(0, boxx, samples)
    X, Y = np.meshgrid(xlim, xlim)  # Create 2D grids
    w = np.array([X.flatten(), Y.flatten(), np.zeros(samples**2)]).T

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))

    # use periodic PeriodicKDTree to search neigbours between P and dots
    cnt_u = leaflet_void(pos_u, w, ts.dimensions, X, Y, axs[0], rad=rad)
    cnt_l = leaflet_void(pos_l, w, ts.dimensions, X, Y, axs[1], rad=rad)

    axs[0].set_title(f"Upper leaflet, void dots: {cnt_u}, radius {rad}")
    axs[1].set_title(f"Lower leaflet, void dots: {cnt_l}, radius {rad}")

    plt.savefig(out, bbox_inches='tight', dpi=300)
    exit()
