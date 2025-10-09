import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import glob
from MDAnalysis.analysis.distances import distance_array
#from dztools.funcs.plotter import COLS
#import scienceplots
#plt.style.use('science')
#plt.savefig('dump.pdf')
#from inftools.misc.data_helper import data_reader
#from inftools.misc.infinit_helper import read_toml

BOX = [12.4138, 12.4138, 12.4138, 90.0, 90.0, 90.0]


def get_frame(xyz, traj):
    """Finds the closest traj frame that matches xyz and returns it"""
    for frame in traj.trajectory[::-1]:
        delta = np.sum(np.abs(xyz[0] - frame.positions[0]))
        if delta < 10e-5:
            return traj.select_atoms("name X").positions
    assert False, "did not find correct frame!"


def calc_ruop(orig, homos):
    # place ru and o in an array and define count list
    at_ruo = orig.select_atoms("name Ru or name O")
    ru_dist = distance_array(at_ruo[0].position, at_ruo[1].position, box=BOX)[0][0]
    count = np.zeros(len(at_ruo))

    # calculate atom and homo distances, and populate count list
    dists = distance_array(at_ruo, homos, box=BOX)
    nearest_idxs = np.argmin(dists, axis=0)

    # iterate through the nearest neighbour
    at_exes = {i: [] for i in range(len(at_ruo))}
    for x_idx, at_idx in enumerate(nearest_idxs):
        count[at_idx] += 1
        at_exes[at_idx].append(dists[at_idx, x_idx])

    # find the ru1-x ru2-x pair distances
    status = 0 if 6 not in (count[0], count[1]) else int(count[1]-count[0])
    if status != 0:
        xinwhichru = 0 if count[0] == 6 else 1
        emptyru = int(not xinwhichru)

        dist1 = np.max(at_exes[xinwhichru])
        x_idx = np.where(dist1 == dists)[1][0]
        dist2 = dists[emptyru, x_idx]
    else:
        o_idx = np.argmax(count)
        ox_dist = np.max(at_exes[o_idx])
        x_idx = np.where(ox_dist == dists)[1][0]
        dist1 = dists[0, x_idx]
        dist2 = dists[1, x_idx]

    # calculate ru_op
    ru_op = (dist1 - dist2)/ru_dist
    return ru_op, x_idx, status


def calc_hop(orig):
    """"""
    at_r = orig.select_atoms("name Ru")
    at_o = orig.select_atoms("name O")
    at_h = orig.select_atoms("name H")
    oh_dists = distance_array(at_o, at_h, box=BOX)
    
    # get shortest o-h-o difference, min == 1
    sort =  np.sort(oh_dists, axis=0)
    oho_list = sort[1, :]/sort[0, :]

    # find hydrogen index
    nearest_idxs = np.argmin(oh_dists, axis=0)
    at_exes = {i: [] for i in range(len(at_o))}
    count = np.zeros(len(at_o))
    for h_idx, o_idx in enumerate(nearest_idxs):
        at_exes[o_idx].append(oh_dists[o_idx, h_idx])
        count[o_idx] += 1
    
    assert 1 in count

    # get h_idx and o_idx based on the complete atom list
    o_idx = np.argmin(count)
    h_idx = np.where(at_exes[o_idx][0] == oh_dists)[1][0]
    o_idx = at_o[o_idx].index
    h_idx = at_h[h_idx].index

    # calculate ru-h-ru
    rh_dists = distance_array(at_r, orig.select_atoms(f"index {h_idx}"), box=BOX)
    hop = rh_dists[1][0] - rh_dists[0][0]

    # calculate solvation shell
    ro_dists = distance_array(at_r, at_o, box=BOX)
    sort = np.sort(ro_dists, axis=1)
    solv1 = np.average(sort[0][:6])
    solv2 = np.average(sort[1][:6])

    return hop, h_idx, o_idx, np.sort(oho_list)[0], solv1, solv2


