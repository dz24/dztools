import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})

BOX = [12.4138, 12.4138, 12.4138, 90.0, 90.0, 90.0]
AT, O, RU = 193, 64, 2
NAMES = ['Ru']*RU + ['O']*O + ['H']*(AT-O-RU)

def finder_a(pos):
    """Calculate something ."""
    atom_array = [0] * (RU+O)
    dic = {i: [] for i in range(RU+O)}
    for i in range(len(pos[AT:][:, 0])):
        dist_arr = distances.distance_array(pos[AT+i],
                                            pos[:RU+O],
                                            box=BOX)[0]
        loc = np.argmin(dist_arr)
        atom_array[loc] += 1
        dic[loc].append({'x_loc': i+AT,
                         'dist': dist_arr[loc],
                         'a_loc': loc})
    return atom_array, dic

def calc_orderp(pos):
    # orderps = []
    # uni = mda.Universe(inp)
    # for pos in uni.trajectory:
    at_ar, dic = finder_a(pos)

    # Excess electron exists in one of the Ru atoms
    if 6 in at_ar[0:2]:
        a_idx = 0 if at_ar[0] == 6 else 1
        x_idx = np.argmax([i['dist'] for i in dic[a_idx]])
    # Excess electron exists in one of the oxygen atoms.
    else:
        # Get list for oxygen with max count(x):
        o_l = [i for i, j in enumerate(at_ar) if j == np.max(at_ar)]
        a_idx0, x_idx0, dist = 0, 0, 0
        for a_idx0 in o_l:
            x_idx0 = np.argmax([i['dist'] for i in dic[a_idx0]])
            if dist < dic[a_idx0][x_idx0]['dist']:
                dist = dic[a_idx0][x_idx0]['dist']
                a_idx, x_idx = a_idx0, x_idx0

    if not at_ar[0] == 523:
        loc = dic[a_idx][x_idx]['x_loc']
        dist1 = distances.distance_array(pos[0], pos[loc], box=BOX)
        dist2 = distances.distance_array(pos[1], pos[loc], box=BOX)
        dist3 = distances.distance_array(pos[0], pos[1], box=BOX)
        orderp = (dist1[0][0]-dist2[0][0])/dist3[0][0]
    else:
        orderp = -1
        loc = -1
    # orderps.append(orderp)
    # return orderps
    return orderp, loc
