
import numpy as np
from MDAnalysis.analysis.distances import distance_array as dar 

# Assign global variables:
BOX = [12.4138, 12.4138, 12.4138, 90.0, 90.0, 90.0]
AT, O, RU = 193, 64, 2
NAMES = ['RU']*RU + ['O']*O + ['H']*(AT-O-RU)


def locate(xyzs):
    # ru1, ru2 = xyz[0], xyz[1]
    # center = (ru1 + ru2)/2
    three = np.array([xyzs[0], xyzs[1], (xyzs[0] + xyzs[1])/2])
    dlist = []
    for idx, xyz in enumerate(xyzs[0:O+RU]):
        dist_arr = dar(xyz, three, box=BOX)[0]
        dlist.append(f'{min(dist_arr):.06f}-{idx}')
    # for i in sorted(dlist)[:17]:
    #     print(i)
    # o_idxes = ' '.join([i.split('-')[1] for i in sorted(dlist)[:17]])
    o_idxes = [int(i.split('-')[1]) for i in sorted(dlist)[:17]]
    o_xyz = np.array([xyzs[i] for i in o_idxes])
    dlist = []
    for idx, xyz in enumerate(xyzs[:AT]):
        if idx < O+RU:
            continue
        dist_arr = dar(xyz, o_xyz, box=BOX)[0]
        dlist.append(f'{min(dist_arr):.06f}-{idx}')
        print(dlist[-1])
    h_idxes = ' '.join([i.split('-')[1] for i in sorted(dlist)[:len(o_idxes)*2-6]])
    print(' '.join([str(i) for i in o_idxes]) + ' ' + h_idxes)
    print(h_idxes)
        # print(dlist[-1])
                     # print(dist_arr)



