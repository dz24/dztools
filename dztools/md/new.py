# -*- coding: utf-8 -*-
# Copyright (c) 2021, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""The order parameter for the hysteresis example."""
from copy import deepcopy
import logging
import glob
import os
import numpy as np
import MDAnalysis.coordinates.XYZ as md2
from MDAnalysis.analysis import distances

# Assign global variables:
box = [12.4138, 12.4138, 12.4138, 90.0, 90.0, 90.0]
at, o, ru = 193, 64, 2
NAMES = ['Ru']*ru + ['O']*o + ['H']*(at-o-ru)


def calculate(xyz1):
    """Calculate the order parameter.

    Here, the order parameter is just the distance between two
    particles.

    Parameters
    ----------
    system : object like :py:class:`.System`
        The object containing the positions and box used for The
        calculation.

    Returns
    -------
    out : list of floats
        The rate-of-change of the distance order parameter.
   """
    # If key == True, we found the correct HOMO files.
    traj1 = md2.XYZReader(latest_s1, box=box)
    traj2 = md2.XYZReader(latest_s2, box=box)
    OP, dist1, spin, loc, bad, ruel_dist = finder(traj1, traj2, f1, f2)

    # We managed to select the correct X orbital:
    print('loc:', loc, 'spin:', spin, 'OP:', OP, 'dist:', dist1)
    print('f1:', f1, 'out of', len(traj1) - 1, 'file:', latest_s1)
    print('f2:', f2, 'out of', len(traj2) - 1, 'file:', latest_s2)
    print('Ru1 xyz:', p1, 'orb s2-1', traj2[f2][-1])

    OPfake = OP
    OPtrue = OP
    complexCount = 0

    if (OPtrue > 0.990):
        stable, complexCount = checkStable(traj1)
        if(stable == False):
            OPfake = 0.989
    elif (OPtrue < -0.990):
        stable, complexCount = checkStable(traj1)
        if (stable == False):
            OPfake = -0.989

    return [OPfake, OPtrue, complexCount, ruel_dist]
    :
    # We did not manage to select the correct X orbital:
    return [9000]


def checkStable(trajectory):
    """Calculate something ."""
    traj_s1 = trajectory
    L = 12.4138
    dH = 0.91 / 2.
    dO = 1.81 / 2.
    cutoffs = {"H": dH, "O": dO}

    def DISTANCE(c1, c2, L):
        vector = c1 - c2
        if L is not None:
            vector -= L * np.around(vector / L)
        d = np.sqrt(sum(vector * vector))
        return d

    def CONNECTED(at1, at2, cutoffs, L):
        c1, el1 = at1[0], at1[1]
        cutoff1 = cutoffs[el1]
        c2, el2 = at2[0], at2[1]
        cutoff2 = cutoffs[el2]
        d = DISTANCE(c1, c2, L)
        return d < cutoff1 + cutoff2

    def EXTRACTNEIGHBORSFROMLIST(atom, leftover, cutoffs, L):
        indexleftover = 0
        extract = []
        while indexleftover < len(leftover):
            secatom = leftover[indexleftover]
            if CONNECTED(atom, secatom, cutoffs, L):
                extract += [secatom]
                del leftover[indexleftover]
            else:
                indexleftover += 1
        return extract, leftover

    def MOLECLIST(atomlist, L, cutoffs):
        moleclist = []
        leftover = deepcopy(atomlist)
        while len(leftover) > 0:
            mol = []
            mol += [leftover[0]]
            del leftover[0]
            iat = 0
            while iat < len(mol):
                atom = mol[iat]
                neighbors, leftover = EXTRACTNEIGHBORSFROMLIST(atom, leftover,
                                                               cutoffs, L)
                mol += neighbors
                iat += 1
            moleclist += [mol]
        return moleclist

    atomList = []
    # Make atomList
    for i in range(2, 193):
        if (i < 66):
            atomList.append([traj_s1[len(traj_s1) - 1][i], "O"])
        else:
            atomList.append([traj_s1[len(traj_s1) - 1][i], "H"])

    myList = MOLECLIST(atomList, L, cutoffs)

    # Make molecules
    molecList = []
    # molecString = ""
    for i in myList:
        tempMolecule = ""
        for j in range(len(i)):
            tempMolecule += i[j][1]
            # molecString += i[j][1]
        # molecString += " | "
        molecList.append(tempMolecule)

    flag = True
    complexCount = 0
    for i in molecList:
        if (i != "OH" and i != "OHH"):
            flag = False
            complexCount += 1
    return flag, complexCount


def picker(HOMOs, p1, p2):
    """Calculate something ."""
    HOMO_s1 = [s for s in HOMOs if 'HOMO_centers_s1' in s]
    HOMO_s2 = [s for s in HOMOs if 'HOMO_centers_s2' in s]
    f1, f2 = None, None
    if len(HOMO_s1) != len(HOMO_s2):
        exit('len homo_s1 != len homo_s2')

    for k in range(len(HOMO_s1)):
        latest_s1 = max(HOMO_s1, key=os.path.getctime)
        latest_s2 = max(HOMO_s2, key=os.path.getctime)

        s1_key, f1, len1 = picker_a(latest_s1, p1, p2)
        s2_key, f2, len2 = picker_a(latest_s2, p1, p2)

        if all((s1_key, s2_key)) or not (HOMO_s1 and HOMO_s2):
            print('bongo', s1_key, s2_key,
                  latest_s1, len1,
                  latest_s2, len2)
            break

        print('not in vvv')
        print(latest_s1)
        print(latest_s2)
        HOMO_s1.remove(latest_s1)
        HOMO_s2.remove(latest_s2)
        print('not in ^^^')

    if all((s1_key, s2_key)):
        return f1, f2, latest_s1, latest_s2, True
    else:
        return None, None, None, None, False


def picker_a(ltst_s0, p1, p2):
    """Calculate something ."""
    traj_s0 = md2.XYZReader(ltst_s0)

    Found, idx = False, -1
    for i in list(range(len(traj_s0))):
        if np.max(np.abs(np.float32(p1) - traj_s0[i][0])) < 10e-7 and \
                np.max(np.abs(np.float32(p2) - traj_s0[i][1])) < 10e-7:
            Found = True
            idx = i
    if Found:
        return True, idx, len(traj_s0)
    return False, None, None



def finder(traj1, traj2, f1, f2):
    """Calculate something ."""
    at_ar = [0] * (ru+o)
    dic = {i: [] for i in range(ru+o)}
    bad = False

    # Generate atom_array and dict:
    at_ar, dic = finder_a(at_ar, dic, traj1, f1, 's1')
    at_ar, dic = finder_a(at_ar, dic, traj2, f2, 's2')

    # Find the atom a_idx and excess electron x_idx:
    if 6 in at_ar[0:2]:
        # Excess electron exists in one of the Ru atoms
        a_idx = 0 if at_ar[0] == 6 else 1
        x_idx = np.argmax([i['dist'] for i in dic[a_idx]])
    elif (at_ar[0] == at_ar[1] == 5):
        # Excess electron exists in one of the oxygen atoms.
        # Get list for oxygen with max count(x):
        o_l = [i for i, j in enumerate(at_ar) if j == np.max(at_ar)]
        a_idx0, x_idx0, dist = 0, 0, 0
        for a_idx0 in o_l:
            x_idx0 = np.argmax([i['dist'] for i in dic[a_idx0]])
            if dist < dic[a_idx0][x_idx0]['dist']:
                dist = dic[a_idx0][x_idx0]['dist']
                a_idx, x_idx = a_idx0, x_idx0
    else:
        print('It did not seem like we converged..')
        bad = True

    if not bad:
        loc = dic[a_idx][x_idx]['x_loc']
        spin = dic[a_idx][x_idx]['spin']
        traj = traj1 if spin == 's1' else traj2
        f = f1 if spin == 's1' else f2

        # Ru1-el, Ru2-el, Ru1-Ru2
        dist1 = distances.distance_array(traj[f][0], traj[f][loc], box=box)
        dist2 = distances.distance_array(traj[f][1], traj[f][loc], box=box)
        dist3 = distances.distance_array(traj[f][0], traj[f][1], box=box)
        OP = (dist1[0][0]-dist2[0][0])/dist3[0][0]

        return OP, dist1[0][0], spin, loc, bad, dist1[0][0]
    else:
        return None, None, None, None, bad, None


def finder_a(atom_array, dic, traj0, f0, spin):
    """Calculate something ."""
    for i in range(len(traj0[f0][at:][:, 0])):
        dist_arr = distances.distance_array(traj0[f0][at+i],
                                            traj0[f0][:ru+o],
                                            box=box)[0]
        loc = np.argmin(dist_arr)
        atom_array[loc] += 1
        dic[loc].append({'x_loc': i+at,
                         'dist': dist_arr[loc],
                         'spin': spin,
                         'a_loc': loc})
    return atom_array, dic
