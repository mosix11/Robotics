import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from modules import ElemetaryRotations as ER

from sympy import *

def EULER_to_rotmat(phi, theta, psi, seq):
    rotmat = eye(3)
    angles = [phi, theta, psi]
    for i, axis in enumerate(seq):
        rotmat = ER.do_elem_rot(rotmat, axis, angles[i], post_multi=True)
    return trigsimp(factor(rotmat))
    
def rotmat_to_EULER(rotmat, seq):
    R11, R12, R13 = rotmat[0, :]
    R21, R22, R23 = rotmat[1, :]
    R31, R32, R33 = rotmat[2, :]
    
    phi, theta, psi = symbols('phi, theta, psi', real=True)
    euler_rotmat = EULER_to_rotmat(phi, theta, psi, seq)

    
    eq11 = R11 - euler_rotmat[0,0]
    eq12 = R12 - euler_rotmat[0,1]
    eq13 = R13 - euler_rotmat[0,2]
    eq21 = R21 - euler_rotmat[1,0]
    eq22 = R22 - euler_rotmat[1,1]
    eq23 = R23 - euler_rotmat[1,2]
    eq31 = R31 - euler_rotmat[2,0]
    eq32 = R32 - euler_rotmat[2,1]
    eq33 = R33 - euler_rotmat[2,2]


    res = solve([eq11, eq12, eq13, eq21, eq22, eq23, eq31, eq32, eq33],
                [phi, theta, psi],
                domain=Interval(-pi, pi, True, False))

    angles = []
    for tup in res:
        angles.append((tup[0].evalf(), tup[1].evalf(), tup[2].evalf()))
    if len(angles) == 0:
        return 'Singularity'
    else:
        return angles
    
    
def RPY_to_rotmat(phi, theta, psi, seq):
    rotmat = eye(3)
    angles = [psi, theta, phi] # reverse angles
    seq_euler = seq[::-1] # reverse axes
    for i, axis in enumerate(seq_euler):
        rotmat = ER.do_elem_rot(rotmat, axis, angles[i], post_multi=True)
    return trigsimp(factor(rotmat))

def rotmat_to_RPY(rotmat, seq):
    angles = rotmat_to_EULER(rotmat, seq[::-1])
    if angles == 'Singularity':
        return 'Singularity'
    else:
        for i, tup in enumerate(angles):
            angles[i] = tup[::-1]
        return angles


# phi, theta, psi = symbols('phi, theta, psi')
# # print(latex(EULER_to_rotmat(phi, theta, psi, 'ZXZ')))
# # print(latex(EULER_to_rotmat(phi, pi, psi, 'ZXZ')))

# rotmat = EULER_to_rotmat(pi/3, pi/4, pi/6, 'ZYZ')
# angles = rotmat_to_EULER(rotmat, 'ZYZ')
# print(latex(angles))

# print(latex(RPY_to_rotmat(psi, theta, phi, 'XYZ')))

# angles = rotmat_to_RPY(RPY_to_rotmat(pi/3, pi/4, pi/6, 'XYZ'), 'XYZ')
# print(latex(angles))