import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import DKModels as DiK
import DifferetialKinematiks as DK
import InverseDifferentialKinematics as IDK
import Utils as Utils
from sympy import *

def print_latex(x):
    print(latex(x), '\n \n \n')

def simp(exp):
    return simplify(factor(exp))

def get_dq_command_jacobian_based(J, dr, W=None, DLS=False, mu=None):
    weighted = True if W != None else False
    
    if weighted:
        if DLS:
            print('not implemented')
            # not implemented yet
            pass
        else:
            J_inv = IDK.get_weithed_pinv(J, W)
    else:
        if DLS:
            J_inv = IDK.get_J_DLS(J, mu)
        else:
            J_inv = IDK.get_inverse(J)
    
    dq = J_inv * dr
    
    return simp(dq)


def get_ddq_command_jacobian_based(J, q, dq, ddr, t, W=None, DLS=False, mu=None):
    
    weighted = True if W != None else False
    
    if weighted:
        if DLS:
            print('not implemented')
            # not implemented yet
            pass
        else:
            J_inv = IDK.get_weithed_pinv(J, W)
    else:
        if DLS:
            J_inv = IDK.get_J_DLS(J, mu)
        else:
            J_inv = IDK.get_inverse(J)
    
    dJ = diff(J, t)
    ddq = J_inv * (ddr - dJ*dq) 
    
    return simp(ddq)



def get_dq_command_nullspace(J, dr, dq_0, W=None):
    (m, n) = shape(J)
    weighted = True if W != None else False
    
    if weighted:
        J_inv = IDK.get_weithed_pinv(J,W)
    else:
        J_inv = IDK.get_inverse(J)
    
    dq = J_inv*dr + (eye(m) - J_inv*J)*dq_0
    
    return simp(dq)


def get_ddq_command_nullspace(J, q, dq, ddr, t, ddq_0, W=None):
    (m, n) = shape(J)
    weighted = True if W != None else False
    
    if weighted:
        J_inv = IDK.get_weithed_pinv(J,W)
    else:
        J_inv = IDK.get_inverse(J)

    
    dJ = diff(J, t)
    
    ddq = J_inv * (ddr - dJ*dq) + (eye(m) - J_inv*J)*ddq_0
    
    return simp(ddq)


def get_projector(J):
    J_sharp = IDK.get_inverse(J)
    m, n = shape(J_sharp)
    P = eye(m) - J_sharp * J
    return simp(P)

def get_dq_command_TP_2task(J1, J2, dr1, dr2, q2_0=0):
    J1_sharp = IDK.get_inverse(J1)
    J2_sharp = IDK.get_inverse(J2)
    P1 = get_projector(J1)
    P2 = get_projector(J2)
    J2P1 = Utils.filter_matrix_of_small_values(J2*P1)
    dq = J1_sharp*dr1 + IDK.get_inverse(J2P1) * (dr2 - J2 * J1_sharp * dr1)
    if q2_0 != 0:
        dq += P1 * get_projector(J2P1) * q2_0
    return simp(dq)

# def check_if_argorithmic_sigularity_exist(J, J_y):
#     R_JT = DK.get_image(J.T)
#     R_J_yT = DK.get_image(J_y.T)
#     # print(type(R_JT))
#     # print(type(R_J_yT))
#     # print_latex(R_JT)
#     # print_latex(R_J_yT)
#     for elem1 in R_JT:
#         for elem2 in R_J_yT:
#             intersect = Intersection(elem1, elem2)
#             print_latex(intersect)
#     # intersection = Intersection(R_JT, R_J_yT)
#     # print_latex(R_JT)

#*******************************************************************************************************************************
#*******************************************************************************************************************************

## Dynamic Redundancy Resolution

def get_tau_command_MTN(J, q, dq, ddr, t, M, c, g): #minimum torque norm
    (m, _) = shape(J)
    n = simp(c+g)
    
    dJ = diff(J, t)
    
    tau = IDK.get_inverse(J*(M**-1)) * ((ddr - dJ*dq) - J*(M**-1)*n)
    return simp(tau)


def get_tau_command_MTN_SIIW(J, q, dq, ddr, t, M, c, g): #minimum (squared inverse inertia weighted) torque norm
    (m, _) = shape(J)
    n = simp(c+g)
    
    dJ = diff(J, t)
    
    tau = M * IDK.get_inverse(J) * ((ddr - dJ*dq) - J*(M**-1)*n)
    return simp(tau)

def get_tau_command_MTN_IIW(J, q, dq, ddr, t, M, c, g): #minimum (inverse inertia weighted) torque norm
    (m, _) = shape(J)
    n = simp(c+g)
    
    dJ = diff(J, t)
    
    tau = J.T * IDK.get_inverse(J*(M**-1)*J.T) * ((ddr - dJ*dq) - J*(M**-1)*n)
    return simp(tau)
    



#*******************************************************************************************************************************
#*******************************************************************************************************************************


