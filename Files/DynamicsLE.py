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


def get_CoM_velocities(dk_model, CoM, t):
    omega_velocities = []
    linear_velocities = []
    DOF = dk_model.get_DOF()
    for k in range(0, DOF):
        _, rot, pos = dk_model.compute_DK(0, k+1)

        omega = simp(DK.get_omega_from_rotmat(rot, t))
        CoM_RF0 = simp(rot * CoM[:,k])
        lin_vel = simp((diff(pos, t) + Utils.cross_product(omega, CoM_RF0)))
        omega = simp(rot.T * omega)
        lin_vel = simp(rot.T * lin_vel)
        omega_velocities.append(omega)
        linear_velocities.append(lin_vel)

    return omega_velocities, linear_velocities

def get_kinetic_energy(dk_model, CoM, masses, inertia_matrices, t):
    DOF = dk_model.get_DOF()
    omega, l_vel = get_CoM_velocities(dk_model, CoM, t)
    kinetic_energies = []
    
    for i in range(0, DOF):
        T = Rational(1, 2) * masses[i] * l_vel[i].T * l_vel[i] + Rational(1, 2) * omega[i].T * inertia_matrices[i] * omega[i]
        T = simp(T)
        kinetic_energies.append(T)
    
    total_T = kinetic_energies[0]
    for k in range(1, DOF):
        total_T += kinetic_energies[k]
        total_T = simp(total_T)

    return total_T

def get_potential_energy(dk_model, CoM, masses, g_vector, t):
    U_i = []
    DOF = dk_model.get_DOF()
    
    for k in range(0, DOF):
        _, rot, pos = dk_model.compute_DK(0, k+1)
        CoM_RF0 = simp(rot * CoM[:,k])
        CoM_RF0 = CoM_RF0 + pos
        U = simp(-masses[k] * g_vector.T * CoM_RF0)
        U_i.append(U)
    
    total_U = U_i[0]
    for k in range(1, DOF):
        total_U += U_i[k]
        total_U = simp(total_U)

    return total_U


def factorize_M_from_kinetic_energy(T, q_vector, t):
    dq = diff(q_vector, t)
    M = T.jacobian(dq).jacobian(dq)
    M = simp(M)
    return M
    
def get_christoffel_matrices(M, q_vector, t):
    DOF, _ = shape(q_vector)
    dq_vector = diff(q_vector, t)
    C = []
    for k in range(DOF):
        C_k = Rational(1, 2) * (M[:,k].jacobian(q_vector) + (M[:,k].jacobian(q_vector)).T - diff(M, q_vector[k]))
        C_k = simp(C_k)
        C.append(C_k)
    
    c = []
    for C_k in C:
        c_k = dq_vector.T * C_k * dq_vector
        c_k = simp(c_k)
        c.append(c_k)
    
    c = Matrix(c)
    return c, C

def factorize_B_from_motor_kinetic_energy(DOF, motor_inertia_matrices, q_vector, nr, t):
    dq = diff(q_vector, t)
    T_m_matx = []
    for k in range (DOF):
        T_m_i = Rational(1, 2) * motor_inertia_matrices[k] * nr[k] * dq[k]**2
        T_m_i = simp(T_m_i)
        T_m_matx[k] = T_m_i
        
    T_m = T_m_matx[0]
    for k in range(1, DOF):
        T_m += T_m_matx[k]
        T_m = simp(T_m)
    
    B = T_m.jacobian(dq).jacobian(dq)
    B = simp(B)
    
    return B
    

def compute_forward_dynamics(DH_table, q_vector, CoM, masses, g_vector, inertia_matrices, motor_inertia_matrices, nr, F_v, F_c, t):
    dk_model = DiK.DireckKinematic(DH_table)
    
    # mat, rot, pos = dk_model.compute_DK(0,3)
    # RB = Matrix([
    #     [1, 0, 0],
    #     [0, 0, 1],
    #     [0, 1, 0]
    # ])
    # print_latex(pos)
    # return 0
    
    DOF = dk_model.get_DOF()

    dq = diff(q_vector, t)
    ddq = diff(dq, t) 

    U = get_potential_energy(dk_model, CoM, masses, g_vector, t)
    T = get_kinetic_energy(dk_model, CoM, masses, inertia_matrices, t)
    
    M = factorize_M_from_kinetic_energy(T, q_vector, t)
    c, C = get_christoffel_matrices(M, q_vector, t)
    g = simp((U.jacobian(q_vector)).T)
    
   
    
    B = None
    U_v = None
    U_c = None
    
    
    if motor_inertia_matrices != None and nr != None:
        B = factorize_B_from_motor_kinetic_energy(DOF, motor_inertia_matrices, q_vector, nr, t)
             
    if F_v != None:
        U_v = -F_v * dq
    
    if F_c != None:
        U_c = -F_c * sign(dq)
    
    

    
    dynamics = simp(((M + B) if B != None else M) * ddq + c + g)
    if U_v != None : dynamics -= U_v
    if U_c != None : dynamics -= U_c
    dynamics = simp(dynamics)
    
    return dynamics, M, c, g, B, U_v, U_c
    
    
    
def transform_dynamic_model_generalized_coordinate(p, q, M, c, g, n, t):
    J = p.jacobian(q)
    J_inv = IDK.get_inverse(J)
    dJ = diff(J, t)
    dp = diff(p, t)
    ddp = diff(dp, t)
    Mp = simp(J_inv.T * M * J_inv)
    cp = simp(J_inv.T * c - Mp * dJ * J_inv * dp)
    gp = simp(J_inv.T * g)
    
    up = simp(J_inv.T * M * J_inv * ddp + J_inv.T * (n - M * J_inv * dJ * J_inv * dp))
    
    
    return up, Mp, cp, gp

# t = symbols('t', real=True, nonnegative=True)

# q1, q2 = Function('q_1', real=True)(t), Function('q_2', real=True)(t)
# l1, l2 = symbols('l_1, l_2', real=True, positive=True)
# m1, m2 = symbols('m_1, m_2', real=True, positive=True)
# d1, d2 = symbols('d_1, d_2', real=True, positive=True)
# gravity = symbols('g', real=True, nonnegative=True)

# q_vector = Matrix([q1, q2])
# DH_table = [
#     [0, l1, q1, 0],
#     [0, l2, q2, 0]
# ]

# DOF = 2

# masses = [m1, m2]
# g_vector = Matrix([0, -gravity, 0])
# inertia_matrices = []
# # for i in range(DOF):
# #     I = MatrixSymbol('I_' + str(i+1), 3, 3)
# #     inertia_matrices.append(I)
# for i in range(DOF):
#     I = Matrix([
#         [symbols('I_' + str(i+1) + '.11', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.12', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.13', real=True, nonnegative=True)],
#         [symbols('I_' + str(i+1) + '.21', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.22', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.23', real=True, nonnegative=True)],
#         [symbols('I_' + str(i+1) + '.31', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.32', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.33', real=True, nonnegative=True)]
#     ])
#     inertia_matrices.append(I)



# CoM = Matrix([
#     [-l1+d1, -l2+d2],
#     [0, 0],
#     [0, 0]
# ])

# u, M, c, g, B, U_v, U_c = compute_forward_dynamics(DH_table, q_vector, CoM, masses, g_vector, inertia_matrices, None, None, None, None, t)

# p = Matrix([
#     q1, q1+q2
# ])
# up, Mp, cp, gp = transform_dynamic_model_generalized_coordinate(p, q_vector, M, c, g, c+g, t)

# t1 = u.evalf(subs=dict(zip(q_vector, [pi/4, pi/4])))
# t2 = up.evalf(subs=dict(zip(q_vector, [pi/4, pi/2])))
# t1 = t1.evalf(subs=dict(zip([m1, m2, l1, l2, d1, d2], [1, 1, 1, 1, 0.5, 0.5])))
# t2 = t2.evalf(subs=dict(zip([m1, m2, l1, l2, d1, d2], [1, 1, 1, 1, 0.5, 0.5])))

# print(simp(t1))
# print(simp(t2))

# print_latex(M)
# print_latex(Mp)


#*******************************************************************************************************************************
#*******************************************************************************************************************************

# t = symbols('t', real=True, nonnegative=True)
# q1, q2, q3 = Function('q_1', real=True, nonnegative=True)(t), Function('q_2', real=True, nonnegative=True)(t), Function('q_3', real=True, nonnegative=True)(t)
# l1, l2, l3 = symbols('l_1, l_2, l_3', real=True, positive=True)
# m1, m2, m3 = symbols('m_1, m_2, m_3', real=True, positive=True)
# d1, d2, d3 = symbols('d_1, d_2, d_3', real=True, positive=True)
# gravity = symbols('g', real=True, nonnegative=True)

# q_vector = Matrix([q1, q2, q3])

# Dh_table = [
#     [0, l1, q1, 0],
#     [l2, 0, q2, pi/2],
#     [0, l3, q3, 0]
# ]

# DOF = 3

# masses = [m1, m2, m3]
# g_vector = Matrix([0, -gravity, 0])
# inertia_matrices = []

# for i in range(DOF):
#     I = Matrix([
#         [symbols('I_' + str(i+1) + '.11', real=True, nonnegative=True), 0, 0],
#         [0, symbols('I_' + str(i+1) + '.22', real=True, nonnegative=True), 0],
#         [0, 0, symbols('I_' + str(i+1) + '.33', real=True, nonnegative=True)]
#     ])
#     inertia_matrices.append(I)
    
# CoM = Matrix([
#     [-l1+d1, 0, -l3+d3],
#     [0, -l2+d2, 0],
#     [0, 0, 0]
# ])

# u, M, c, g, B, U_v, U_c = compute_forward_dynamics(Dh_table, q_vector, CoM, masses, g_vector, inertia_matrices, None, None, None, None, t)

# print_latex(M)


#*******************************************************************************************************************************
#*******************************************************************************************************************************

# t = symbols('t', real=True, nonnegative=True)

# q1, q2 = Function('q_1', real=True, nonnegative=True)(t), Function('q_2', real=True, nonnegative=True)(t)
# l1, l2 = symbols('l_1, l_2', real=True, positive=True)
# m1, m2 = symbols('m_1, m_2', real=True, positive=True)
# gravity = symbols('g', real=True, nonnegative=True)

# rc_1_x, rc_1_y, rc_2_x, rc_2_y = symbols('rc_1x, rc_1y, rc_2x, rc_2y', real=True)

# q_vector = Matrix([q1, q2])
# DH_table = [
#     [0, l1, q1, 0],
#     [0, l2, q2, 0]
# ]

# DOF = 2

# masses = [m1, m2]
# g_vector = Matrix([0, -gravity, 0])
# inertia_matrices = []
# # for i in range(DOF):
# #     I = MatrixSymbol('I_' + str(i+1), 3, 3)
# #     inertia_matrices.append(I)
# for i in range(DOF):
#     I = Matrix([
#         [symbols('I_' + str(i+1) + '.11', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.12', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.13', real=True, nonnegative=True)],
#         [symbols('I_' + str(i+1) + '.21', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.22', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.23', real=True, nonnegative=True)],
#         [symbols('I_' + str(i+1) + '.31', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.32', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.33', real=True, nonnegative=True)]
#     ])
#     inertia_matrices.append(I)



# CoM = Matrix([
#     [rc_1_x, rc_2_x],
#     [rc_1_y, rc_2_y],
#     [0, 0]
# ])

# u, M, c, g, B, U_v, U_c = compute_forward_dynamics(DH_table, q_vector, CoM, masses, g_vector, inertia_matrices, None, None, None, None, t)

# # m11 = M[0,0]
# # m11 -= 2 * m2 * l1 * (l2 + rc_2_x) * cos(q2)
# # m11 -= 2 * -m2 * l1 * rc_2_y * sin(q2)
# # m11 -= inertia_matrices[0][2,2] + inertia_matrices[1][2,2] + m1*((l1 + rc_1_x)**2 + rc_1_y**2) + m2*l1**2 + m2*((l2 + rc_2_x)**2 + rc_2_y**2)
# # m11 = simp(m11)
# print_latex(M)
# print_latex(c)
# print_latex(g)

#*******************************************************************************************************************************
#*******************************************************************************************************************************


# t = symbols('t', real=True, nonnegative=True)
# q1, q2, q3 = Function('q_1', real=True)(t), Function('q_2', real=True)(t), Function('q_3', real=True)(t)
# l3 = symbols('l_3', real=True, positive=True)
# m1, m2, m3 = symbols('m_1, m_2, m_3', real=True, positive=True)
# d1, d2, d3 = symbols('d_1, d_2, d_3', real=True, positive=True)
# gravity = symbols('g', real=True, nonnegative=True)

# q_vector = Matrix([q1, q2, q3])

# Dh_table = [
#     [q1, 0, pi/2, pi/2],
#     [q2, 0, -pi/2, pi/2],
#     [0, l3, q3+pi/2, 0] 
# ]

# DOF = 3

# masses = [m1, m2, m3]
# g_vector = Matrix([0, -gravity, 0])
# inertia_matrices = []

# for i in range(DOF):
#     I = Matrix([
#         [symbols('I_' + str(i+1) + '.11', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.12', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.13', real=True, nonnegative=True)],
#         [symbols('I_' + str(i+1) + '.21', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.22', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.23', real=True, nonnegative=True)],
#         [symbols('I_' + str(i+1) + '.31', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.32', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.33', real=True, nonnegative=True)]
#     ])
#     inertia_matrices.append(I)
    
# CoM = Matrix([
#     [0, 0, d3-l3],
#     [d1-q1, d2-q2, 0],
#     [0, 0, 0]
# ])

# # compute_forward_dynamics(Dh_table, q_vector, CoM, masses, g_vector, inertia_matrices, None, None, None, None, t)
# u, M, c, g, B, U_v, U_c = compute_forward_dynamics(Dh_table, q_vector, CoM, masses, g_vector, inertia_matrices, None, None, None, None, t)
# p = Matrix([
#     q2+l3*cos(q3), q1+l3*sin(q3), q3
# ])
# up, Mp, cp, gp = transform_dynamic_model_generalized_coordinate(p, q_vector, M, c, g, c+g, t)

# print_latex(Mp)








#*******************************************************************************************************************************
#*******************************************************************************************************************************

# t = symbols('t', real=True, nonnegative=True)
# DOF = 3
# q1, q2, q3 = Function('q_1', real=True)(t), Function('q_2', real=True)(t), Function('q_3', real=True)(t)
# q = Matrix([q1, q2, q3])
# l1, l2, l3 = symbols('l_1, l_2, l_3', real=True, positive=True)
# m1, m2, m3 = symbols('m_1, m_2, m_3', real=True, positive=True)
# d1, d2, d3 = symbols('d_1, d_2, d_3', real=True, positive=True)
# gravity = symbols('g', real=True, nonnegative=True)

# masses = [m1, m2, m3]
# g_vector = Matrix([0, -gravity, 0])
# # g_vector = Matrix([-gravity, 0, 0])


# inertia_matrices = []

# for i in range(DOF):
#     I = Matrix([
#         [symbols('I_' + str(i+1) + '.11', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.12', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.13', real=True, nonnegative=True)],
#         [symbols('I_' + str(i+1) + '.21', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.22', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.23', real=True, nonnegative=True)],
#         [symbols('I_' + str(i+1) + '.31', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.32', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.33', real=True, nonnegative=True)]
#     ])
#     inertia_matrices.append(I)

# CoM = Matrix([
#     [d1, 0, d3],
#     [0, d2, 0],
#     [0, 0, 0]
# ])


#*******************************************************************************************************************************
#*******************************************************************************************************************************


# t = symbols('t', real=True, nonnegative=True)
# DOF = 2
# q1, q2 = Function('q_1', real=True)(t), Function('q_2', real=True)(t)
# q = Matrix([q1, q2])
# l1, l2 = symbols('l_1, l_2', real=True, positive=True)
# m1, m2 = symbols('m_1, m_2', real=True, positive=True)
# d1, d2 = symbols('d_1, d_2', real=True, positive=True)
# gravity = symbols('g', real=True, nonnegative=True)

# masses = [m1, m2]
# # g_vector = Matrix([0, -gravity, 0])
# # g_vector = Matrix([gravity, 0, 0])
# g_vector = Matrix([0, 0, -gravity])

# inertia_matrices = []

# for i in range(DOF):
#     I = Matrix([
#         [symbols('I_' + str(i+1) + '.11', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.12', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.13', real=True, nonnegative=True)],
#         [symbols('I_' + str(i+1) + '.21', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.22', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.23', real=True, nonnegative=True)],
#         [symbols('I_' + str(i+1) + '.31', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.32', real=True, nonnegative=True), symbols('I_' + str(i+1) + '.33', real=True, nonnegative=True)]
#     ])
#     inertia_matrices.append(I)

# CoM = Matrix([
#     [0, -l2+d2],
#     [-l1+d1, 0],
#     [0, 0]
# ])


#*******************************************************************************************************************************
#*******************************************************************************************************************************


## R
# DH_table = [
#     [0, l1, q1, 0]
# ]

## RR
# DH_table = [
#     [0, l1, q1, 0],
#     [0, l2, q2, 0]
# ]

## RRR
# DH_table = [
#     [0, l1, q1, 0],
#     [0, l2, q2, 0],
#     [0, l3, q3, 0]
# ]

# ## PR with P on x axis of world frame
# DH_table = [
#     [q1, 0, pi, pi/2],
#     [0, l2, pi/2+q2, 0]
# ]
# HTM_RF0_to_RFw = Matrix([
#     [0, 0, 1, 0],
#     [1, 0, 0, 0],
#     [0, 1, 0, 0],
#     [0, 0, 0, 1]
# ])


# PR with P on y axis of world frame
# DH_table = [
#     [q1, 0, 0, pi/2],
#     [0, l2, q2, 0]
# ]
# HTM_RF0_to_RFw = Matrix([
#     [1, 0, 0, 0],
#     [0, 0, 1, 0],
#     [0, -1, 0, 0],
#     [0, 0, 0, 1]
# ])

## RP with pi/2 angle between links
# DH_table = [
#     [0, l1, q1, pi/2],
#     [q2, 0, 0, 0],
# ]
# HTM_RF0_to_RFw = Matrix([
#     [1, 0, 0, 0],
#     [0, 1, 0, 0],
#     [0, 0, 1, 0],
#     [0, 0, 0, 1]
# ])


## RP with 0 angle between links
# DH_table = [
#     # [0, 0, q1+pi/2, pi/2],
#     # [q2+l1, 0, 0, 0],
#     [0, 0, q1+pi/2, pi/2],
#     [q2, 0, 0, 0],
# ]

# HTM_RF0_to_RFw = Matrix([
#     [1, 0, 0, 0],
#     [0, 1, 0, 0],
#     [0, 0, 1, 0],
#     [0, 0, 0, 1]
# ])

## RRP
# DH_table = [
#     [0, l1, q1, 0],
#     [0, l2, q2, pi/2],
#     [q3, 0, 0, 0],
# ]
# HTM_RF0_to_RFw = Matrix([
#     [1, 0, 0, 0],
#     [0, 1, 0, 0],
#     [0, 0, 1, 0],
#     [0, 0, 0, 1]
# ])


## RPR
# DH_table = [
#     [0, l1, q1, -pi/2],
#     [q2, 0, 0, pi/2], 
#     [0, l3, q3, 0],
# ]
# HTM_RF0_to_RFw = Matrix([
#     [1, 0, 0, 0],
#     [0, 1, 0, 0],
#     [0, 0, 1, 0],
#     [0, 0, 0, 1]
# ])

## RPP


## PRR
# with P on x axis of world frame
# DH_table = [
#     [q1, 0, pi, pi/2],
#     [0, l2, pi/2+q2, 0],
#     [0, l3, q3, 0]
# ]
# HTM_RF0_to_RFw = Matrix([
#     [0, 0, 1, 0],
#     [1, 0, 0, 0],
#     [0, 1, 0, 0],
#     [0, 0, 0, 1]
# ])

## PRR
# with P on y axis of world frame
# DH_table = [
#     [q1, 0, 0, pi/2],
#     [0, l2, q2, 0],
#     [0, l3, q3, 0]
# ]
# HTM_RF0_to_RFw = Matrix([
#     [1, 0, 0, 0],
#     [0, 0, 1, 0],
#     [0, 1, 0, 0],
#     [0, 0, 0, 1]
# ])

## PPR
# DH_table = [
#     [q1, 0, pi/2, pi/2],
#     [q2, 0, pi/2, pi/2],
#     [0, l3, q3, 0]
# ]
# HTM_RF0_to_RFw = Matrix([
#     [0, 0, 1, 0],
#     [1, 0, 0, 0],
#     [0, 1, 0, 0],
#     [0, 0, 0, 1]
# ])

## PRP
# DH_table = [
#     [q1, 0, pi, pi/2],
#     [0, 0, q2+pi, pi/2],
#     [q3, 0, 0, 0]
# ]
# HTM_RF0_to_RFw = Matrix([
#     [0, 0, 1, 0],
#     [1, 0, 0, 0],
#     [0, 1, 0, 0],
#     [0, 0, 0, 1]
# ])



# dk_model = DiK.DireckKinematic(DH_table)
# HTM, rot, pos = dk_model.compute_DK(0, 2)

# print_latex(simp(HTM_RF0_to_RFw * HTM))
# # print(DK.get_analytic_jacobian(pos, [q1, q2, q3]))


# u, M, c, g, B, U_v, U_c = compute_forward_dynamics(DH_table, q, CoM, masses, g_vector, inertia_matrices, None, None, None, None, t)
# print_latex(M)
# print_latex(c)
# print_latex(g)
# print_latex(g.jacobian(q).norm(2))
