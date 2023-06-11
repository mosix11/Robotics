import numpy as np
import matplotlib.pyplot as pyplot
from sympy import *
from sympy.physics.control.lti import TransferFunctionMatrix, TransferFunction

import roboticstoolbox as rtb
import DKModels as DiK
import DifferetialKinematiks as DK
import InverseDifferentialKinematics as DIK
import RedundancyResolution as RR
import EnvironmentInteraction as EI
import MatrixInversion as MI
import Control as Control
import DynamicsLE as DLE
import Utils as Utils


def print_latex(expr):
    print(latex(expr), '\n \n \n')

def simp(exp):
    return simplify(factor(exp))

# t = symbols('t', real=True, nonnegative=True)
# q1, q2, q3 = Function('q_1', real=True)(t), Function('q_2', real=True)(t), Function('q_3', real=True)(t)
# q = Matrix([q1, q2, q3])
# v = symbols('v', real=True, positive=True)

# q0 = Matrix([pi/4, 0, 0])
# r = Matrix([
#     cos(q1) + cos(q2) + cos(q3),
#     sin(q1) + sin(q2) + sin(q3)
# ])
# r_y = Matrix([
#     q2
# ])
# r_aux = Matrix.vstack(r, r_y)

# J = DK.get_analytic_jacobian(r, q)
# J_y = DK.get_analytic_jacobian(r_y, q)
# J_aux = DK.get_analytic_jacobian(r_aux, q)



# dr = diff(r, t)
# dr_y = diff(r_y, t)
# dr_aux = diff(r_aux, t)

# r0 = r.evalf(subs=dict(zip(q, q0)))

# r_desired = r0 + Matrix([
#     0, v*t
# ])
# r_y_desired = Matrix([
#     0
# ])

# r_aux_desired = Matrix.vstack(r_desired, r_y_desired)
# dr_aux_desired = diff(r_aux_desired, t)
# # J_aux_desired = DK.get_analytic_jacobian(r_aux_desired, q)

# # dr = J*dq -> dq = J**-1 * dr

# # dq = simp(DIK.get_inverse(J_aux) * dr_aux_desired)


# singularities_J_aux = DK.get_singularities_of_j(J_aux, q)
# det_J_aux = DK.get_det_j(J_aux)
# # print_latex(det_J_aux)

# eq1 = Eq(q1-q3, 0)
# eq2 = Eq(2*cos(q1) - (1+sqrt(2)/2), 0)
# res = solve([eq1, eq2], [q1, q3], domain=Interval(-pi, pi, True, False))
# tmp = acos((sqrt(2)+2)/4).evalf()
# qs = [tmp, 0, tmp]

# Js = J.evalf(subs=dict(zip(q, qs)))
# Js_y = J_y.evalf(subs=dict(zip(q, qs)))
# Js_aux = J_aux.evalf(subs=dict(zip(q, qs)))
# dr_aux_desired = dr_aux_desired.subs(v, 1)

# dq_ps = RR.get_dq_command_jacobian_based(Js_aux, dr_aux_desired)
# dq_dls = RR.get_dq_command_jacobian_based(Js_aux, dr_aux_desired, DLS=True, mu=sqrt(0.25))
# dq_tp = RR.get_dq_command_TP_2task(Js, Js_y, diff(r_desired, t).subs(v, 1), diff(r_y_desired, t))
# print_latex(dq_tp)

#*******************************************************************************************************************************
#*******************************************************************************************************************************

# t = symbols('t', real=True, nonnegative=True)
# q1, q2, q3 = Function('q_1', real=True)(t), Function('q_2', real=True)(t), Function('q_3', real=True)(t)
# q = Matrix([q1, q2, q3])

# r = Matrix([
#     cos(q1) + cos(q1+q2) + cos(q1+q2+q3),
#     sin(q1) + sin(q1+q2) + sin(q1+q2+q3)
# ])
# r_y = Matrix([
#     q1+q2+q3
# ])
# r_aux = Matrix.vstack(r, r_y)

# J = DK.get_analytic_jacobian(r, q)
# J_y = DK.get_analytic_jacobian(r_y, q)
# J_aux = DK.get_analytic_jacobian(r_aux, q)


# qd = Matrix([pi/4, 0, pi/4])
# vd = Matrix([2, -1])
# vd_y = Matrix([0])
# vd_aux = Matrix([vd, vd_y])

# FP_precision = 4
# J_eval = J.evalf(FP_precision, subs=dict(zip(q, qd)), chop=True)
# J_y_eval = J_y.evalf(FP_precision, subs=dict(zip(q, qd)), chop=True)
# J_aux_eval = Matrix.vstack(J_eval, J_y_eval)

# dq = RR.get_dq_command_jacobian_based(J_aux_eval, vd_aux)
# print(dq)
# print(vd_aux - J_aux_eval*dq)

# dq_tp = RR.get_dq_command_TP_2task(J_y_eval, J_eval, vd_y, vd)
# print(dq_tp)
# print(vd_aux - J_aux_eval*dq_tp)
#*******************************************************************************************************************************
#*******************************************************************************************************************************



# t = symbols('t', real=True, nonnegative=True)
# q1, q2 = Function('q_1', real=True)(t), Function('q_2', real=True)(t)
# q = Matrix([q1, q2])

# m11, m12, M22 = Function('m11', real=True)(t), Function('m12', real=True)(t), Function('M22', real=True)(t)

# M = Matrix([
#     [m11, m12],
#     [m12, M22]
# ])

# n = Function('n', real=True)(t)
# taw1, taw2 = symbols('taw_1, taw_2', real=True)
# u = Matrix([taw1, taw2])


# k = symbols('k', real=True, positive=True)
# constraints = Matrix([k, q2])
# h = q - constraints
# A = EI.compute_J_of_constraints_A(h, q)

# print_latex(A)

#*******************************************************************************************************************************
#*******************************************************************************************************************************


# t = symbols('t', real=True, nonnegative=True)
# e = Function('e', real=True)(t)
# de = diff(e, t)
# dde = diff(de, t)
# eq = Eq(dde + 4*de + 4 * e, 0)
# res = dsolve(eq, e, ics={e.subs(t, 0):1.8, de.subs(t, 0):0})
# print(res)
# et = res.rhs
# error = et.subs(t, 1)
# print(error.evalf())
# # print_latex(res)


#*******************************************************************************************************************************
#*******************************************************************************************************************************

# q1, q2 = symbols('q_1, q_2', real=True)
# q = Matrix([q1, q2])
# m1, m2, g0, dc1, l1 = symbols('m_1, m_2, g_0, d_c1, l_1', real=True)

# gq = Matrix([(m1*dc1 + m2*l1)*g0*sin(q1) , 0])
# dgq = gq.jacobian(q)
# print_latex(simp(dgq))
# ind_norm_dgq = dgq.norm(2)
# print_latex(simp(ind_norm_dgq))

#*******************************************************************************************************************************
#*******************************************************************************************************************************

# t = symbols('t', real=True, positive=True)
# s = symbols('s')
# e = Function('e', real=True)(t)

# de = e.diff(t)
# dde = de.diff(t)

# error_dynamics = dde + 100*de + 20*e

# lap_error_dynamics = laplace_transform(error_dynamics, t, s)
# print_latex(lap_error_dynamics)

#*******************************************************************************************************************************
#*******************************************************************************************************************************


# t = symbols('t', real=True, positive=True)
# s = symbols('s')
# q1, q2 = Function('q_1', real=True)(t), Function('q_2', real=True)(t)
# qs1, qs2 = Function('q_s1', real=True)(t), Function('q_s2', real=True)(s)
# q1d, q2d = symbols('q_1d, q2_d', real=True)
# q = Matrix([q1, q2])
# qs = Matrix([qs1, qs2])
# dq = q.diff(t)
# ddq = dq.diff(t)
# qd = Matrix([q1d, q2d])

# m1, m2 = symbols('m_1, m_2', real=True, positive=True)
# M = Matrix([
#     [m1, 0],
#     [0, m2]
# ])
# Ks = symbols('K_s', real=True, positive=True)

# Kp1, Kp2, Kd1, Kd2 = symbols('K_P1, K_P2, K_D1, K_D2', real=True, nonnegative=True)
# Kp = Matrix([
#     [Kp1, 0],
#     [0, Kp2]
# ])
# Kd = Matrix([
#     [Kd1, 0],
#     [0, Kd2]
# ])

# tmp = Matrix([
#     [Ks*(q1-q2) + Kp1*q1],
#     [Ks*(q2-q1) + Kp2*q2]
# ])
# Kp_bar = tmp.jacobian(q)

# closed_loop_eq = simp(M*ddq + Kd*dq + Kp_bar*q - Kp*qd)


# x = Matrix([q1, q2, q1.diff(t), q2.diff(t)])
# dx = x.diff(t)
# A_lower = MI.get_inverse(M) * (Kp*qd - Kd*dq - Kp_bar*q)

# A_upper = Matrix([q1.diff(t), q2.diff(t)])
# A = Matrix.vstack(A_upper, A_lower)
# A = simp(A.jacobian(x))

# # print(A.det())
# # print(A.trace())
# lamda = symbols('lambda')
# A = A - eye(4)*lamda

# # ll = simp(A.trace()**2 - 4*A.det())

# print_latex(roots(simp(A.det()), lamda))

# laplace_cle = laplace_transform(closed_loop_eq, t, s, legacy_matrix=False)[0]

# lap_eq1 = Utils.laplace_simp(laplace_cle[0], q, qs, t, s)
# lap_eq2 = Utils.laplace_simp(laplace_cle[1], q, qs, t, s)
# lap_eq = Matrix([
#     lap_eq1,
#     lap_eq2
# ])

# lhs = (s*lap_eq + Kp*qd).jacobian(qs)
# rhs = Kp*qd
# W = simp(MI.get_inverse(lhs)*rhs)
# W_tf = TransferFunctionMatrix.from_Matrix(W, s)
# W = W[0].as_numer_denom()[1].collect(s**5) #.collect(s).collect(s**2).collect(s**3)
# print_latex(W)
# denom = W[0].as_numer_denom()[1]
# print_latex(roots(denom, s))


#*******************************************************************************************************************************
#*******************************************************************************************************************************


## Pendulum under gravity

# t = symbols('t', real=True, nonnegative=True)
# q = Function('theta', real=True)(t)
# I0, m, d, g0, L = symbols('I_0, m, d_c, g_0, L', real=True, positive=True)

# dq = q.diff(t)
# ddq = dq.diff(t)

# dynamics = (I0 + m*d**2)*ddq + m*d*g0*sin(q)

# dynamics = Matrix([dynamics.subs({m:10, d:1, g0:9.81})])

# Control.iterative_scheme_regulation(dynamics, 500, 45, Matrix([q]), t, Matrix([pi/3]), Matrix([0]))

#*******************************************************************************************************************************
#*******************************************************************************************************************************


## Reduced Dynamics

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

# # PR with P on y axis of world frame
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
# RF0_to_RFw = Matrix([
#     [1, 0, 0],
#     [0, 0, 1],
#     [0, -1, 0],
# ])

# dk_model = DiK.DireckKinematic(DH_table)
# HTM, rot, pos = dk_model.compute_DK(0, 2)
# rot, pos = RF0_to_RFw*rot, RF0_to_RFw*pos
# r = Matrix([pos[0], pos[1]])
# # print_latex(simp(HTM_RF0_to_RFw * HTM))
# # print(DK.get_analytic_jacobian(pos, [q1, q2, q3]))

# dq = q.diff(t)
# ddq = dq.diff(t)
# J = DK.get_analytic_jacobian(r, q)


# u, M, c, g, B, U_v, U_c = DLE.compute_forward_dynamics(DH_table, q, CoM, masses, g_vector, inertia_matrices, None, None, None, None, t)
# # print_latex(M)
# # print_latex(c)
# # print_latex(g)


# Ax, Ay, Bx, By = symbols('A_x, A_y, B_x, B_y', real=True)
# hq = (pos[1] - By)/(Ay - By) - (pos[0] - Bx)/(Ax - Bx)
# hq = Matrix([hq])
# A = hq.jacobian(q)
# D = Matrix([0 , Ay - By]).T
# # A_new = Matrix.vstack(A, D)
# # print(EI.is_sigular(A_new))
# E, F = EI.compute_E_F_from_A_D(A, D)

# pseudo_v = D*dq


# print_latex(E)
# print_latex(F)

#*******************************************************************************************************************************
#*******************************************************************************************************************************