import numpy as np
import matplotlib.pyplot as pyplot
from sympy import *

import roboticstoolbox as rtb
import DKModels as DiK
import DifferetialKinematiks as DK
import InverseDifferentialKinematics as DIK
import RedundancyResolution as RR

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


t = symbols('t', real=True, nonnegative=True)
q1, q2, q3 = Function('q_1', real=True)(t), Function('q_2', real=True)(t), Function('q_3', real=True)(t)
q = Matrix([q1, q2, q3])

r = Matrix([
    cos(q1) + cos(q1+q2) + cos(q1+q2+q3),
    sin(q1) + sin(q1+q2) + sin(q1+q2+q3)
])
r_y = Matrix([
    q1+q2+q3
])
r_aux = Matrix.vstack(r, r_y)

J = DK.get_analytic_jacobian(r, q)
J_y = DK.get_analytic_jacobian(r_y, q)
J_aux = DK.get_analytic_jacobian(r_aux, q)


qd = Matrix([pi/4, 0, pi/4])
vd = Matrix([2, -1])
vd_y = Matrix([0])
vd_aux = Matrix([vd, vd_y])

FP_precision = 4
J_eval = J.evalf(FP_precision, subs=dict(zip(q, qd)), chop=True)
J_y_eval = J_y.evalf(FP_precision, subs=dict(zip(q, qd)), chop=True)
J_aux_eval = Matrix.vstack(J_eval, J_y_eval)

# dq = RR.get_dq_command_jacobian_based(J_aux_eval, vd_aux)
# print(dq)
# print(vd_aux - J_aux_eval*dq)

dq_tp = RR.get_dq_command_TP_2task(J_y_eval, J_eval, vd_y, vd)
print(dq_tp)
print(vd_aux - J_aux_eval*dq_tp)











