import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import *

import MinimalOrientation as MO




def diff_variable_elem_rot(rotmat, t):
    res = diff(rotmat, t)
    return simplify(factor(res))

def get_omega_from_rotmat(rotmat, t):
    d_rotmat = diff_variable_elem_rot(rotmat, t)
    sw = d_rotmat*rotmat.T
    w = Matrix([-sw[1,2], sw[0,2], -sw[0,1]])
    return simplify(factor(w))


def find_d_rotmat_with_omega(rotmat, w):
    sw = Matrix([
        [0, -w[2], w[1]],
        [w[2], 0, -w[0]],
        [-w[1], w[0], 0]
    ])
    d_rotmat = sw*rotmat
    return simplify(factor(d_rotmat))


def find_omega_from_EULER_angles(phi, theta, psi, t, seq):
    rotmat = MO.EULER_to_rotmat(phi, theta, psi, seq)
    w = get_omega_from_rotmat(rotmat, t)
    return simplify(factor(w))

def find_omega_from_RPY_angles(phi, theta, psi, t, seq):
    rotmat = MO.RPY_to_rotmat(phi, theta, psi, seq)
    w = get_omega_from_rotmat(rotmat, t)
    return simplify(factor(w))    


def get_analytic_jacobian(r, variables):
    j_r = r.jacobian(variables)
    return simplify(factor(j_r))


def get_det_j(j):
    (m, n) = shape(j)
    if m == n:
        det_j = j.det()
    elif m<n:
        det_j = (j*j.T).det()
    elif m>n:
        det_j = (j.T*j).det()
    return simplify(factor(det_j))

def get_singularities_of_j(mat, variables):
    det = get_det_j(mat)
    eq = Eq(det, 0)
    res = solve(eq, variables, dict=True, domain=Interval(-pi, pi, True, False))
    return res


def get_nullspace(mat):
    null_space = mat.nullspace()

    for i, v in enumerate(null_space):
        null_space[i] = simplify(factor(v))
    return null_space

def get_columnspace(mat):
    cs = mat.columnspace()
    for i, v in enumerate(cs):
        cs[i] = simplify(factor(v))
    return cs

def get_image(mat):
    rref_mat, idx = mat.T.rref()
    im_mat = []
    for i in range(len(idx)):
        im_mat.append(rref_mat[i,:].T)
    return im_mat


def get_taw(J, F):
    ## τ = 𝐽.T(𝑞)F 
    taw = J.T*F
    return simplify(factor(taw))




# def get_range(mat):
#     vectors = [mat[:,i] for i in range(shape(mat)[1])]
#     range_mat = Matrix.orthogonalize(*vectors, normalize=True)
#     for i, v in enumerate(range_mat):
#         range_mat[i] = simplify(factor(v))
#     return range_mat
    

# t = symbols('t', real=True)

# # phi = Function('phi')(t)
# # theta = Function('theta')(t)
# # psi = Function('psi')(t)

# # px = Function('px')(t)
# # py = Function('py')(t)
# # pz = Function('pz')(t)
# # task = Matrix([px, py, pz, phi, theta, psi])

# q1, q2, q3, q4 = Function('q1')(t), Function('q2')(t), Function('q3')(t), Function('q4')(t)

# r = Matrix([
#     q2*cos(q1) + q4*cos(q1+q3),
#     q2*sin(q1) + q4*sin(q1+q3),
#     q1+q3
# ])

# j_r = get_analytic_jacobian(r, [q1, q2, q3, q4])
# det_j_r = get_det_j(j_r)
# singls = get_singularities_of_j(j_r, [q1, q2, q3, q4])

# j_r_singularity = simplify(factor(j_r.subs([(q2, 0), (q3, 0)])))

# null_jrs = get_nullspace(j_r_singularity.T)

# print(latex(null_jrs))

# print(det_j_r)
# print(simplify(factor(det_j_r.subs(q2, singls[1][q2]))))
# a = j_r*get_null_space(j_r)[0]
# a = j_r*ppp
# print(latex(simplify(factor(a))))
# print(latex(get_nullspace(j_r_singularity)))

# j_task = get_analytic_jacobian(task, )


# w = find_omega_from_RPY_angles(phi, theta, psi, t, 'XYZ')
# rx, ry, rz = symbols('rx, ry, rz', real=True)
# axis = Matrix([rx, ry, rz])

# rotmat, _ = ER.angle_axis_to_rotmat(axis, phi)


# d_rotmat = find_d_rotmat_with_omega(rotmat, w)
# print(latex(d_rotmat))