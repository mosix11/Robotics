import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import *


q = symbols('q')
elm_rot_Z = Matrix([
    [cos(q), -sin(q), 0],
    [sin(q), cos(q), 0],
    [0, 0, 1]
])

elm_rot_X = Matrix([
    [1, 0, 0],
    [0, cos(q), -sin(q)],
    [0, sin(q), cos(q)]
])

elm_rot_Y = Matrix([
    [cos(q), 0, sin(q)],
    [0, 1, 0],
    [-sin(q), 0, cos(q)]
])


def get_variable_elem_rot(angle_func, axis):
    if axis == 'Z':
        rotmat = elm_rot_Z.subs([(q, angle_func)])
        return rotmat
    elif axis == 'X':
        rotmat = elm_rot_X.subs([(q, angle_func)])
        return rotmat
    elif axis == 'Y':
        rotmat = elm_rot_Y.subs([(q, angle_func)])
        return rotmat
    else:
        return 'axis undefined'
    


def do_elem_rot(item, axis, theta, post_multi=False):
    if axis == 'Z':
        rot = elm_rot_Z.subs([(q, theta)])
        return item*rot if post_multi else rot*item
    elif axis == 'X':
        rot = elm_rot_X.subs([(q, theta)])
        return item*rot if post_multi else rot*item
    elif axis == 'Y':
        rot = elm_rot_Y.subs([(q, theta)])
        return item*rot if post_multi else rot*item
    else:
        return 'axis undefined'


def angle_axis_to_rotmat(axis, angle):
    q = angle
    rx, ry, rz = axis[0], axis[1], axis[2]
    r = Matrix([rx, ry, rz])
    S_r = Matrix([
        [0, -rz, ry],
        [rz, 0, -rx],
        [-ry, rx, 0]
    ])
    rotmat = r*r.T + (eye(3) - r*r.T)*cos(q) + S_r*sin(q)
    rotmat = trigsimp(factor(rotmat))
    
    rotmat = rotmat.subs([(rx, axis[0]), (ry, axis[1]), (rz, axis[2]), (q, angle)])
    # trace = trigsimp(factor(rot_angle_axis.trace()))
    trace_raa = 1 + 2*cos(q)
    
    return rotmat, trace_raa

def rotmat_to_angle_axis(rot_mat):
    R11, R12, R13 = rot_mat[0, :]
    R21, R22, R23 = rot_mat[1, :]
    R31, R32, R33 = rot_mat[2, :]
    
    sinx_1 = sqrt((R12 - R21)**2 + (R13 - R31)**2 + (R23 - R32)**2) / 2
    sinx_2 = -sinx_1
    cosx = (R11 + R22 + R33 - 1) / 2
    
    angle1 = atan2(sinx_1, cosx)
    angle2 = atan2(sinx_2, cosx)
    
    if sinx_1 == 0 or sinx_2 == 0:
        if angle1 == 0:
            return (angle1, 'undefined'), (angle2, 'undefined')
        elif angle1 == pi:
            
            rx, ry, rz = symbols('rx, ry, rz')
            eq1 = rx*ry - R12/2
            eq2 = rx*rz - R13/2
            eq3 = ry*rz - R23/2
            eq4 = rx**2 - (R11+1)/2
            eq5 = ry**2 - (R22+1)/2
            eq6 = rz**2 - (R33+1)/2
            res = nonlinsolve([eq1, eq2, eq3, eq4, eq5, eq6], [rx, ry, rz])
            axes = []
            for axis in res:
                axes.append(simplify(factor(Matrix(axis))))
                
            axis1_norm = simplify(axes[0].norm())
            axis2_norm = simplify(axes[1].norm())
            print('norms', axis1_norm, axis2_norm)
            return (-pi, axes[0]), (pi, axes[1])
            
    else:
        axis1 = (Rational(1, 2*sinx_1))*Matrix([R32-R23, R13-R31, R21-R12])
        axis2 = (Rational(1, 2*sinx_2))*Matrix([R32-R23, R13-R31, R21-R12])
        axis1_norm = simplify(axis1.norm())
        axis2_norm = simplify(axis2.norm())
        print('norms', axis1_norm, axis2_norm)
        return (angle1, axis1), (angle2, axis2)
    





# rot, _ = get_rot_angle_axis([0, 0, 1], pi)
# (angle1, axis1), (angle2, axis2) = find_angle_axis(rot)

# rot_mat = Matrix([
#     [-1, 0, 0],
#     [0, -1/sqrt(2), -1/sqrt(2)],
#     [0, -1/sqrt(2), 1/sqrt(2)]
# ])
# print(latex(find_angle_axis(rot_mat)))




