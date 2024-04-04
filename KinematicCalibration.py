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


def callibrate(DH_table, experiments, q_vec, nominal_values, t, callibrate_d=False, callibrate_a=False, callibrate_theta=False, callibrate_alpha=False):
    
    dk_model = DiK.DireckKinematic(DH_table)
    DOF = dk_model.get_DOF()
    mat, rot, pos = dk_model.compute_DK(0, DOF)
    
    d_vec = Matrix([col[0] for col in DH_table])
    a_vec = Matrix([col[1] for col in DH_table])
    theta_vec = Matrix([col[2] for col in DH_table])
    alpha_vec = Matrix([col[3] for col in DH_table])

    f = pos[:2] # set based on the question - it could be 2d postion, 2d position and orientation or full cartesian 3d positiona and orientation
    f = Matrix(f) # slicing returns list
    # print(type(f))
    # print(shape(f))
    # print_latex(f)

    f_nom = f.evalf(subs=nominal_values)
    
    phi = None
    if callibrate_alpha:
        phi_alpha = simp(f.jacobian(alpha_vec))
        phi = phi_alpha if phi==None else Matrix.hstack(phi, phi_alpha)
    if callibrate_a:
        phi_a = simp(f.jacobian(a_vec))
        phi = phi_a if phi==None else Matrix.hstack(phi, phi_a)
    if callibrate_d:
        phi_d = simp(f.jacobian(d_vec))
        phi = phi_d if phi==None else Matrix.hstack(phi, phi_d)
    if callibrate_theta:
        phi_theta = simp(f.jacobian(theta_vec))
        phi = phi_theta if phi==None else Matrix.hstack(phi, phi_theta)
    
    

    dr_bar = None
    phi_bar = None
    
    for conf, r in experiments:

        subs = dict(zip(q_vec, conf))
        
        conf = Matrix(conf)
        r = Matrix(r)
        
        f_conf = f_nom.evalf(subs=subs)
        delta_r = r - f_conf
        phi_conf = phi.evalf(subs=subs)
        
        dr_bar = delta_r if dr_bar == None else Matrix.vstack(dr_bar, delta_r)
        phi_bar = phi_conf if phi_bar == None else Matrix.vstack(phi_bar, phi_conf)
        
        
    # print_latex(dr_bar)
    # print_latex(phi_bar)
    
    delta_params = IDK.get_inverse(phi_bar) * dr_bar
    
    print_latex(delta_params)
    
        



t = symbols('t', real=True, nonnegative=True)

q1, q2 = Function('q_1', real=True, nonnegative=True)(t), Function('q_2', real=True, nonnegative=True)(t)
q_vec = Matrix([q1, q2])


# a1, a2 = symbols('a_1, a_2', real=True, nonnegative=True) 
# a_vec = Matrix([a1, a2])
# d1, d2 = symbols('d_1, d_2', real=True)
# d_vec = Matrix([d1, d2])
# alpha1, alpha2 = symbols('alpha_1, alpha_2', real=True)
# alpha_vec = Matrix([alpha1, alpha2])

# DH_table = [
#     [d1, a1, q1, alpha1],
#     [d2, a2, q2, alpha2]
# ]

l1, l2 = symbols('l_1, l_2', real=True, nonnegative=True) # a1, a2

DH_table = [
    [0, l1, q1, 0],
    [0, l2, q2, 0]
]

nominal_values = {l1: 1, l2: 1}

experiments = [
    ([0,0], [2, 0]),
    ([pi/2, 0], [0, 2]),
    ([pi/4, -pi/4], [1.6925, 0.7425]),
    ([0, pi/4], [1.7218, 0.6718])
]

callibrate(DH_table, experiments, q_vec, nominal_values, t, callibrate_a=True)