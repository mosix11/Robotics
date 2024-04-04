import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import DKModels as DiK
import DifferetialKinematiks as DK
import InverseDifferentialKinematics as IDK
import Utils as Utils
import DynamicsLE as DLE
from sympy import *

def print_latex(expr):
    print(latex(expr), '\n \n \n')

def simp(exp):
    return simplify(factor(exp))

def iterative_scheme_regulation(dynamics, Kp, Kd, q, t, q_desired, q_initial):
    
    n, _ = shape(q_desired)
    dq = q.diff(t)
    ddq = dq.diff(t)
    
    q_i = [q_initial]
    zero_vector = Matrix([0]*n)
    u_0 = zero_vector
    u_i = [u_0]
    error = [q_desired - q_initial]
    for i in range(1, 10):
        control_effort = Kp*(q_desired - q) + u_i[i-1]
        dynamics_equilibrium_i = dynamics
        for idx, elem in enumerate(dynamics_equilibrium_i):
            dynamics_equilibrium_i[idx] = dynamics_equilibrium_i[idx].subs({ddq[idx]: 0, dq[idx]: 0})
        
        closed_loop_equilibrium_i = dynamics_equilibrium_i - control_effort
        q_i_next = []
        for idx, elem in enumerate(closed_loop_equilibrium_i):
            x = symbols('x', real=True)
            q_i_next.append(nsolve(closed_loop_equilibrium_i[idx].subs(q[idx], x), x, 0))
            
        q_i_next = Matrix(q_i_next)
        q_i.append(q_i_next)
        u_i_next = u_i[i-1] + Kp*(q_desired - q_i_next)
        u_i.append(u_i_next)
        error.append(q_desired - q_i[i])

    for idx, err in enumerate(error):
        print(err.evalf(), u_i[idx].evalf())
        
        