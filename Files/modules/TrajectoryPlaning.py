import numpy as np
import matplotlib.pyplot as plt
from DKModels import DireckKinematic
from modules import ElemetaryRotations as ER
from modules import MinimalOrientation as MO
from modules import DifferetialKinematiks as DK
from modules import InverseDifferentialKinematics as IDK
from sympy import *


def get_velocity(p, s, t):
    vel = diff(p, t)
    return simplify(factor(vel))

def get_tangent_of_v(p, s, t):
    tv = diff(p, s)
    return simplify(factor(tv))

def get_acceleration(p, s, t):
    vel = get_velocity(p, s, t)
    acc = diff(vel, t)
    return simplify(factor(acc))

def get_tangent_of_acc(p, s, t):
    tv = get_tangent_of_v(p, s, t)
    ta = diff(tv, s)
    return simplify(factor(ta))

def get_Fernet_frame(p, s, t):
    # p: path --> p(s), s = parameter --> s(t), t = time
    ts = diff(p, s)/diff(p, s).norm() # unit tangent vector
    ts = simplify(factor(ts))
    ns = diff(ts, s)/diff(ts, s).norm() # unit normal vector
    ns = simplify(factor(ns))
    bs = ts.cross(ns) # unit binormal vector
    bs = simplify(factor(bs))
    return ts, ns, bs
    
    
def get_curvature(p, s, t):
    # p: path --> p(s), s = parameter --> s(t), t = time
    p_first_der = diff(p, s)
    p_second_der = diff(p_first_der, s)
    p_third_der = diff(p_second_der, s)
    ks = (p_first_der.cross(p_second_der).norm()) / (p_first_der.norm())**3
    return simplify(factor(ks))
    
def get_torsion(p, s, t):
    # p: path --> p(s), s = parameter --> s(t), t = time
    p_first_der = diff(p, s)
    p_second_der = diff(p_first_der, s)
    p_third_der = diff(p_second_der, s)
    ts = (p_first_der.T*(p_second_der.cross(p_third_der))) / (p_first_der.cross(p_second_der).norm())**2
    return simplify(factor(ts))



def trapezoidal_get_Ts(v_max, a_max):
    return simplify(factor(v_max/a_max))

def trapezoidal_get_T(v_max, a_max, L):
    T = (L*a_max + v_max**2) / (a_max*v_max)
    return simplify(factor(T))

def trapezoidal_has_coast(v_max, a_max, L):
    x = (v_max**2)/a_max
    return L>x, x


def cubic_profile_joint():
    t = symbols('t', real=True, nonnegative=True)
    T = symbols('T', real=True, positive=True)
    a, b, c, d = symbols('a, b, c, d', real=True)
    q_in, q_fin, v_in, v_fin = symbols('q_in, q_fin, v_in, v_fin', real=True)
    
    # v_in, v_fin = 0, 0, 0, 0
    # T = 1
    # p_in = 0
    # p_fin = 1
    
    q = q_in + (q_fin-q_in)*(a*t**3 + b*t**2 + c*t + d)
    q_vel = diff(q, t)
    
    # Imposing conditions
    q_0 = q.subs(t, 0)
    q_T = q.subs(t, T)
    q_vel_0 = q_vel.subs(t, 0)
    q_vel_T = q_vel.subs(t, T)
    
    # Equations to find a,b,c,d
    eq1 = Eq(q_0, q_in) # First position at t=0
    eq2 = Eq(q_T, q_fin) # Last position at t=T
    eq3 = Eq(q_vel_0, v_in) # First velocity at t=0
    eq4 = Eq(q_vel_T, v_fin) # Last velocity at t=T
    
    res = solve([eq1, eq2, eq3, eq4], [a, b, c, d])
    print(res)
    
    
    
def cubic_profile_cartesian():
    t = symbols('t', real=True, nonnegative=True)
    T = symbols('T', real=True, positive=True)
    a, b, c, d = symbols('a, b, c, d', real=True)
    p_in, p_fin, v_in, v_fin = symbols('p_in, p_fin, v_in, v_fin', real=True)

    # v_in, v_fin = 0, 0, 0, 0
    # T = 1
    # q_in = 0
    # q_fin = 1

    p = p_in + (p_fin-p_in)*(a*t**3 + b*t**2 + c*t + d)
    p_vel = diff(p, t)

    
    
    # Imposing conditions
    p_0 = p.subs(t, 0)
    p_T = p.subs(t, T)
    p_vel_0 = p_vel.subs(t, 0)
    p_vel_T = p_vel.subs(t, T)
    
    # Equations to find a,b,c,d
    eq1 = Eq(p_0, p_in) # First position at t=0
    eq2 = Eq(p_T, p_fin) # Last position at t=T
    eq3 = Eq(p_vel_0, v_in) # First velocity at t=0
    eq4 = Eq(p_vel_T, v_fin) # Last velocity at t=T
    
    res = solve([eq1, eq2, eq3, eq4], [a, b, c, d])
    print(res)


def quintic_profile():
    t = symbols('t', real=True, nonnegative=True)
    T = symbols('T', real=True, positive=True)
    a, b, c, d, e, f = symbols('a, b, c, d, e, f', real=True)
    q_in, q_fin, v_in, v_fin, a_in, a_fin = symbols('q_in, q_fin, v_in, v_fin, a_in, a_fin', real=True)
    
    # v_in, v_fin, a_in, a_fin = 0, 0, 0, 0
    # T = 1
    # q_in = 0
    # q_fin = 1
    
    q = q_in + (q_fin-q_in)*(a*t**5 + b*t**4 + c*t**3 + d*t**2 + e*t + f)
    q_vel = diff(q, t)
    q_acc = diff(q_vel, t)
    
    # Imposing conditions
    q_0 = q.subs(t, 0)
    q_T = q.subs(t, T)
    q_vel_0 = q_vel.subs(t, 0)
    q_vel_T = q_vel.subs(t, T)
    q_acc_0 = q_acc.subs(t, 0)
    q_acc_T = q_acc.subs(t, T)
    
    # Equations to find a,b,c,d,e,f
    eq1 = Eq(q_0, q_in) # First position at t=0
    eq2 = Eq(q_T, q_fin) # Last position at t=T
    eq3 = Eq(q_vel_0, v_in) # First velocity at t=0
    eq4 = Eq(q_vel_T, v_fin) # Last velocity at t=T
    eq5 = Eq(q_acc_0, a_in)
    eq6 = Eq(q_acc_T, a_fin)

    
    res = solve([eq1, eq2, eq3, eq4, eq5, eq6], [a, b, c, d, e, f])
    print(latex(res))



def profile_4_3_4():
    t = symbols('t', real=True, nonnegative=True)
    T = symbols('T', real=True, positive=True)
    
    
def bang_cost_bang_profile():
    t = symbols('t', real=True, nonnegative=True)
    T = symbols('T', real=True, positive=True)
    ts = symbols('ts', real=True, nonnegative=True)
    
    q_in, q_fin = symbols('q_in, q_fin', real=True)
    
    v_max, a_max = symbols('v_max, a_max', real=True)
    
    L = symbols('L', real=True, positive=True)
    delta = Function('delta')(t)
    
    ts = v_max/a_max
    T = (L*a_max + v_max**2) / (a_max*v_max)
    
    a_max = v_max/ts
    v_max = ts*a_max
    
    if t < ts:
        delta = (a_max**2)/2
    elif t >= ts and t <= (T-ts):
        delta = (v_max*t) - ((v_max**2)/(2*a_max))
    else:
        delta = -((a_max*(t-T)**2)/2) + v_max*T - ((v_max**2)/a_max)
    
    q = q_in + (q_fin-q_in)*(delta/L)
    
    
    

    


def spline_profile_2_cubic():
    t = symbols('t', real=True, nonnegative=True)
    T = symbols('T', real=True, positive=True)
    
    a, b, c, d = symbols('a1, b1, c1, d1', real=True)
    q_in, q_fin, v_in, v_fin = symbols('q_in, q_fin, v_in, v_fin', real=True)
    
    # v_in, v_fin = 0, 0, 0, 0
    # T = 1
    # p_in = 0
    # p_fin = 1
    
    q = q_in + (q_fin-q_in)*(a*t**3 + b*t**2 + c*t + d)
    q_vel = diff(q, t)
    
    # Imposing conditions
    q_0 = q.subs(t, 0)
    q_T = q.subs(t, T)
    q_vel_0 = q_vel.subs(t, 0)
    q_vel_T = q_vel.subs(t, T)
    
    # Equations to find a,b,c,d
    eq1 = Eq(q_0, q_in) # First position at t=0
    eq2 = Eq(q_T, q_fin) # Last position at t=T
    eq3 = Eq(q_vel_0, v_in) # First velocity at t=0
    eq4 = Eq(q_vel_T, v_fin) # Last velocity at t=T
    
    res = solve([eq1, eq2, eq3, eq4], [a, b, c, d])

