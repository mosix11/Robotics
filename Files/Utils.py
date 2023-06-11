import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import *


def print_latex(x):
    print(latex(x), '\n \n \n')

def simp(exp):
    return simplify(factor(exp))

def degree_to_radian(degree):
    rad = degree*(pi/180)
    return rad

def radian_to_degree(radian):
    degree = rad/(pi/180)
    return degree

def get_norm(vec):
    return simplify(factor(vec.norm()))



def vector_to_S(v):
    S = Matrix([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])
    return S
    
def cross_product(r1, r2):
    S_r1 = vector_to_S(r1)
    prod = S_r1*r2
    return prod
    
def filter_matrix_of_small_values(J):
    (m, n) = shape(J)
    mat = []
    for i in range(m):
        row = []
        for j in range(n):
            elem = 0
            if J[i,j] > 0.0001 or J[i,j] < -0.0001:
                elem=float(N(J[i,j]))
                elem=round(elem, 4)
            row.append(elem)
        mat.append(row)
    
    mat = Matrix(mat)
    return mat

def singular_value_decomposition(A):
    AH = A.H
    m, n = A.shape
    if m >= n:
        V, S = (AH * A).diagonalize()
        
        S_diag = [S[i, i] for i in range(S.rows)]
        ranked = []
        for i, x in enumerate(S.diagonal()):
            if not x.is_zero:
                ranked.append(i)
        
        V = V[:, ranked]
        S = Matrix.diag([sqrt(S[i, i]) for i in range(S.rows) if i in ranked])
        
        V, _ = V.QRdecomposition()
        U = A * V * S.inv()
    else:
        U, S = (A * AH).diagonalize()

        S_diag = [S[i, i] for i in range(S.rows)]
        ranked = []
        for i, x in enumerate(S.diagonal()):
            if not x.is_zero:
                ranked.append(i)
                
        U = U[:, ranked]
        S = Matrix.diag([sqrt(S[i, i]) for i in range(S.rows) if i in ranked])
        
        U, _ = U.QRdecomposition()
        V = AH * U * S.inv()

    return U, S, V


def get_matrix_norm(A):
    return simp(A.norm(2)) 


def laplace_simp(eq, qt, qs, t, s):
    lq = []
    for q_i in qt:
        lq.append(laplace_transform(q_i, t, s))
    
    dic = {}
    for i, lq_i in enumerate(lq):
        dic[lq_i] = qs[i]
    
    for q_i in qt:
        dic[q_i.subs(t, 0)] = 0
        
    for q_i in qt:
        dic[q_i.diff(t).subs(t, 0)] = 0
    
    eq = eq.subs(dic)
    
    return eq