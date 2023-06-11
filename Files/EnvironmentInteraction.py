from sympy import *
import Utils as Utils
import MatrixInversion as MI

def print_latex(x):
    print(latex(x), '\n \n \n')

def simp(exp):
    return simplify(factor(exp))


def compute_J_of_constraints_A(h, q):
    A = h.jacobian(q)
    return simp(A)


def is_sigular(mat):
    det = mat.det()
    if det == 0:
        return True
    else:
        return False

def compute_lambda_from_dynamic_and_A(A, M, n, u, q, t):
    dq = diff(q, t)
    M_inv = MI.get_inverse(M)
    lamda = MI.get_inverse(A*M_inv*A.T) * (A*M_inv*(n-u) - diff(A, t)* dq)
    return simp(lamda)


## TODO implement if needed
# def compute_constrained_dynamics(A, lamda, M, n, q, t):
#     u = 
#     return 
    

def compute_E_F_from_A_D(A, D):
    sA,c = shape(A)
    sD,c = shape(D)
    T = Matrix.vstack(A, D)
    T_inv = MI.get_inverse(T)
    E = T_inv.extract(range(0,c), range(0, sA))
    F = T_inv.extract(range(0,c), range(sA, sA+sD))
    E = simp(E)
    F = simp(F)
    return E, F


def compute_lambda_from_dynamic_and_A_D_E_F(A, D, E, F, M, n, u, q, t):
    dq = diff(q, t)
    ddq = diff(dq, t)
    v = D*dq
    dv = D*ddq + diff(D, t)*dq
    
    lamda = E.T * (M*F*dv - M*(E*diff(A, t) + F*diff(D, t))*F*v + n - u)
    return simp(lamda)

def compute_reduced_inertia(M, F):
    return simp(F.T * M * F)

## TODO implement if needed
def compute_reduced_dynamics(A, D, E, F, M , n , q, u , t):
    
    return None


