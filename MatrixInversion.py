import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import *

import ElemetaryRotations as ER
import MinimalOrientation as MO



def get_inverse(J):
    (m, n) = shape(J)
    is_square = m == n
    rank_j = J.rank()
    if is_square:
        is_fullrank = rank_j == m
        if is_fullrank:
            return simplify(J**-1)
        else:
            return simplify(J.pinv())
            
    else:
        if m<n:
            is_fullrank = rank_j == m
            if is_fullrank:
                return simplify(J.T*((J*J.T)**-1))
            else:
                return simplify(J.pinv())
        else:
            is_fullrank = rank_j == n
            if is_fullrank:
                return simplify(((J.T*J)**-1)*J.T)
            else:
                return simplify(J.pinv())


def get_weithed_pinv(J, W):
    (m, n) = shape(J)
    is_square = m == n
    rank_j = J.rank()
    
    if is_square:
        print('alllleeerrrttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt')

        pass
        # is_fullrank = rank_j == m
        # if is_fullrank:
        #     return simplify(J**-1)
        # else:
        #     return simplify(J.pinv())
            
    else:
        if m<n:
            is_fullrank = rank_j == m
            if is_fullrank:
                return simplify((W**-1)*J.T*((J*(W**-1)*J.T)**-1))
            else:
                J_aux = J * W**(-1/2)
                Jw_pinv = W**(-1/2) * J_aux.pinv()
                return simplify(Jw_pinv)
        else:
            is_fullrank = rank_j == n
            if is_fullrank:
                return simplify(((J.T*(W**-1)*J)**-1)*J.T*(W**-1))
            else:
                # not sure
                J_aux = J * W**(-1/2)
                Jw_pinv = W**(-1/2) * J_aux.pinv()
                return simplify(Jw_pinv)



def get_J_DLS(J, mu):
    (m, n) = shape(J)
    is_square = m == n
    rank_j = J.rank()
    
    if is_square:
        is_fullrank = rank_j == m
        if is_fullrank:
            return simplify(J.T*((J*J.T + (mu**2)*eye(m))**-1))
        else:
            # implement with svd in utils and slide 21 of kinematic redundancy
            print('alllleeerrrttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt')
            return simplify(J.T*(get_inverse((J*J.T + (mu**2)*eye(m)))))
            
    else:
        if m<n:
            is_fullrank = rank_j == m
            if is_fullrank:
                return simplify(J.T*((J*J.T + (mu**2)*eye(m))**-1))
            else:
                # implement with svd in utils and slide 21 of kinematic redundancy
                print('alllleeerrrttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt')
                return simplify(J.T*(get_inverse((J*J.T + (mu**2)*eye(m)))))
        else:
            is_fullrank = rank_j == n
            if is_fullrank:
                return simplify(((J.T*J + (mu**2)*eye(n))**-1)*J.T)
            else:
                # implement with svd in utils and slide 21 of kinematic redundancy
                print('alllleeerrrttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt')
                return simplify((get_inverse((J.T*J + (mu**2)*eye(n))))*J.T)

def get_F(J, taw):
    ## F = 𝐽**-1.T(𝑞).τ
    J_inv = get_inverse(J)
    F = J_inv.T*taw
    return F