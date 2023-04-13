import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import *

from modules import ElemetaryRotations as ER
from modules import MinimalOrientation as MO


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
            
def get_F(J, taw):
    ## F = 𝐽**-1.T(𝑞).τ
    J_inv = get_inverse(J)
    F = J_inv.T*taw
    return F