import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import *

def degree_to_radian(degree):
    rad = degree*(pi/180)
    return rad

def radian_to_degree(radian):
    degree = rad/(pi/180)
    return degree

def get_norm(vec):
    return simplify(factor(vec.norm()))
