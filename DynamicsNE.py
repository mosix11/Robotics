import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import DKModels as DiK
import DifferetialKinematiks as DK
import Utils as Utils
from sympy import *


def print_latex(x):
    print(latex(x), '\n \n \n')

def simp(exp):
    return simplify(factor(exp))

