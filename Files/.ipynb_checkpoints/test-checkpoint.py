import numpy as np
import matplotlib.pyplot as pyplot
import sympy as smp
import swift
import roboticstoolbox as rtb
import spatialmath as sm
import spatialmath.base as smb
import modules.DKModels 
import modules.DifferetialKinematiks as DK

def print_latex(expr):
    print(smp.latex(expr), '\n \n \n')

q1, q2, q3, q4 = smp.symbols('q_1, q_2, q_3, q_4', real=True)
a4 = smp.symbols('a_4', real=True, positive=True)

robot = rtb.DHRobot(
    [rtb.RevoluteDH(alpha=smp.pi/2),
     rtb.RevoluteDH(alpha=smp.pi/2),
     rtb.PrismaticDH(alpha=-smp.pi/2),
     rtb.RevoluteDH(a=a4)],
     name='Random Robot'
)


panda = rtb.models.DH.Panda()




