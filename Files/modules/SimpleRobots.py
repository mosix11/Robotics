import sys
import numpy as np
import matplotlib.pyplot as plt
from modules import DKModels as DiK
from modules import ElemetaryRotations as ER
from modules import MinimalOrientation as MO
from modules import DifferetialKinematiks as DK
from modules import InverseDifferentialKinematics as IDK
from modules import TrajectoryPlaning as TP
from modules import Utils as Utils
from sympy import *

def print_latex(x):
    print(latex(x), '\n \n \n')

def simp(exp):
    return simplify(factor(exp))
#*******************************************************************************************************************************
#*******************************************************************************************************************************

## RR Planar

t = symbols('t', real=True, nonnegative=True)
q1, q2 = Function('q1', real=True, nonnegative=True)(t), Function('q2', real=True, nonnegative=True)(t)
l1, l2 = symbols('l1, l2', real=True, nonnegative=True)

# q1 = 0
# q2 = pi
# l1 = l2 = 1

## Direct Kinematic
px = l1*cos(q1) + l2*cos(q1+q2)
py = l1*sin(q1) + l2*sin(q1+q2)
phi = q1 + q2
r = Matrix([px, py, phi])

## Inverse Kinematik
# if q2 = 0 only 1 solution, if q2 = pi infinite solutions
c2 = (px**2 + py**2 - (l1**2 + l2**2)) / (2*l1*l2)
s2_p = sqrt(1-c2**2)
s2_n = -s2_p
q2_p = simp(atan2(s2_p, c2))
q2_n = simp(atan2(s2_n, c2))
q1_p = simp(atan2(py, px) - atan2(l2*s2_p, l1+l2*cos(q2_p)))
q1_n = simp(atan2(py, px) - atan2(l2*s2_n, l1+l2*cos(q2_n)))
# Elbow up and down
solution1 = (q1_p, q2_p)
solution2 = (q1_n, q2_n)

# print_latex(solution1)
# print_latex(solution2)

#*******************************************************************************************************************************
#*******************************************************************************************************************************

## PR Planar (for RP simpy change q1 and q2 in equations of direct kinematic)

t = symbols('t', real=True, nonnegative=True)
q1, q2 = Function('q1', real=True, nonnegative=True)(t), Function('q2', real=True, nonnegative=True)(t)
L = symbols('L', real=True, nonnegative=True)

px = q1 + L*cos(q2)
py = L*sin(q2)
phi = q2
r = Matrix([px, py])

## Inverse Kinematik
q1_p = px + sqrt(L**2 - py**2)
q1_n = px - sqrt(L**2 - py**2)
q2_p = simp(atan2(py, px - q1_p))
q2_n = simp(atan2(py, px - q1_n))

solution1 = (q1_p, q2_p)
solution2 = (q1_n, q2_n)

#*******************************************************************************************************************************
#*******************************************************************************************************************************

## PP Planar
t = symbols('t', real=True, nonnegative=True)
q1, q2 = Function('q1', real=True, nonnegative=True)(t), Function('q2', real=True, nonnegative=True)(t)

# Direct Kinematics
px = q2
py = q1
r = Matrix([px, py])

# Inverse Kinematik
q1 = px 
q2 = py
solution = (q1, q2)

#*******************************************************************************************************************************
#*******************************************************************************************************************************

### RRR Planar

t = symbols('t', real=True, nonnegative=True)
q1, q2, q3 = Function('q1', real=True, nonnegative=True)(t), Function('q2', real=True, nonnegative=True)(t), Function('q3', real=True, nonnegative=True)(t)
l1 , l2, l3= symbols('l1, l2, l3', real=True, nonnegative=True)

# Dh_table = [
#     [0, l1, q1, 0],
#     [0, l2, q2, 0],
#     [0, l3, q3, 0],
# ]
# dk = DireckKinematic(Dh_table)
# HTM, rot, pos = dk.compute_DK(0, 3)

px = l1*cos(q1) + l2*cos(q1+q2) + l3*cos(q1+q2+q3)
py = l1*sin(q1) + l2*sin(q1+q2) + l3*sin(q1+q2+q3)
# rot = Matrix([
#     [cos(q1+q2+q3), -sin(q1+q2+q3)],
#     [sin(q1+q2+q3), cos(q1+q2+q3)]
# ])
phi = q1+q2+q3
r = Matrix([px, py, phi])

# Inverse Kinematic
# In general if the task r is only the position it will have redundancy (m<n) so it will have oo**n-m solutions
# for Inverse Kinematics. But if we impose the orientation of the E-E as well to the task it will be n=m.
# Then we devide the Inverese Kinematic task to two part. One for the position and one for the orientation




#*******************************************************************************************************************************
#*******************************************************************************************************************************

### PRR Planar 
### Exam Robotics1_21.09.10.pdf
t = symbols('t', real=True, nonnegative=True)
q1, q2, q3 = Function('q1', real=True, nonnegative=True)(t), Function('q2', real=True, nonnegative=True)(t), Function('q3', real=True, nonnegative=True)(t)
L = symbols('L', real=True, nonnegative=True)
# if different link length for two revolut joins use l1 and l2 instead of L
# L = 1
# q1 = 1
# q2 = pi/4
# q3 = pi/2

# Direct Kinematics
px = q1 + L*cos(q2) + L*cos(q2+q3)
py = L*sin(q2) + L*sin(q2+q3)
phi = q2 + q3
r = Matrix([px, py, phi])
px = 2*L
py = 0
phi = symbols('alpha', real=True)

# Inverse Kinematics
q1_p = simp(px - L*cos(phi) + sqrt(L**2*cos(phi)**2 + 2*L*sin(phi)*py - py**2))
q1_n = simp(px - L*cos(phi) - sqrt(L**2*cos(phi)**2 + 2*L*sin(phi)*py - py**2))
q2_p = simp(atan2((py/L)-sin(phi), ((px-q1_p)/L)-cos(phi)))
q2_n = simp(atan2((py/L)-sin(phi), ((px-q1_n)/L)-cos(phi)))
q3_p = phi - q2_p
q3_n = phi - q2_n

solution1 = (q1_p, q2_p, q3_p)
solution2 = (q1_n, q2_n, q3_n)

# print_latex(solution1)
# print_latex(solution2)


#*******************************************************************************************************************************
#*******************************************************************************************************************************

### RPR Planar 
### Exam Robotics1_21.06.11.pdf 
t = symbols('t', real=True, nonnegative=True)
q1, q2, q3 = Function('q1', real=True, nonnegative=True)(t), Function('q2', real=True, nonnegative=True)(t), Function('q3', real=True, nonnegative=True)(t)
L = symbols('L', real=True, nonnegative=True)

L = 0.6
px = q2*cos(q1) + L*cos(q1+q3)
py = q2*sin(q1) + L*sin(q1+q3)
phi = q1+q3
px = 2
py = 0.4
phi = -pi/2
r = Matrix([px, py, phi])
pos = Matrix([px, py])
# J_pos = DK.get_analytic_jacobian(pos, [q1, q2, q3])
# v = Matrix([0, -2.5])
# j_pos_inv = IDK.get_inverse(J_pos)
# q_dot = j_pos_inv*v

# Inverse Kinematics
q2 = simp(sqrt((px-L*cos(phi))**2 + (py-L*sin(phi))**2))
q1 = simp(atan2(py-L*sin(phi), px-L*cos(phi)))
q3 = phi - q1
solution1 = (q1, q2, q3)

print_latex(solution1)
sys.exit()

#*******************************************************************************************************************************
#*******************************************************************************************************************************

## RRP Polar

t = symbols('t', real=True, nonnegative=True)
q1, q2, q3 = Function('q1', real=True, nonnegative=True)(t), Function('q2', real=True, nonnegative=True)(t), Function('q3', real=True, nonnegative=True)(t)
d1 = symbols('d1', real=True, nonnegative=True)

q1 = pi/2
q2 = 0
q3 = 1
## Direct Kinematic
px = q3*cos(q2)*cos(q1)
py = q3*cos(q2)*sin(q1)
pz = d1 + q3*sin(q2)
r = Matrix([px, py, pz])

## Inverse Kinematic
# we consider q3 >= 0
# if q3 = 0 then q1 and q2 are undefined
q3 = simp(sqrt(px**2 + py**2 + (pz-d1)**2))
q2_p = simp(atan2((pz-d1)/q3, sqrt(px**2+py**2)/q3))
q2_n = simp(atan2((pz-d1)/q3, -sqrt(px**2+py**2)/q3))
#if px**2 + py**2 = 0 the q1 is undefined
q1_p = simp(atan2(py/cos(q2_p), px/cos(q2_p)))
q1_n = simp(atan2(py/cos(q2_n), px/cos(q2_n)))

solution1 = (q1_p, q2_p, q3)
solution2 = (q1_n, q2_n, q3)

# print_latex(solution1)
# print_latex(solution2)

#*******************************************************************************************************************************
#*******************************************************************************************************************************

## RRR Elbow Type

t = symbols('t', real=True, nonnegative=True)
q1, q2, q3 = Function('q1', real=True, nonnegative=True)(t), Function('q2', real=True, nonnegative=True)(t), Function('q3', real=True, nonnegative=True)(t)
d1 , l2, l3= symbols('d1, l2, l3', real=True, nonnegative=True)

d1 = l2 = l3 = 1
q1 = 0
q2 = pi/3
q3 = pi/6

## Direct Kinematic
px = cos(q1)*(l2*cos(q2) + l3*cos(q2+q3))
py = sin(q1)*(l2*cos(q2) + l3*cos(q2+q3))
pz = d1 + l2*sin(q2) + l3*sin(q2+q3)
r = Matrix([px, py, pz])

## Inverse Kinematic
c3 = (px**2 + py**2 + (pz-d1)**2 - l2**2 - l3**2) / (2*l2*l3)
s3_p = sqrt(1-c3**2)
s3_n = -s3_p
q3_p = simp(atan2(s3_p, c3))
q3_n = simp(atan2(s3_n, c3))

# if px**2 + py**2 = 0 then q1 is undefined and infinite solutions exist
c1_p = px/(sqrt(px**2+py**2))
c1_n = -c1_p
s1_p = py/(sqrt(px**2+py**2))
s1_n = -s1_p
q1_p = simp(atan2(py, px))
q1_n = simp(atan2(-py, -px))

# if px**2 + py**2 + (pz-d1)**2 = 0 then q2 is undefined
## s1_p & s3_p --> pp
c2, s2 = symbols('c2, s2', real=True)
A = Matrix([
    [l2+l3*c3, -l3*s3_p],
    [l3*s3_p, l2+l3*c3]
])
b = Matrix(
    [c1_p*px+s1_p*py, pz-d1]
)
system = A,b
res = linsolve(system, [c2, s2])
temp = []
for tup in res:
    temp.append(tup)
res = temp[0]
c2, s2 = res[0], res[1]
q2_pp = simp(atan2(s2, c2))
solution_p_p_pp = (q1_p, q2_pp, q3_p)

## s1_n & s3_p --> np
c2, s2 = symbols('c2, s2', real=True)
A = Matrix([
    [l2+l3*c3, -l3*s3_p],
    [l3*s3_p, l2+l3*c3]
])
b = Matrix(
    [c1_n*px+s1_n*py, pz-d1]
)
system = A,b
res = linsolve(system, [c2, s2])
temp = []
for tup in res:
    temp.append(tup)
res = temp[0]
c2, s2 = res[0], res[1]
q2_np = simp(atan2(s2, c2))
solution_n_p_np = (q1_n, q2_np, q3_p)
# print_latex(solution_n_p_np)

## s1_n & s3_n --> nn
c2, s2 = symbols('c2, s2', real=True)
A = Matrix([
    [l2+l3*c3, -l3*s3_n],
    [l3*s3_n, l2+l3*c3]
])
b = Matrix(
    [c1_n*px+s1_n*py, pz-d1]
)
system = A,b
res = linsolve(system, [c2, s2])
temp = []
for tup in res:
    temp.append(tup)
res = temp[0]
c2, s2 = res[0], res[1]
q2_nn = simp(atan2(s2, c2))
solution_n_n_nn = (q1_n, q2_nn, q3_n)
# print_latex(solution_n_n_nn)

## s1_p & s3_n --> pn
c2, s2 = symbols('c2, s2', real=True)
A = Matrix([
    [l2+l3*c3, -l3*s3_n],
    [l3*s3_n, l2+l3*c3]
])
b = Matrix(
    [c1_p*px+s1_p*py, pz-d1]
)
system = A,b
res = linsolve(system, [c2, s2])
temp = []
for tup in res:
    temp.append(tup)
res = temp[0]
c2, s2 = res[0], res[1]
q2_pn = simp(atan2(s2, c2))
solution_p_m_pn = (q1_n, q2_pn, q3_n)
# print_latex(solution_p_m_pn)


#*******************************************************************************************************************************
#*******************************************************************************************************************************

