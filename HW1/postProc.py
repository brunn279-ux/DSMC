import numpy as np
from matplotlib import pyplot as plt

cx, cy, cz = np.loadtxt('vel.txt', delimiter=',', usecols=(0,1,2), unpack=True, skiprows=2)

#print(np.isnan(cz).any())
print(cx)

k = 1.380649e-23 #m^2*kg*s^-2*K^-1
M_ar = 39.95/1000 # kg/mol
m_ar = M_ar/6.02e23
R_u = 8.314 # J/mol/K
R_ar = R_u/M_ar

u = np.mean(cx)
v = np.mean(cy)
w = np.mean(cz)

print(f"C = ({u}, {v}, {w})")

T = 0.0

for i in range(cx.size):
    T += 1/3/k*m_ar*(cx[i-1]*cx[i-1]+cy[i-1]*cy[i-1]+cz[i-1]*cz[i-1]-u*u-v*v-w*w)
    #print(f"T = {T/i}, C = ({cx[i-1]}, {cy[i-1]}, {cz[i-1]})")

T = T/cx.size
print(f"T = {T}")
#print(cx)
