import numpy as np
import math
from matplotlib import pyplot as plt

plt.rcParams['text.usetex'] = True

with open("vel.txt") as f:
    lines = f.readlines()

N = int(float(lines[0].split(":")[1].strip()))
TrueTemp = float(lines[1].split(":")[1].strip())

TrueCx = float(lines[2].split(":")[1].strip())
TrueCy = float(lines[3].split(":")[1].strip())
TrueCz = float(lines[4].split(":")[1].strip())

#print(f"{N} {TrueTemp} {TrueCx}")

cx, cy, cz = np.loadtxt('vel.txt', delimiter=',', usecols=(0,1,2), unpack=True, skiprows=7)

#print(np.isnan(cz).any())
#print(cx)

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

velBins, edges = np.histogram(cx, bins=N//10, density=True)

beta = (m_ar/2/k/TrueTemp)**0.5
norm = (m_ar/2/k/TrueTemp/3.1415926)**0.5

Cdist = np.linspace(-TrueCx, 3*TrueCx, num=10000)

velVar = cx.var()
#thermDist = beta/(3.141592)**0.5*np.exp(-beta*beta*Cdist*Cdist)
#velDist = norm*np.exp(-m_ar/2/k/TrueTemp*Cdist*Cdist)
velDist = (m_ar/2/3.1415926/k/TrueTemp)**0.5*np.exp(-m_ar/2/k/TrueTemp*(Cdist-TrueCx)**2)

#print(N)

fig, ax = plt.subplots()
plt.plot(Cdist, velDist, label=r'True $C_x$ distribution',color='b')
plt.hist(cx, bins='auto',density=True, label=r'$C_x$ distribution with' + '\n' + r'{} Particles'.format(N), color='#f55742')
plt.grid()
ax.legend(loc='upper right', fontsize=10)
ax.set_xlabel(r'$C_x$')
ax.set_ylabel(r'$f_0(C_x)$')
plt.savefig(f'/Users/abbieb/Documents/School/Grad/Spring2026/DSMC/HW1/Figures/N{N}_V{TrueCx}_T{TrueTemp}.pdf')
#plt.show()
