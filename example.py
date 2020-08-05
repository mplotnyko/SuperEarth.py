import matplotlib.pyplot as plt
import numpy as np
from SuperEarth import *

data = np.loadtxt('data.csv',skiprows=1,delimiter=',').T

#1 Plot M-R with uncompressed density
fig=plt.figure(dpi=150)
cmf = f_cmf(data[0],data[1])
rho_0 = f_rho0(cmf)
cb = plt.scatter(data[0],data[1],c=rho_0/1e3)
plt.colorbar(cb,label=r'$\rho_0$ [g/cm$^3$]')

plt.xscale('log')
plt.yscale('log')
plt.xlim(1,15)
plt.ylim(1,2.5)
plt.ylabel(r'R/R$_\oplus$')
plt.xlabel(r'M/M$_\oplus$')
plt.savefig('MR_rho0.jpg',fig=fig,bbox_inches='tight')

#2 Plot FeSi of the star in M-R
fig=plt.figure(dpi=150)
FeSi = np.random.normal(2,1,500)
cmf = star_to_planet(FeSi)
X,Y,Z = plot_cont(cmf)
cb = plt.contourf(X,Y,Z,100,cmap='Reds')
plt.colorbar(cb,label='Normalized Probability')
plt.scatter(data[0],data[1])

plt.xscale('log')
plt.yscale('log')
plt.xlim(1,15)
plt.ylim(1,2.5)
plt.ylabel(r'R/R$_\oplus$')
plt.xlabel(r'M/M$_\oplus$')
#plot the estimates for sample compositions
c = [0,0.32,0.67,0.999]
x = np.linspace(1,15,30)
for i,item in enumerate(c):
    R = np.zeros(len(x))
    for j,item in enumerate(x):
        R[j] = guess_R(item,0,cmf=c[i])
    plt.plot(x,R,c='k',zorder=10)
plt.savefig('FeSi_star.jpg',fig=fig,bbox_inches='tight')

