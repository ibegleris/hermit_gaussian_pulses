# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 13:02:53 2015
Finds the widths of gaussian pulses for the first two modes within an optical fibre
that has a parabolic index profile. 
@author: john
"""
from __future__ import division
import numpy as np
from scipy.special import hermite as h
import matplotlib.pylab as plt
from scipy.integrate import dblquad
from scipy.optimize import fsolve

def psi(y,x,omega,l,m,n = 2):
    return np.abs(h(l)(2**0.5*x/omega)*h(m)(2**0.5*y/omega)*np.exp(-(x**2+y**2)/omega**2))**n

def effective(omega,xmin,xmax,ymin,ymax,l,m):
    omega1,omega2 = omega    
    n =2
    a = dblquad(psi,xmin,xmax,lambda x : ymin,lambda x: ymax,args = (omega1,0,0,n))[0]
    n = 4
    eff1 = a**2 / dblquad(psi,x.min(),x.max(),lambda x : ymin,lambda x: ymax,args = (omega1,0,0,n))[0]

    n =2
    a = dblquad(psi,xmin,xmax,lambda x : ymin,lambda x: ymax,args = (omega2,0,1,n))[0]
    n = 4
    eff2 = a**2 / dblquad(psi,x.min(),x.max(),lambda x : ymin,lambda x: ymax,args = (omega2,0,1,n))[0]
    
    return (eff1 - 1.61e-10,eff2 - 1.70e-10)




def plotting(xmin,xmax,ymin,ymax,omega,l,m):
    x = np.linspace(xmin,xmax,256)
    y = np.linspace(ymin,ymax,256)
    z = np.zeros([len(x),len(y)])
    for i,xx in enumerate(x):
        for j, yy in enumerate(y):
            z[i,j] = psi(xx,yy,omega,l,m)
    
    plt.pcolormesh(x, y, z/np.max(z))
    plt.title('pcolor')
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.colorbar()
    plt.show()
    return
    
    
def main():
    f_diam = 124.9e-6
    f_r = f_diam*0.5
    xmin = -f_r/3
    xmax =f_r /3
    ymin = -f_r/3
    ymax = f_r/3
    omega_0 =  7.15876961e-06
    l=0
    m=0
    omega = fsolve(effective,(omega_0,omega_0),args = (xmin,xmax,ymin,ymax,l,m))
    plotting(xmin,xmax,ymin,ymax,omega[0],0,0)
    plotting(xmin,xmax,ymin,ymax,omega[1],0,1)
if __name__ == "__main__":
    main()