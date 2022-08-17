#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 14:18:40 2022

@author: akimlavrinenko
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

A = 0.1
epsilon = 0.25
omega = 2*np.pi/10

def a(t):
    return epsilon*np.sin(omega*t)

def b(t):
    return 1-2*epsilon*np.sin(omega*t)

def f(t,x):
    return a(t)*x**2+b(t)*x


def dpsidy(t,x,y):
    return -np.pi*A*np.sin(np.pi*f(t,x))*np.cos(np.pi*y)

def dpsidx(t,x,y):
    return np.pi*A*np.cos(np.pi*f(t,x))*np.sin(np.pi*y)*(2*a(t)*x+b(t))

n = 0.01
t0 = 0
t = -15
deltaT = -0.05


xi = np.arange(0,2+n,n)
yi = np.arange(0,1+n,n)
xi,yi = np.meshgrid(xi,yi,sparse=False)

u = np.zeros(np.shape(xi))
v = np.zeros(np.shape(xi))

xp = np.zeros(np.shape(xi))
yp = np.zeros(np.shape(xi))

#forward
time = np.linspace(t0,t,int((t-t0)/deltaT))
#backward
# time = np.linspace(t0,t,int((t0-t)/deltaT))

x = xi
y = yi

c = '#FFD700'

for t in time:
    for i in range(len(xi[:,0])):
        for j in range(len(yi[0,:])):
            u[i,j] = dpsidy(t,xi[i,j], yi[i,j])
            v[i,j] = dpsidx(t,xi[i,j], yi[i,j])
   
    
    ut = griddata((xi.ravel(), yi.ravel()),(u.ravel()),(x,y),method='linear')
    vt = griddata((xi.ravel(), yi.ravel()),(v.ravel()),(x,y),method='linear')
    xp = np.full(np.shape(xi), deltaT)*ut + x
    yp = np.full(np.shape(xi), deltaT)*vt + y
    print(np.round(t,3))
    x = xp
    y = yp

    

FTLEArr = np.zeros((len(xi[:,0]),len(yi[0,:])))

xy = np.asarray((x,y))
xyi = np.asarray((xi,yi))
DF = np.zeros((2,2))

for i in range(1,len(yi[0,:])-1):
    for j in range(1,len(xi[:,0])-1):
        
        xdeltaT2i = (xy[0,j,i+1] - xy[0,j,i-1])
        xdeltaT2j = (xy[0,j+1,i] - xy[0,j-1,i])
        xdeltaT1i = (xyi[0,j,i+1] - xyi[0,j,i-1])
        
        ydeltaT2i = (xy[1,j,i+1] - xy[1,j,i-1])
        ydeltaT2j = (xy[1,j+1,i] - xy[1,j-1,i])
        ydeltaT1j = (xyi[1,j+1,i] - xyi[1,j-1,i])
        
        s1 = (xdeltaT1i * ydeltaT1j)/2
        s2 = (xdeltaT2i * ydeltaT2j)/2
        
        DF[0,0] = xdeltaT2i/xdeltaT1i
        DF[0,1] = xdeltaT2j/ydeltaT1j
        DF[1,0] = ydeltaT2i/xdeltaT1i
        DF[1,1] = ydeltaT2j/ydeltaT1j
        
        delta = np.dot(DF.T,DF)
        lambdaMax = max(np.linalg.svd(delta)[1])
        FTLE = 1/abs(t - t0) * np.log((lambdaMax**0.5))
        FTLEArr[j,i] = FTLE

FTLEArr[np.isnan(FTLEArr)] = 0
ridge = FTLEArr[FTLEArr > -0.2] = 0

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
plt.contourf(xi,yi,FTLEArr)
plt.colorbar()
plt.xlabel('xi',fontsize=16)
plt.ylabel('yi',fontsize=16)
# plt.xlim(-0.5, 0.5)
# plt.ylim(-0.5, 0.5)
# plt.clim(0,1)
plt.grid()
# plt.savefig('FTLE_doubleGear.png',dpi=150)
plt.show()
