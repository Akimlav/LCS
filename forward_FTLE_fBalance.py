#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 14:18:40 2022

@author: akimlavrinenko
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from NEK5000_func_lib import particleCoordsAndVel
from os import walk
import time
startTime = time.time()


path = '/Users/akimlavrinenko/Documents/coding/data/room_data/roomBackUp000002/'

files = next(walk(path), (None, None, []))[2]  # [] if no file
fileName = '.3D'
filesList = [word for word in files if fileName in word]
filesList.sort()

n = 0.01
ps = 0
fileStep = 1

filesList = filesList[::fileStep]

t0, a0 = particleCoordsAndVel(path, filesList[0])

xi = yi = np.arange(-0.5,0.5 + n,n)
xi,yi = np.meshgrid(xi,yi)

x = xi
y = yi

for file in range(0,len(filesList)):

    t1, a1 = particleCoordsAndVel (path, filesList[file])
    deltaT = np.round((t1 - t0), 3)
    a01 = np.asarray(a1[ps])
    ac1 = a01[:,0,:]
    av1 = a01[:,1,:]
    
    xv1 = griddata((ac1[:,0],ac1[:,1]),(av1[:,0]),(x, y),method='nearest')
    yv1 = griddata((ac1[:,0],ac1[:,1]),(av1[:,1]),(x, y),method='nearest')
    
    xp2int = np.full(np.shape(xi), deltaT)*xv1 + x
    yp2int = np.full(np.shape(yi), deltaT)*yv1 + y

    x = xp2int
    y = yp2int
    t = t1
    print(t)
    
FTLEArr = np.zeros((len(xi[:,0]),len(yi[0,:])))
xy = np.asarray((x,y))
xyi = np.asarray((xi,yi))
DF = np.zeros((2,2))

print('trajectory integrated')
print('Elapsed time:', np.round((time.time() - startTime),3))

for i in range(1,len(yi[0,:])-1):
    for j in range(1,len(xi[:,0])-1):
        
        xdeltaT2i = (xy[0,j,i+1] - xy[0,j,i-1])
        xdeltaT2j = (xy[0,j+1,i] - xy[0,j-1,i])
        xdeltaT1i = (xyi[0,j,i+1] - xyi[0,j,i-1])
        
        ydeltaT2i = (xy[1,j,i+1] - xy[1,j,i-1])
        ydeltaT2j = (xy[1,j+1,i] - xy[1,j-1,i])
        ydeltaT1j = (xyi[1,j+1,i] - xyi[1,j-1,i])
                
        DF[0,0] = xdeltaT2i/xdeltaT1i
        DF[0,1] = xdeltaT2j/ydeltaT1j
        DF[1,0] = ydeltaT2i/xdeltaT1i
        DF[1,1] = ydeltaT2j/ydeltaT1j
        
        delta = np.dot(DF,DF.T)
        lambdaMax = max(np.linalg.svd(delta)[1])
        FTLE = 1/(t - t0) * np.log((lambdaMax**0.5))
        FTLEArr[j,i] = FTLE

print('plotting FTLE')

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
plt.contourf(xi,yi,FTLEArr)
plt.colorbar()
# plt.xlabel('x',fontsize=14)
# plt.ylabel('y',fontsize=14)
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
# plt.clim(0,1)
plt.grid()
# plt.savefig('FTLE_doubleGear.png',dpi=150)
plt.show()

print('Total elapsed time:', np.round((time.time() - startTime),3))
