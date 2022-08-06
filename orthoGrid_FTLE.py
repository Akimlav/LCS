#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 12:31:23 2022

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

n = 0.0025
deltaT = 1e-3
ps = 0
fileStep = 1

filesList = filesList[::fileStep]
fileList = filesList[:2]

t0, a0 = particleCoordsAndVel (path, filesList[-1])
a00 = np.asarray(a0[ps])
ac0 = a00[:,0,:]
av0 = a00[:,1,:]

# target grid to interpolate to
xi = yi = np.arange(-0.5,0.5 + n,n)
xi,yi = np.meshgrid(xi,yi)

# interpolate
zx0 = griddata((ac0[:,0],ac0[:,1]),(av0[:,0]),(xi,yi),method='linear')
zy0 = griddata((ac0[:,0],ac0[:,1]),(av0[:,1]),(xi,yi),method='linear')

xp1int = np.full(np.shape(xi), deltaT)*zx0 + xi
yp1int = np.full(np.shape(yi), deltaT)*zy0 + yi

# fig = plt.figure(figsize=(7,7))
# ax = fig.add_subplot(111)
# # plt.contourf(xi,yi,zy)
# # plt.colorbar()
# plt.plot(xi,yi,'k.')
# # plt.plot(xp2int,yp2int,'r.')
# plt.plot(xp1int,yp1int,'b.')
# plt.xlabel('xi',fontsize=16)
# plt.ylabel('yi',fontsize=16)
# plt.xlim(-0.5, -0.4)
# plt.ylim(-0.5, -0.4)
# plt.grid()
# plt.show()

x = xp1int
y = yp1int
t = t0


for file in reversed(range(0,len(filesList)-1)):
# for file in range(4):
    
    t1, a1 = particleCoordsAndVel (path, filesList[file])
    deltaT = np.round((t1 - t), 2)
    print(deltaT)
    a01 = np.asarray(a1[ps])
    ac1 = a01[:,0,:]
    av1 = a01[:,1,:]
    
    xv1 = griddata((ac1[:,0],ac1[:,1]),(av1[:,0]),(x, y),method='linear')
    yv1 = griddata((ac1[:,0],ac1[:,1]),(av1[:,1]),(x, y),method='linear')
    
    xp2int = np.full(np.shape(xi), deltaT)*xv1 + x
    yp2int = np.full(np.shape(yi), deltaT)*yv1 + y

    x = xp2int
    y = yp2int
    t = t1
    
    
    # plt.contourf(xi,yi,zy)
    # plt.colorbar()

x = np.nan_to_num(x)
y = np.nan_to_num(y)
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)
plt.plot(xi,yi,'k.', markersize=1)
plt.plot(x,y,'r.', markersize=1)
# plt.plot(xp1int,yp1int,'b.')
plt.xlabel('xi',fontsize=16)
plt.ylabel('yi',fontsize=16)
plt.xlim(-0.5, -0.4)
plt.ylim(-0.5, -0.4)
plt.grid()
plt.show()

ftleList = []
FTLEArr = np.zeros((len(x)-2,len(y)-2))
FTLEArr0 = np.zeros((len(x),len(y)))

for i in range(1,len(x)-1):
    for j in range(1,len(y)-1):
        # print(i,j)
        # print(x[i,j])
        xy = np.asarray((x,y))
        xyi = np.asarray((xi,yi))
        
        xdeltaT2 = xy[0,j,i+1] - xy[0,j,i-1]
        xdeltaT1 = xyi[0,j,i+1] - xyi[0,j,i-1]
        
        ydeltaT2 = xy[1,i+1,j] - xy[1,i-1,j]
        ydeltaT1 = xyi[1,i+1,j] - xyi[1,i-1,j]
        
        s1 = (xdeltaT1 * ydeltaT1)/2
        s2 = (xdeltaT2 * ydeltaT2)/2
        
        DF = np.zeros((2,2))
        DF[0,0] = xdeltaT2/xdeltaT1
        DF[0,1] = xdeltaT2/ydeltaT1
        DF[1,0] = ydeltaT2/xdeltaT1
        DF[1,1] = ydeltaT2/ydeltaT1
        
        delta = np.dot(DF,DF.T)
        lambdaMax = max(np.linalg.svd(delta)[1])
        FTLE = 1/(t1 - t0) * np.log((lambdaMax**0.5))
        # FTLEArr[i,j] = FTLE
        FTLEArr0[i,j] = FTLE
        # ftleList.append(FTLE)
    


# resArr = np.zeros((len(np.ravel(xi)),3))
# resArr[:,0] = np.ravel(xi)
# resArr[:,1] = np.ravel(yi)
# resArr[:,2] = np.asarray(ftleList)


# target grid to interpolate to
# xi = yi = np.arange(-0.5,0.5,0.0001)
# xi,yi = np.meshgrid(xi,yi)

# # set mask
# # mask = (xi > 0.5) & (xi < 0.6) & (yi > 0.5) & (yi < 0.6)

# # interpolate
# zi = griddata((resArr[:,0],resArr[:,1]),resArr[:,2],(xi,yi),method='linear')

# # mask out the field
# # zi[mask] = np.nan

# plot
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
plt.contourf(xi,yi,FTLEArr0)
plt.colorbar()
# plt.plot(resArr[:,0],resArr[:,1],'k,')
plt.xlabel('xi',fontsize=16)
plt.ylabel('yi',fontsize=16)
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
plt.grid()
# plt.savefig('FTLE_' + file1[:-3] + '_' + file2[:-3] + '_' + str(step) +'.png',dpi=150)
# plt.close(fig)
plt.show()

print("--- %s seconds ---" % (time.time() - startTime))