#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 19:48:19 2022

@author: akimlavrinenko
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
from NEK5000_func_lib import particleCoordsNew
import time
startTime = time.time()
# data = np.genfromtxt('/Users/akimlavrinenko/Documents/coding/data/room_data/sphere_trajectory_0.dat')
path = '/home/akim/room_production_pp/roomBackUp00001/'
path = '/Users/akimlavrinenko/Documents/coding/data/room_data/roomBackUp000001/'

# stepList = [2,1]
# fileList = ['fbalance01000.3D', 'fbalance03000.3D', 'fbalance000010.3D']
# for step in stepList:
    # for f in fileList:
        # print(f, step)
file1 = 'fbalance00001.3D'   
file2 = 'fbalance00002.3D'
t0, a0 = particleCoordsNew (path, file1)
t1, a1 = particleCoordsNew (path, file2)
a00 = np.asarray(a0[4])
a11 = np.asarray(a1[4])
step = 1

x0 = a00[::step,0]
y0 = a00[::step,1]
x1 = a11[::step,0]
y1 = a11[::step,1]

y00 = y0[-0.4 > y0 ]
y00 = y00[-0.5 < y00 ]

y0Ind = np.where(np.isin(y0, y00))[0]

x0 = x0[y0Ind]
y0 = y0[y0Ind]
x1 = x1[y0Ind]
y1 = y1[y0Ind]
x = np.asarray([x0,x1])
y = np.asarray([y0,y1])
# step = 1
# x = data[::step,1::3]
# y = data[::step,2::3]
# z = data[::step,3::3]
# time = data[::step, 0]

coordsList = []
tList = [t1,t0]

for t in range(len(tList)):
    for i in range(len(x[0,:])):
            # print(t, i)
            xTemp = x[t,:]+1
            yTemp = y[t,:]+1
            xyTemp = np.vstack((xTemp, yTemp)).T
            #x less than x0
            xLessThan0 = xyTemp[xyTemp[:,0] < xyTemp[i,0]]
            xIndLessThan0 = np.where(np.isin(xTemp, xLessThan0))
            xq02 = xyTemp[xIndLessThan0[0]]
            xyq0 = xq02[xq02[:,1]>xyTemp[i,1]]#y less than y0
            xyq2 = xq02[xq02[:,1]<xyTemp[i,1]]#y more than y0
            #x more than x0
            xMoreThan0 = xyTemp[xyTemp[:,0] > xyTemp[i,0]]
            xIndMoreThan0 = np.where(np.isin(xTemp, xMoreThan0))
            xq13 = xyTemp[xIndMoreThan0[0]]
            xyq1 = xq13[xq13[:,1]>xyTemp[i,1]]#y less than y0
            xyq3 = xq13[xq13[:,1]<xyTemp[i,1]]#y more than y0
            
            #quarter 0
            xdq0 = [(number - xyTemp[i,0])**2 for number in xyq0[:,0]]
            ydq0 = [(number - xyTemp[i,1])**2 for number in xyq0[:,1]]
            ldq0 = [xdq0,ydq0]
            sdq0 = [sum(s) for s in zip(*ldq0)] #x+y+z distance
            dq0 = [x**0.5 for x in sdq0]
            dq0 = np.asarray(dq0)
            try:
                inddminq0 = dq0.argmin()
                xcoordsq0 = xyq0[inddminq0,0]
                ycoordsq0 = xyq0[inddminq0,1]
            except:
                inddminq0 = i
                xcoordsq0 = xTemp[i]
                ycoordsq0 = yTemp[i]
            #quarter 1
            xdq1 = [(number - xyTemp[i,0])**2 for number in xyq1[:,0]]
            ydq1 = [(number - xyTemp[i,1])**2 for number in xyq1[:,1]]
            ldq1 = [xdq1,ydq1]
            sdq1 = [sum(s) for s in zip(*ldq1)] #x+y+z distance
            dq1 = [x**0.5 for x in sdq1]
            dq1 = np.asarray(dq1)
            try:
                inddminq1 = dq1.argmin()
                xcoordsq1 = xyq1[inddminq1,0]
                ycoordsq1 = xyq1[inddminq1,1]
            except:
                inddminq1 = i
                xcoordsq1 = xTemp[i]
                ycoordsq1 = yTemp[i]
            #quarter 2
            xdq2 = [(number - xyTemp[i,0])**2 for number in xyq2[:,0]]
            ydq2 = [(number - xyTemp[i,1])**2 for number in xyq2[:,1]]
            ldq2 = [xdq2,ydq2]
            sdq2 = [sum(s) for s in zip(*ldq2)] #x+y+z distance
            dq2 = [x**0.5 for x in sdq2]
            dq2 = np.asarray(dq2)
            try:
                inddminq2 = dq2.argmin()
                xcoordsq2 = xyq2[inddminq2,0]
                ycoordsq2 = xyq2[inddminq2,1]
            except:
                inddminq2 = i
                xcoordsq2 = xTemp[i]
                ycoordsq2 = yTemp[i]
            #quarter 3
            xdq3 = [(number - xyTemp[i,0])**2 for number in xyq3[:,0]]
            ydq3 = [(number - xyTemp[i,1])**2 for number in xyq3[:,1]]
            ldq3 = [xdq3,ydq3]
            sdq3 = [sum(s) for s in zip(*ldq3)] #x+y+z distance
            dq3 = [x**0.5 for x in sdq3]
            dq3 = np.asarray(dq3)
            try:
                inddminq3 = dq3.argmin()
                xcoordsq3 = xyq3[inddminq3,0]
                ycoordsq3 = xyq3[inddminq3,1]
            except:
                inddminq3 = i
                xcoordsq3 = xTemp[i]
                ycoordsq3 = yTemp[i]
                
            coordsList.append(np.asarray([xcoordsq0, ycoordsq0, xcoordsq1, ycoordsq1, xcoordsq2, ycoordsq2, xcoordsq3, ycoordsq3]))


    plt.plot(xTemp, yTemp, 'r,')
    # plt.plot(xTemp[i], yTemp[i], 'ko')
    # plt.plot(xyq0[:,0], xyq0[:,1], 'bo')
    # plt.plot(xyq1[:,0], xyq1[:,1], 'yo')
    # plt.plot(xyq3[:,0], xyq3[:,1], 'go')
    # plt.plot(xyq2[:,0], xyq2[:,1], 'ro')
    # plt.plot(xcoordsq0, ycoordsq0, 'mo')
    # plt.plot(xcoordsq1, ycoordsq1, 'mo')
    # plt.plot(xcoordsq2, ycoordsq2, 'mo')
    # plt.plot(xcoordsq3, ycoordsq3, 'mo')

    plt.grid()
    plt.show()


tCoordsList = [coordsList[z:z+(len(tList))] for z in range(0, len(coordsList), len(tList))]

t0 = 1
t1 = 0
ftleList = []

for i in range(len(tCoordsList)):
    xdeltaT0q12 = abs(tCoordsList[i][t0][2] - tCoordsList[i][t0][4])
    xdeltaT1q12 = abs(tCoordsList[i][t1][2] - tCoordsList[i][t1][4])
    
    xdeltaT0q03 = abs(tCoordsList[i][t0][6] - tCoordsList[i][t0][0])
    xdeltaT1q03 = abs(tCoordsList[i][t1][6] - tCoordsList[i][t1][0])
    
    ydeltaT0q12 = abs(tCoordsList[i][t0][3] - tCoordsList[i][t0][5])
    ydeltaT1q12 = abs(tCoordsList[i][t1][3] - tCoordsList[i][t1][5])
     
    ydeltaT0q03 = abs(tCoordsList[i][t0][1] - tCoordsList[i][t0][7])
    ydeltaT1q03 = abs(tCoordsList[i][t1][1] - tCoordsList[i][t1][7])
    
    DF = np.zeros((2,2))
    DF[0,0] = xdeltaT1q12/xdeltaT0q12
    DF[0,1] = xdeltaT1q03/ydeltaT0q03
    DF[1,0] = ydeltaT1q12/xdeltaT0q12
    DF[1,1] = ydeltaT1q12/ydeltaT0q12
    delta = np.dot(DF,DF.T)
    lambdaMax = max(np.linalg.svd(delta)[1])
    FTLE = 1/(tList[0]-tList[1]) * np.log((lambdaMax**0.5))
    ftleList.append(FTLE)
    

resArr = np.zeros((len(x[0,:]),3))
resArr[:,0] = x[0]
resArr[:,1] = y[0]
resArr[:,2] = np.asarray(ftleList)


# target grid to interpolate to
xi = yi = np.arange(-0.5,0.5,0.0001)
xi,yi = np.meshgrid(xi,yi)

# set mask
# mask = (xi > 0.5) & (xi < 0.6) & (yi > 0.5) & (yi < 0.6)

# interpolate
zi = griddata((resArr[:,0],resArr[:,1]),resArr[:,2],(xi,yi),method='linear')

# mask out the field
# zi[mask] = np.nan

# plot
fig = plt.figure(figsize=(10,2))
ax = fig.add_subplot(111)
plt.contourf(xi,yi,zi)
plt.colorbar()
# plt.plot(resArr[:,0],resArr[:,1],'k,')
plt.xlabel('xi',fontsize=16)
plt.ylabel('yi',fontsize=16)
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, -0.4)
plt.grid()
plt.savefig('FTLE_' + file1[:-3] + '_' + file2[:-3] + '_' + str(step) +'.png',dpi=150)
# plt.close(fig)
plt.show()

print("--- %s seconds ---" % (time.time() - startTime))
