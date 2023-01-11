#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 14:38:40 2023

@author: deo010
"""

#%%

import numpy as np
import math
import os
from matplotlib import pyplot as plt
from scipy.signal import find_peaks

#%% Plotting graphs

def plotgraphs(filename):
    
    
    file = np.genfromtxt("{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,2:].T
    phase = file[0,:]
    flux = file[1,:]
    
    plt.figure(figsize=(20,12))
    plt.style.use('dark_background')
    plt.plot(phase,flux,'w')
    plt.xlabel('Time (ms)')
    plt.ylabel('Flux Density (mJy)')
    plt.title('{}'.format(filename))


#%% Find the noisy channels

def noise(filename):
    file = np.genfromtxt("{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,2:].T
    phase = file[0,:]
    flux = file[1,:]
    
    length = len(phase)
    stds = []
    stdbins = []
    bins = 3
    
    for i in range(bins):                                                         
        fluxbins = flux[(i*int(length/bins)):((i+1)*int(length/bins))]         # create bins
        phasebins = phase[(i*int(length/bins)):((i+1)*int(length/bins))]       # pick out bins of phase    
        std = np.std(fluxbins)                                                 # take std of each bin
                                       
        stds += [std]
        stdbins += [phasebins]                  
        
    return stdbins[np.argmin(stds)]

#%%


noise('t160114_192200.sfARCH.ar.zap.Fp')

#%% Calculate experimental rms

def rmsprac(filename,noise):
    
    ######### Finding RMS #########
    
    # the equation used is RMS = sqrt [ 1/n * sum (x_i) for i values ]  
    
    file = np.genfromtxt("{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,3:].T
    I = file[0,:]
    
    start = int(min(noise))
    end = int(max(noise))
    obs = end-start
    
    
    rmsp = math.sqrt( (1/obs)*np.sum( (I[start:end])**2 ) )
    
    print("The practical rms is {}".format(rmsp))
    
    return rmsp, obs

#%% Calculate theoretical rms

def rmstheory(beta, SEFD, np, t_int, del_f):
    
    ######### Calculating RMS #########
    
    # Equation used is radiometer equation
    
    beta = float(beta)                                                         # correction factor
    SEFD = float(SEFD)                                                         # milliJansky, calculated using a fold-mode observation of J1550-5418 taken on 2015-03-12-22:14:02 UTC
    np = float(np)                                                             # polarisation number summed
    t_int = float(t_int)                                                       # integration time per phase bin
    del_f = float(del_f)                                                       # bandwidth
    
    rmst = (beta*SEFD)/(math.sqrt(np*t_int*del_f))
    print("The theoretical rms is {}".format(rmst))
    
    return rmst

#%% Find the correction factor

def correction(filename, rmst, rmsp):
    factor = rmst/rmsp
    print("The correction factor is {}.".format(factor))
    os.system("awk '{{if (NR!=1) {{{{ $4=$4*{} }}}} print }}' {} > cal_{}".format(factor,filename,filename))

#%% Find a window for the width of the fluence

def find_nearest(array,value):
    array = np.asarray(array)
    idx = np.searchsorted(array, value, side="left")
    if array[idx] < value:
        return array[idx], array[idx+1]
    else:
        return array[idx-1],array[idx]
    
    # if value < idx:
        # return array[idx], array[idx]

#%%

def checkless(list1, val):
    index = []
    value = []
    
    for i,x in enumerate(list1):
        if x < val:
            index.append(i)
            value.append(x)
    
    return index, value

#%% Finding peaks

def peakfinds(filename,prominence):
    file = np.genfromtxt("cal_{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,2:].T
    phase = file[0,:]
    flux = file[1,:]
    
    trough = -flux
    
    peaks, _ = find_peaks(trough,prominence = 100)
    
    print(trough)

    # maxi = max(flux)
    # arr = np.argmax(flux)
    
    # left,right = find_nearest(peaks,arr)

    
    plt.figure(figsize=(20,12))
    plt.plot(phase,-flux,'w')
    plt.plot(peaks,-flux[peaks],'o',color='red')
    # plt.plot(arr,maxi,'o',color='red')
    
    # return left, right
    

#%% Finding troughs
def troughfinds(filename,sig):
    file = np.genfromtxt("cal_{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,2:].T
    phase = file[0,:]
    flux = file[1,:]
    
    minind,minval = checkless(flux,-sig)
    maxi = max(flux)
    arr = np.argmax(flux)
    
    minim, maxim = find_nearest(minind,arr)
    print(minim)
    print(maxim)
    
    plt.figure(figsize=(20,12))
    plt.plot(phase,flux,'w')
    plt.plot(minind,minval,'o',color='red')
    plt.plot(arr,maxi,'o',color='red')
    plt.xlabel('Time (ms)')
    plt.ylabel('Flux Density (mJy)')
    plt.title('{}'.format(filename))
    
    return minim, maxim
    #plt.plot(peaks,flux[peaks],'o',color='red')
    
#%%

troughfinds('t160117_031649.sfARCH.ar.zap.Fp.txt',35)

#%% Find fluence

def window(filename,start,end):
    file = np.genfromtxt("cal_{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,2:].T
    phase = file[0,:]
    flux = file[1,:]

    
    burst = end-start
    fluence = (np.sum(flux))/1000
    
    phase = file[0,:]
    
    
    plt.figure(figsize=(20,12))
    plt.style.use('dark_background')
    plt.plot(phase,flux,'w')
    plt.axvline(start,color='red')
    plt.axvline(end,color='red')
    plt.xlabel('Time (ms)')
    plt.ylabel('Flux Density (mJy)')
    plt.title('{}'.format(filename))
    
    print("Burst width is {} ms".format(burst))
    print("Fluence is {} Jy ms".format(fluence))
    
    return fluence

#%% Big function
def cal_fluence(filename):
    
    ##### CALIBRATION #####
    
    beta = 1.1
    SEFD = 79000
    np = 2
    t = 1.5
    del_f = 1024e6
    
    
    # Find which channels do not contain the peak
    noisetime = noise(filename)
    # Find the practical RMS 
    rmsp,obs = rmsprac(filename,noisetime)
    # Find the theoretical RMS
    rmst = rmstheory(beta, SEFD, np, t/obs, del_f)
    # Create text file with calibrated data
    correction(filename, rmst,rmsp)
    # Plot calibrated plot
    plotgraphs('cal_{}'.format(filename))
    
    ##### FLUENCE #####
    # Find the peak window
    start, end = troughfinds(filename,rmst)
    # Find the fluence of the burst
    window(filename,start,end)
    
    
    
#%% 
cal_fluence('t160114_192200.sfARCH.ar.zap.Fp.txt')
    
#%%

cal_fluence("t160117_031649.sfARCH.ar.zap.Fp.txt")

#%%

cal_fluence("t160120_184739.sfARCH.ar.zap.Fp.txt")

# maxphase = peaks[np.argmax(flux[peaks])]
# maxflux = flux[maxphase]

# print(peaks)
# print(maxflux)
# print(maxphase)

# troughs, _ = find_peaks(-flux,prominence = prominence)
#%%




























  
