# Nisqually Glacier Model, finite difference scheme
# Max Stevens
# version 0.1, 29 July 2015


### In order to run transient, you must put the files from firnmodel.py output into DataImport folder.
### All units should be kilograms, meters, and seconds (MKS system)

import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.sparse.linalg as splin
from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr
from scipy.integrate import cumtrapz
import math
import csv
import json

# Set path to find all files to import and set output location
spot = os.path.dirname(sys.argv[0]) #Add Folder
print spot 
# os.chdir(spot) #change pwd to location of firnmodel.py
sys.path.insert(1,spot)
# ResultsPlace=os.path.join(spot,'Results')
# sys.path.insert(3,os.path.join(spot,'DataImport'))
# DataPath = os.path.join(spot,'DataImport')

#DataPathUser = '/Users/maxstev/documents/grad_school/pire/CFM/CommunityFirnModel/code/firnmodel/layers'

np.set_printoptions(linewidth=1200) #set output reading to be wider. A bit easier to read

#Set Globals
rho_i = 917.0 #kg/m^3
R = 8.314 #J/mol/K
M_air = 28.97e-3 #kg/mol
rho_bco = 815. #kg/m^3
p_0 = 1.01325e5 # Standard Amtmospheric Pressure, Pa
T_0 = 273.15 # Standard Temp, K
sPerYear = 365.25*24*3600 #seconds per year
g=9.8 #m/s^2

#def glacier():
    
length=10000.0 #m
dx = 50.0
nodes = np.arange(0,length+dx,dx)
nx = np.size(nodes)

dref=1.0
stabdt=dx**2/(4*dref)

tstart=0.0 #years
speryear=365.25*24*3600
tendyrs=100.0
tend=tendyrs*speryear 
dt=stabdt/2
t=np.arange(tstart,tend+dt,dt)
nt=np.size(t)

delt = 15e3/np.log(3) #creating bed profile
bed =  length*np.exp(-nodes/delt)    # bed profile in m  
#bedgrad = np.gradient(bed,dx)    

H_0 = 1.*np.ones(len(nodes))
B_yrs = (22.0-(40.0/5000)*nodes)
B=B_yrs/speryear
n=3.0
fd=1.9e-24
fs=5.7e-20


for i_time in range(0,nt):
    #print i_time
    tt=np.arange(10,100,10)
    
    if i_time==0:
        s=(nt,nx)
        H=np.zeros(s)
        diffu_hold=np.zeros(s)
        H_t=H_0
        glac=np.nonzero(H_t>0)
        bb=np.nonzero(H_t<=0)
        terminus=np.max(glac)
        ii=np.zeros(nt)
        H[i_time,:] = H_t
        H_t[0]=H_t[1]
        H_t[H_t<0]=0
        
    ii[i_time]=i_time

    h = bed + H_t
    hup=np.append(h[0],h[0:-1])
    hdown=np.append(h[1:],h[-1])
    
    C=(fd*H_t**2+fs)*(rho_i*g)**n
    
    #C=2*A/(n+2)*(rho_i*g)**n
    #hgrad=(hup+hdown)
    #D = C * H_t**n * (np.abs(hgrad)/(2*dx))**n-1
    D = C * H_t**n * (np.abs(np.gradient(h,dx)))**n-1
    print i_time,D[0:10]
    D[terminus+1:]=D[terminus]
    D[D>0]=0.0
    
    diffu_hold[i_time,:]=D
    
    Dup=np.append(D[0],D[0:-1])
    Ddown=np.append(D[1:],D[-1])
    
    
    Fdown = 0.5*(D+Ddown)*((hdown-h)/dx)
    Fup = 0.5*(D+Dup)*((h-hup)/dx)
    
    H_t = H_t + dt/dx*(Fdown-Fup)+B*dt
    
    glac=np.nonzero(H_t>0)
    bb=np.nonzero(H_t<=0)
    terminus=np.max(glac)
    H_t[0]=H_t[1]
    H_t[H_t<0]=0
    
    H[i_time,:] = H_t
        
        #if i_time in tt:
        #    print "i_time=",i_time
        #    plt.figure(1)
        #    plt.clf()
        #    plt.plot(nodes,bed,'k')
        #    plt.plot(nodes,h,'b')
        #    plt.show()
            
            

        
    
    #return H, nodes, bed, diffu_hold, ii
    
    
    
#if __name__ == "__main__":
#    import time
#    
#    tic=time.time()
#
#    #d = glacierFV1D()
#    
#    H, nodes, bed, diffu_hold, ii = glacier()
#                
#    elapsed=time.time()-tic
#    elapsed_min=elapsed/60.
#    mins=np.floor(elapsed_min)
#    secs=(elapsed_min-mins)*60
#    print mins, 'min', secs, 'sec elapsed'
#    
#    
#    plt.figure(1)
#    plt.plot(nodes,bed)
#    plt.show()
#
