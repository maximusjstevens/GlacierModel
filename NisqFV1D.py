# Nisqually Glacier Model
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

np.set_printoptions(linewidth=600) #set output reading to be wider. A bit easier to read

#Set Globals
rho_i = 917.0 #kg/m^3
R = 8.314 #J/mol/K
M_air = 28.97e-3 #kg/mol
rho_bco = 815. #kg/m^3
p_0 = 1.01325e5 # Standard Amtmospheric Pressure, Pa
T_0 = 273.15 # Standard Temp, K
sPerYear = 365.25*24*3600 #seconds per year
g=9.8 #m/s^2

# Downward Advection (of air)        
def w(): # Function for downward advection of air and also calculates total air content. 
    pass
    
def A(P): # Power-law scheme, Patankar eq. 5.34
    A = np.maximum( (1 - 0.1 * np.abs( P ) )**5, np.zeros(np.size(P) ) )
    return A    
    
def F_upwind(F): # Upwinding scheme
    F_upwind = np.maximum( F, 0 )
    return F_upwind

def solver(a_U,a_D,a_P,b): #routine to solve Ax=b
    nz=np.size(b)
    Diags = (np.append([a_U,-a_P],[a_D],axis=0))
    cols=np.array([1, 0, -1])      
    big_A=spdiags(Diags,cols,nz,nz,format='csc')
    big_A=big_A.T        
    rhs=-b    
    H_t=splin.spsolve(big_A,rhs)    
    return H_t
    
def gamma_phi(phi,phi_ref,gamma_ref,gamma_grad):
    Gamma_phi = ( gamma_ref * phi_ref ) / phi
    return Gamma_phi
    
def glacier():
    
    conv_test = 0.01
    max_iter = 20
    
    length = 10000.0 # m
    z_res=10.0 #resolution of grid, m

    tend=500.
    time_yr=np.array([0, tend]) # Atmospheric measurements times            
    time_yr_s=time_yr*sPerYear    
    yrs=np.around(time_yr[-1]-time_yr[0])
    time_total=yrs*sPerYear #total model run time in seconds
    stpsperyear=1.0 #If this is for transient this number must (for now) be the same time steps as the input density/depth files. Make sure that it is a float.
    t_steps=np.int(yrs*stpsperyear)
    dt=time_total/t_steps #time step size. 
    model_time=np.arange(np.around(time_yr[0]*sPerYear),np.around(time_yr[-1]*sPerYear),dt) #set model time steps
    model_time_years=model_time/sPerYear
    
    nt=np.size(model_time) #number of time steps
    
    dz, z_edges, z_nodes, nodes, nz_P, nz_fv = space(length,z_res) #call on space function to set up spatial grid
    
    delt = 15e3/np.log(3)
    bed =  length*np.exp(-z_nodes/delt)    # bed profile in m  
    bedgrad = np.gradient(bed,dz)
    print 'bedgrad=',bedgrad

    bc_u_0 = 0.0 #this is at head of glacier
    bc_type = 2
    bc_u   = np.concatenate( ([bc_u_0], [bc_type]))

    bc_d_0 = 0.0 #toe of glacier
    bc_type = 1
    bc_d   = np.concatenate(([ bc_d_0 ], [ bc_type ]))
    
    H_0 = 1.*np.ones(len(z_nodes)) 
#     H_0[:]=bc_u_0   
    
    fd=1.9e-24
    fs=5.7e-20
    B = (22.0-(40.0/5000)*z_nodes)/sPerYear/(t_steps/yrs) #mass balance
    print "B=",B[0]          
    #B = 0.                                                         
    for i_time in range(0,nt): 

#         diffu = diffusivity() #get diffusivity profile
        

        if i_time==0:
            s=(nt,nz_P)
            H=np.zeros(s)
            diffu_hold=np.zeros(s)
            H_t=H_0
#             with open(ConcPath, "w") as f:
#                 writer = csv.writer(f)       
#                 writer.writerow(np.append(model_time_years[i_time],H_t))
#             with open(RhoPath, "w") as f:
#                 writer = csv.writer(f)       
#                 writer.writerow(np.append(model_time_years[i_time],rho_prof))
#             with open(DiffuPath, "w") as f:
#                 writer = csv.writer(f)       
#                 writer.writerow(np.append(model_time_years[i_time],diffu))
#             with open(ZPath, "w") as f:
#                 writer = csv.writer(f)       
#                 writer.writerow(np.append(model_time_years[i_time],z_nodes))
        h = bed + H_t
        hgrad = np.gradient(h,dz)
        Hgrad = np.gradient(H,dz)
        #print 'i_time=',i_time
        if i_time == 0:
            print "nt=",nt
            print "bedgrad=",bedgrad
            print "hgrad=",hgrad
            print "Hgrad=",Hgrad
        
        diffu = (rho_i*g)**3 * H_t**3 * (hgrad)**2 * (fd*H_t**2+fs)
        #diffu = -1 * (rho_i*g)**3 * 100**3 * (bedgrad)**2 * (fd*100**2+fs)
        H[i_time,:] = H_t
        
        
        diffu_P=diffu
        diffu_hold[i_time,:]=diffu_P
        dZ = np.concatenate(([1],np.diff(z_edges),[1]))
        
        dZ_u = np.diff(z_nodes)
        dZ_u = np.append(dZ_u[0], dZ_u)
        #dZ_u = np.gradient(z_nodes)
        
        dZ_d = np.diff(z_nodes)
        dZ_d = np.append(dZ_d,dZ_d[-1])
        #dZ_d = np.gradient(z_nodes)
    
        f_u = np.append(0, (1 -(z_nodes[1:] - z_edges)/dZ_u[1:]))
        f_d = np.append(1 - (z_edges - z_nodes[0:-1])/dZ_d[0:-1], 0)
        
        diffu_U = np.append(diffu_P[0], diffu_P[0:-1] )
        diffu_D = np.append(diffu_P[1:], diffu_P[-1])
    
        diffu_u =  1/ ( (1 - f_u)/diffu_P + f_u/diffu_U )
        diffu_d =  1/ ( (1 - f_d)/diffu_P + f_d/diffu_D )
        
    
        S_C_0 =  - B*dt + (-diffu_d+diffu_u)*np.gradient(bed,dz) #-1*B?
        xx=(-diffu_d+diffu_u)*np.gradient(bed,dz)
        print -B[0]*dt,xx[0]
    
        S_C=S_C_0 
        
        #S_P=(-diffu_d+diffu_u)*np.gradient(bed,dz) #gravity term, S_P is H-dependent source
        S_P=0.0
        
        b_0 = S_C*dZ
        
#         rho_interface=np.interp(z_edges,z_nodes,rho_prof)
        
#         w_edges = w()

#         w_u = np.append(w_edges[0],  w_edges )
#         w_d = np.append(w_edges, w_edges[-1])
        
        D_u = ((diffu_u) / dZ_u) #check signs
        D_d = ((diffu_d) / dZ_d)
    
        
#         F_u =  w_u #Is this correct?
#         F_d =  w_d 
        
#         P_u = F_u/ D_u
#         P_d = F_d/ D_d
        
#         a_U = D_u * A( P_u ) + F_upwind(  F_u )
#         a_D = D_d * A( P_d ) + F_upwind( -F_d )

        a_U = - D_u # *-1 if using dH/dt = -d/dx(D dH/dx)
        a_D = - D_d
    
        a_P_0 = dZ/dt
                                
        a_P = a_P_0 + a_U + a_D + S_P*dZ
                   
        b = b_0 + a_P_0*H_t #+ a_U*-1*bedgrad + a_D*bedgrad      
        
        #Upper boundary
        a_P[0] = 1 
        a_U[0] = 0
        a_D[0] = 1
        b[0]= -1*dZ_u[1]*bc_u[0]
        
        #Down boundary
        a_P[-1] = 1 
        #a_U[-1] = 0
        a_D[-1] = 0
        a_U[-1] = 0
        #b[-1]=dZ_d[-1]*bc_d[0]
        b[-1]=bc_d[0]
        

        H_t = solver(a_U,a_D,a_P,b)
        
        H_t[H_t<0.0]=0.0
        

#         rho_hold[i_time,:]=rho_prof
        
#         if i_time>0:        
#             with open(ConcPath, "a") as f:
#                 writer = csv.writer(f)       
#                 writer.writerow(np.append(model_time_years[i_time],H_t))
#             with open(RhoPath, "a") as f:
#                 writer = csv.writer(f)       
#                 writer.writerow(np.append(model_time_years[i_time],rho_prof))
#             with open(DiffuPath, "a") as f:
#                 writer = csv.writer(f)       
#                 writer.writerow(np.append(model_time_years[i_time],diffu))
        
    return H, z_nodes, bed, diffu_hold
    
def space(length,z_res):
    Lz=length
    nz_fv=np.around(length/z_res)
    nz_P=nz_fv+2
    
    dz=Lz/nz_fv
    z_edges=dz*np.arange(0,nz_fv+1)
        
    z_nodes=np.concatenate(([z_edges[0]],z_edges[0:-1]+np.diff(z_edges)/2,[z_edges[-1]]))
    nodes = np.size(z_nodes)
    
    return dz, z_edges, z_nodes, nodes, nz_P, nz_fv
            
#def glacierFV1D():
#    
#    d={}
#    
#    H, nodes, bed = glacier()
#            
#    d[xx]=H
#    d[yy]=nodes
#            
#    return d
    
    
if __name__ == "__main__":
    

    import time
    
    tic=time.time()

    #d = glacierFV1D()
    
    H, nodes, bed, diffu_hold = glacier()
                
    elapsed=time.time()-tic
    elapsed_min=elapsed/60.
    mins=np.floor(elapsed_min)
    secs=(elapsed_min-mins)*60
    print mins, 'min', secs, 'sec elapsed'
    
    
    plt.figure(1)
    plt.plot(nodes,bed)
    plt.show()
        