import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time


def glacier():
    D={}
    ## --------------------------------
    #define parameters and constants
    #----------------------------------
    rho=917.0              # density of ice in kg/m^3
    g=9.8                  # gravity in m/s^2
    n=3.0                    # glen's flow law constant
    
    s_per_year = 365.35*24*3600        # number of seconds in a year (s yr^-1)
    A_T = 6.8e-24*s_per_year       # softness parameter for -5 degree ice (yr^-1 Pa^-3)
    fd = 2*A_T/(n+2)                           # constants lumped together, and put into standard units.  (Pa^-3 yr^-1)
    fs = 5.7e-20*s_per_year                # sliding parameter fromOerlemans (Pa^-3 m^2 yr-1)
    
    ## --------------------------------
    #define grid spacing and time
    #----------------------------------
    xmx = 20000.          # max of the domain size in m
    delx = 100              # grid spacing in m
    nxs = np.round(xmx/delx) + 1  # number of grid points
    x = np.arange(0,xmx+delx,delx)                   # x array (each point)
    
    delt = 0.0125 # time step in yrs 
    ts = 0             # starting time in yrs
    tf = 500      # final time in yrs
    nts=np.floor((tf-ts)/delt) + 1 # number of time steps ('round' used because need nts as integer)
    nyrs = tf-ts       # just final time - start time
    
    ## -----------------
    # define accumulation/ablation
    #-----------------
    
    # climate foricing for steady state
    b = 3-(4./5000)*x      # mass balance in /yr
    
    ## ---------------------------------
    # glacier bed geometries
    #---------------------------------
    delb = 15.0e3/np.log(3)
    zb =  3000.*np.exp(-x/delb)     # bed profile in m  
    #zb = (-0.1*x+1000)  
    
    ## -----------------------
    # initialize arrays
    #------------------------
    # load glac_init_class.mat
    # H0  = interp1(x_init,h_init,x,'linear')
    
    H0 = np.ones_like(x)     # initial thickness of ice (in m)--do not confuse with h, ice elevation.  h = H+zb
    Qp=np.zeros_like(x)      # Qp equals jj+1/2 flux
    Qm=np.zeros_like(x)      # Qm equals jj-1/2 flux
    up = np.zeros_like(x)    # velocity at jj +1/2 grid 
    um = np.zeros_like(x)    # velcotiy at jj-1/2 grid
    dHdt=np.zeros_like(x)        #  slope of ice thickness -- do not confuse with dhdx, surface slope.  dhdx = dzdx+dHdx
    H = H0                    # initialize height array
    
    yr = 0         # for counting years
    idx_out = 0    # 
    deltout = 1    # for counting years
    nouts = np.round(nts/1)   # for counting years
    edge_out = np.zeros_like(nouts)    # for counting
    t_out = np.zeros_like(nouts)           # for counting
    
    ## -----------------------------------------
    # begin loop over time
    #-----------------------------------------
    
    for ii in range(int(nts)): 
        
        
        t = delt*(ii) # time in years
        if ii<=2 or ii>=int(nts)-1:
            print ii
            print t 
            
        #-----------------------------------------
        # begin loop over space
        #-----------------------------------------
        for jj in range(int(nxs)-1):  #
        #for jj in range(int(nxs)):  #  
            
            if jj==0:      # at the glacier divide

                H_ave =(H[0] + H[1])/2          # average glacier thickness between divide and second point
                dHdx = (H[1] - H[0])/delx      # linear slope of glacier thickness between divide and second point
                dzbdx = (zb[1] - zb[0])/delx     # linear slope of bed between divide and second point
                dhdx = dHdx+dzbdx             # linear surface slope of ice between divide and second point
            
                Qp[0] = -(dhdx)**3 * H_ave**4*(rho*g)**3*(fd*H_ave+fs/H_ave)       # flux at plus half grid point
                up[0] = Qp[0]/H_ave    # velocity at half grid point
            
                Qm[0] = 0                      # flux at minus half grid point = 0 because we start at the divide--nothing flowing in
                um[0] = Qm[0]/H[0]       # velocity is also zero coming in
            
                dHdt[0] = b[0] - Qp[0]/(delx/2)    # change in thickness at the divide for this timestep is based on snowfall in/out, flux out
        
            elif H[jj]==0 and H[jj-1]>0: # at the glacier toe
        
                Qp[jj] = 0   # nothing to flux out of the end
                up[jj] = 0   # no velocity b/c no ice
        
                H_ave = (0+H[jj-1])/2              # thickness of glacier is average of  last point with thickness and a thickness of 0     
                dHdx = 0-H[jj-1]/delx              # slope is based on last point with thickness and a thickness of 0
                dzbdx =(zb[jj]-zb[jj-1])/delx       # bed slope between toe and last point with thickness
                dhdx = dHdx+dzbdx            # surface slope between toe and last point with thickness
        
                Qm[jj] = -(rho*g)**3*H_ave**4*(dhdx)**3*(fd*H_ave+fs/H_ave)    # flux coming from the last point with thickness and melting out before it reaches a new grid point
                um[jj] = Qm[jj]/H_ave         # velocity at the toe
        
                dHdt[jj] = b[jj] -(Qp[jj]- Qm[jj])/delx     # change in height of the glacier at the toe based on the snowfall in/melt out and the flux in/out
                edge = jj 	#index of glacier toe - used for fancy plotting
        
            elif H[jj]<=0 and H[jj-1]<=0: # beyond glacier toe - no glacier flux

                dHdt[jj] = b[jj]  # thickness change is only from any change in accumulaiton in/melt out
        
                # no flux in, no flux out
                Qp[jj] = 0 
                Qm[jj] = 0
        
                um[jj] =0 
                up[jj] = 0
            
            else:  # within the glacier
                    
                # first for the flux going out
                H_ave = (H[jj+1] + H[jj])/2              # the average thickness between jj and jj+1
                dHdx = (H[jj+1] - H[jj])/delx           # the linear ice thickness slope between jj and jj+1
                dzbdx = (zb[jj+1] - zb[jj])/delx       # the linear bed slope between jj and jj+1
                dhdx = dHdx + dzbdx                      # the surface slope between jj and jj+1
        
                Qp[jj] = -(rho*g)**3*H_ave**4*(dhdx)**3*(fd*H_ave+fs/H_ave)  # flux coming out
                up[jj] = Qp[jj]/H_ave                           # velocity coming out
                
                # now for the flux coming in
                H_ave = (H[jj-1] + H[jj])/2               # the average thickness between jj and jj-1
                dHdx = (H[jj] - H[jj-1])/delx            # the linear ice thickness slope between jj and jj-1
                dzbdx = (zb[jj] - zb[jj-1])/delx       # the linear bed slope between jj and jj-1
                dhdx = dHdx + dzbdx                      # the surface slope between jj and jj-1

                Qm[jj] = -(rho*g)**3*H_ave**4*dhdx**3*(fd*H_ave+fs/H_ave)        # flux coming in
                um[jj] = Qm[jj]/H_ave            # velocity coming in
                
                # change in height at point jj for this timestep
            dHdt[jj] = b[jj] - (Qp[jj] - Qm[jj])/delx       # based on accumulation in/melt out and flux in/out
            
            #end  
            
            dHdt[nxs-1] = 0 # enforce no change at edge of domian
        
        
        #end  # done with each point in space
            
        # For the next timestep, the thickness is equal to the thickness profile from the previous timstep, plus the amount that the thickness has changed at each point
        # over this time period, (dHdt*delt).  The ice thickness at each point
        # needs to be positive, so we choose the max of 0 and the summation.
        # since dHdt can be negative, we want to be sure we don't have values
        # for a negative ice thickness.
        
        H = np.maximum(0 , (H + (dHdt*delt)))
        H[-1]=0.0
            
        #end  # done with each time step
#     D['H']=H
#     D['x'
    Hinit = H
    return Hinit, x, Qm,Qp,um,up,zb

if __name__ == "__main__":
    tic=time.time()
    Hinit, x, Qm,Qp,um,up,zb = glacier()
    H=Hinit
    ###
    #plt.figure(1)        
    #
    #subplot(221)  plot(x,Qm,'co') ylabel('flux')title('western cv') 
    #subplot(222)  plot(x,Qp,'co') title('eastern cv')
    #subplot(223)  plot(x,um,'co') ylabel('velocity') xlabel('distance (m)') 
    #subplot(224)  plot(x,up,'co') xlabel('distance (m)') 
    #
    #
    plt.figure(2)
    plt.plot(x,H,'c*')
    
    #
    plt.figure(3)
    plt.plot(x,H + zb,'c')
    plt.plot(x,zb,'k')
    
    plt.show()

    elapsed=time.time()-tic
    elapsed_min=elapsed/60.
    mins=np.floor(elapsed_min)
    secs=(elapsed_min-mins)*60 

    #print "run took %s minutes %s seconds" % mins, secs

