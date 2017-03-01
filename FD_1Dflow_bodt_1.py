### 1-D Glacier Flow Model based on continuity: dH/dt = dQ/dx + b_dot
### Finite Difference, explicit time steps
### Python version coded by Max Stevens, 8/14/15, based on (copied from) Matlab model from Gerard Roe/Kat Huybers/Summer Rupper
### References: Oerlemans 2001: Glaciers and climate change. Page 61. 

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time
import matplotlib
import os
import sys
import pandas as pd
matplotlib.rcParams['axes.formatter.useoffset']=False

spot = os.path.dirname(sys.argv[0]) #Add Folder

def glacier(D,b_in,delx,nxs,x,zb,bedprof,nyrs,spin,Hinit=None):

    #D={}
    ## -------------------------------------
    # define parameters and constants
    ## -------------------------------------
    rho=917.0              # density of ice in kg/m^3
    g=9.8                  # gravity in m/s^2
    n=3.0                    # glen's flow law constant
    mu = 0.65;      # melt rate in m /yr /degC
    gam = 5.5e-3;     # atmospheric lapse rate in degC/m

    s_per_year = 365.25*24*3600        # number of seconds in a year (s yr^-1)
    A_T = 6.8e-24*s_per_year       # softness parameter for -5 degree ice (yr^-1 Pa^-3). Units given in Oerlemans are per second, so need sPerYear conversion
    fd = 2*A_T/(n+2)                           # constants lumped together, and put into standard units.  (Pa^-3 yr^-1)
    fs = 5.7e-20*s_per_year                # sliding parameter fromOerlemans, page 62 (Pa^-3 m^2 yr-1). Units given in Oerlemans are per second

    #delt = 0.0004 # time step in yrs. Stable size is ~delx^2 / R, where R is ~3-10e5. Seems like the smaller delx is, the smaller R must be.
    dfa=10.e5 
    delt = delx**2/dfa #time step in years, varies based on grid size
    print delt
    ts = 0             # starting time in yrs
    if spin:
        tf = 100      # final time in yrs
    else:
        tf=nyrs #number of years of model run
    nts=np.floor((tf-ts)/delt) + 1 # number of time steps ('round' used because need nts as integer)
    nyrs = tf-ts       # just final time - start time

    ## ------------------------------------
    # initialize arrays
    ## ------------------------------------

    if Hinit is None:
        H0 = 20.0*np.ones_like(x)     # initial thickness of ice (in m)--do not confuse with h, ice elevation.  h = H+zb
        H0[-1]=0.0                 # impose 0 thickness at end of domain.
    else:
        H0 = Hinit
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
    
    #if not spin:
    #    t_hold=np.zeros((nts))
    #    L_hold=np.zeros_like(t_hold)
    #    H_hold=np.zeros((nts,nxs))
    #    Qp_hold=np.zeros_like(H_hold)
    #    Qm_hold=np.zeros_like(H_hold)
    #    up_hold=np.zeros_like(H_hold)
    #    um_hold=np.zeros_like(H_hold)
    #    dHdt_hold=np.zeros_like(H_hold)
    #    b_hold=np.zeros_like(H_hold) 
        
    if not spin:
        t_hold=np.zeros((nyrs+1))
        L_hold=np.zeros_like(t_hold)
        H_hold=np.zeros((nyrs+1,nxs))
        Qp_hold=np.zeros_like(H_hold)
        Qm_hold=np.zeros_like(H_hold)
        up_hold=np.zeros_like(H_hold)
        um_hold=np.zeros_like(H_hold)
        dHdt_hold=np.zeros_like(H_hold)
        b_hold=np.zeros_like(H_hold)
        ELA_hold=np.zeros_like(t_hold)
        bnet_hold=np.zeros_like(t_hold)
        Tyr_hold=np.zeros_like(t_hold) 
        Pyr_hold=np.zeros_like(H_hold)  

    ## ------------------------------------
    # begin loop over time
    ## ------------------------------------
    yy=0 #for tracking
    
    Pmean = 10.0*np.cos(np.pi/180.0*(1/45.0*x-30.0))+1.0
    modelyear=1900.0
    for ii in range(int(nts)):
        seasontime=modelyear-np.floor(modelyear)
  
        t = delt*(ii) # time in years

        if spin:
            #ela=3000.0
            #icsu=H+zb
            #b=(icsu-ela)/116.
            b=b_in
            
        if not spin and t==np.floor(t):
            #b = 0.004*zb-11.236
            b=b_in[yr,:]
            #ela=3000.0
            #icsu=H+zb
            #b=(icsu-ela)/116.
            #if seasontime>0.5 and seasontime<0.83:

            #Pyr=Pmean+Pp[yr] #add noise to mean precip
            #P = np.ones_like(x)*(Pyr) #add that noise over whole glacier
            #Tyr=Tmean+Tp[yr] #add noise to mean temp
            #T_wk    = (Tyr)*np.ones_like(x) - gam*(zb+H);    # add the noisy temp to the whole glacier and adjust for lapse rate
            #melt = np.maximum(0,mu*T_wk);      # convert the temeperature to melt/ablation with meltrate factor--don't "negative" melt
            #b = P-melt;     # accumulation/ablation profile for this year
            yr=yr+1
            
            
        ##### calculate flux and velocity
        ## ------------------------------------
        # divide/top
        ## ------------------------------------
        H_ave_div =(H[0] + H[1])/2          # average glacier thickness between divide and second point
        dHdx_div = (H[1] - H[0])/delx      # linear slope of glacier thickness between divide and second point
        dzbdx_div = (zb[1] - zb[0])/delx     # linear slope of bed between divide and second point
        dhdx_div = dHdx_div+dzbdx_div             # linear surface slope of ice between divide and second point

        Qp[0] = -(dhdx_div)**3 * H_ave_div**4*(rho*g)**3*(fd*H_ave_div+fs/H_ave_div)       # flux at plus half grid point
        up[0] = Qp[0]/H_ave_div    # velocity at half grid point

        Qm[0] = 0                      # flux at minus half grid point = 0 because we start at the divide--nothing flowing in
        um[0] = Qm[0]/H[0]       # velocity is also zero coming in

        dHdt[0] = b[0] - Qp[0]/(delx/2)    # change in thickness at the divide for this timestep is based on snowfall in/out, flux out

        ## ------------------------------------
        # toe
        ## ------------------------------------
        beyond_tu=np.nonzero(H<=0) #index of toe and beyond glacier
        beyond=beyond_tu[0] #extract indicies from tuple
        if np.size(beyond)==0:
            toe=nxs
        else:
            toe=np.min(beyond)       
    
        Qp[toe] = 0   # nothing to flux out of the end
        up[toe] = 0   # no velocity b/c no ice

        H_ave_toe = (0+H[toe-1])/2              # thickness of glacier is average of  last point with thickness and a thickness of 0     
        dHdx_toe = 0-H[toe-1]/delx              # slope is based on last point with thickness and a thickness of 0
        dzbdx_toe =(zb[toe]-zb[toe-1])/delx       # bed slope between toe and last point with thickness
        dhdx_toe = dHdx_toe+dzbdx_toe            # surface slope between toe and last point with thickness

        Qm[toe] = -(rho*g)**3*H_ave_toe**4*(dhdx_toe)**3*(fd*H_ave_toe+fs/H_ave_toe)    # flux coming from the last point with thickness and melting out before it reaches a new grid point
        um[toe] = Qm[toe]/H_ave_toe         # velocity at the toe

        dHdt[toe] = b[toe] -(Qp[toe]- Qm[toe])/delx     # change in height of the glacier at the toe based on the snowfall in/melt out and the flux in/out


        ## ------------------------------------
        # beyond glacier toe: no glacier flux
        ## ------------------------------------
        beyond=beyond[1:] #exclude first point of beyond, which is the toe
        dHdt[beyond] = b[beyond]  # thickness change is only from any change in accumulaiton in/melt out

        # no flux in, no flux out
        Qp[beyond] = 0 
        Qm[beyond] = 0

        um[beyond] =0 
        up[beyond] = 0
    
        ## ------------------------------------
        #### within the glacier
        ## ------------------------------------
        within_tup=np.nonzero(H>0) #indices of nodes where glacier exists
        within=within_tup[0]
        if np.size(within)==np.size(zb):
            within=within[0:-1]
        within=within[1:] #exclude the first, which is the divide
     
        # first for the flux going out
        H_ave_wo = (H[within+1] + H[within])/2              # the average thickness between jj and jj+1
        dHdx_wo = (H[within+1] - H[within])/delx           # the linear ice thickness slope between jj and jj+1
        dzbdx_wo = (zb[within+1] - zb[within])/delx       # the linear bed slope between jj and jj+1
        dhdx_wo = dHdx_wo + dzbdx_wo                      # the surface slope between jj and jj+1

        Qp[within] = -(rho*g)**3*H_ave_wo**4*(dhdx_wo)**3*(fd*H_ave_wo+fs/H_ave_wo)  # flux coming out
        up[within] = Qp[within]/H_ave_wo                           # velocity coming out
    
        # now for the flux coming in
        H_ave_wi = (H[within-1] + H[within])/2               # the average thickness between jj and jj-1
        dHdx_wi = (H[within] - H[within-1])/delx            # the linear ice thickness slope between jj and jj-1
        dzbdx_wi = (zb[within] - zb[within-1])/delx       # the linear bed slope between jj and jj-1
        dhdx_wi = dHdx_wi + dzbdx_wi                      # the surface slope between jj and jj-1

        Qm[within] = -(rho*g)**3*H_ave_wi**4*dhdx_wi**3*(fd*H_ave_wi+fs/H_ave_wi)        # flux coming in
        um[within] = Qm[within]/H_ave_wi            # velocity coming in
        ###### done with flux/velocity
                     
        # change in height within glacier for this timestep
        dHdt[within] = b[within] - (Qp[within] - Qm[within])/delx       # based on accumulation in/melt out and flux in/out 
        dHdt[nxs-1] = 0 # enforce no change at edge of domain

        ### done with space
        
        
        ELAind=np.max(np.nonzero(b[b>0])) #index of the ELA
        L=np.min(x[H<=0]) #glacier length
        
        glacier_ind=np.nonzero(H[H>0]) #indices of where the glacier exists
        bnet=np.sum(b[glacier_ind])/delx #net mass balance for this time step
        
            
        if not spin and t==np.floor(t): #keep the model output from this time step
            t_hold[yy]=t
            L_hold[yy]=L
            H_hold[yy,:]=H 
            Qp_hold[yy,:]=Qp #units m^2 yr^-1
            Qm_hold[yy,:]=Qm
            up_hold[yy,:]=up #units m yr^-1
            um_hold[yy,:]=um
            dHdt_hold[yy,:]=dHdt
            b_hold[yy,:]=b
            bnet_hold[yy]=bnet
            #Tyr_hold[yy]=Tyr
            #Pyr_hold[yy,:]=Pyr
            ELA_hold[yy]=H[ELAind]+zb[ELAind]
            yy=yy+1      
            
        
        # For the next timestep, the thickness is equal to the thickness profile from the previous timstep, plus the amount that the thickness has changed at each point
        # over this time period, (dHdt*delt).  The ice thickness at each point
        # needs to be positive, so we choose the max of 0 and the summation.
        # since dHdt can be negative, we want to be sure we don't have values
        # for a negative ice thickness.
        H_old=H
        H = np.maximum(0 , (H + (dHdt*delt)))
        modelyear=modelyear+delt
        
        if spin and (H_old==H).all(): #stop spin up once glacier has reached steady state
            print 'spin up done at t =',t
            break
      
    ##### end time loop
    
    if spin:
        #Hinit = H
        D['Hinit']=H
        D['x']=x
        D['zb']=zb
        D['binit']=b
        D['Qpinit']=Qp
        D['uinit']=up
        D['Linit']=L
        print "spin"
    if not spin:
        D['t']=t_hold
        D['L']=L_hold
        D['H']=H_hold
        D['Qp']=Qp_hold
        D['Qm']=Qm_hold
        D['up']=up_hold
        D['um']=um_hold
        D['dHdt']=dHdt_hold
        D['x']=x
        D['zb']=zb
        D['b']=b_hold
        D['Pyr']=Pyr_hold
        D['Tyr']=Tyr_hold
        D['ELA']=ELA_hold
        D['bnet']=bnet_hold
        print "no spin"
    
    return D
    
def filesaver(D,savedir):
    np.savetxt(os.path.join(savedir,'Qp.csv'), D['Qp'], delimiter=",")
    np.savetxt(os.path.join(savedir,'xnodes.csv'), D['x'], delimiter=",")
    np.savetxt(os.path.join(savedir,'thickness.csv'), D['H'], delimiter=",")
    np.savetxt(os.path.join(savedir,'bed.csv'), D['zb'], delimiter=",")
    np.savetxt(os.path.join(savedir,'length.csv'), D['L'], delimiter=",")
    np.savetxt(os.path.join(savedir,'time.csv'), D['t'], delimiter=",")
    np.savetxt(os.path.join(savedir,'ELA.csv'), D['ELA'], delimiter=",")
    np.savetxt(os.path.join(savedir,'bnet.csv'), D['bnet'], delimiter=",")
    np.savetxt(os.path.join(savedir,'b_spe.csv'), D['b'], delimiter=",")
    np.savetxt(os.path.join(savedir,'Hinit.csv'), D['Hinit'], delimiter=",")
    np.savetxt(os.path.join(savedir,'binit.csv'), D['binit'], delimiter=",")
    np.savetxt(os.path.join(savedir,'Qpinit.csv'), D['Qpinit'], delimiter=",")
        

if __name__ == "__main__":
    tic=time.time()
    spot = os.path.dirname(sys.argv[0])
    savedir=os.path.join(spot,'resultsFolder')
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    plt.close("all")
    bedtype = 'nisq' #choose 'nisq','lin', or 'poly'
    ######## Bed options: choose which one you like by commenting/uncommenting ######
    
    if bedtype == 'nisq':
        ##### nisqually bed, from kate allstadt's digitization of Kennard (1984?) bed from monopulse radar. Probably very wrong, but good for qualitative studies I reckon.
        DataPath = os.path.join(spot,'NisqBed.csv')
        bedprof=np.genfromtxt(DataPath,delimiter=',')
    
    elif bedtype == 'lin':
        ###### linear bed that is Nisqually-ish
        xx=np.arange(0,10100.,100)
        yy=-0.3*xx+4000
        zz=np.vstack((xx,yy))
        bedprof=zz.T
    elif bedtype == 'poly':
        ##### 2nd order poly fit to Nisqually bed
        xx=np.arange(0,10100.,100)
        yy=3.26188711e-05*xx**2 + -6.48240291e-01*xx +  4.33066591e+03
        zz=np.vstack((xx,yy))
        bedprof=zz.T
    #####
    
    
    ##### use longmire data #####
    Long_indata=np.genfromtxt('/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/DATA/Climate/Longmire_AP_MT.csv',delimiter=',')
    Long_indata=Long_indata[1:,:]
    Lprecip=Long_indata[:,1]
    Ltemp=Long_indata[:,2]
    
    datayears=Long_indata[:,0]
    modelyears=datayears[5:]
    nyrs=len(modelyears) #number of years to run simulation for. Explicit time stepping, so it takes a while if you want high-resolution space (e.g. 2 minutes for 25m grid spacing and 50-year model run on my macbook)
    
    #modelprecip=np.interp(modelyears,datayears,Lprecip)
    #modeltemp=np.interp(modelyears,datayears,Ltemp)
    
    modelprecip1=pd.Series(Lprecip).interpolate()
    modeltemp1=pd.Series(Ltemp).interpolate()
    modelprecip=np.array(modelprecip1[5:])
    modeltemp=np.array(modeltemp1[5:])    
    ##### climate data ######    
    
    ### mean temp and precip - these are values that make reasonable glacier lengths. Or, a reasonable specific mass balance.
    lapsemult=3.*2.7
    LTmean=14.0
    Tmean = LTmean + lapsemult #mean annual temperature in C (at sea level, I think)
    Pmean = 1.7 #mean annual precip (in m w.e. I think)
    #Pmean = 10*np.cos(np.pi/180*(1/45.*x-30))+1 #this is a Pmean that is smaller at the summit. 
    
    sigT = 1.0;     # standard deviation of temperature, in deg C
    sigP = 1.0;     # standard deviation of precipitation in m/yr
        
    #Tp = sigT*np.random.randn(nyrs+1); # temperature forcing (around the mean).  changes every year. Is added to mean in time-stepping loop. 
    #Pp = sigP*np.random.randn(nyrs+1);  # precip forcing (around the mean).  changes every year. 
    
    #Tp = np.zeros(nyrs+1) #use this if you just want no climate variability
    #Pp = np.zeros(nyrs+1)
    
    #Pp[11:]=1 #precip step change
    #Pp[12:14]=1 #precip delta function
    
    Pp=modelprecip-Pmean
    Tp=modeltemp-LTmean
    
    P = np.ones_like(x)*(Pmean) #mean precip
    T_wk    = (Tmean) - gam*(zb+H);    # mean temperature forcing
    melt = np.maximum(0,mu*T_wk);      # convert the temeperature to melt/ablation with meltrate factor--don't "negative" melt
    b = P-melt; # steady-state mass balance
    
    Pyr=Pmean+Pp[yr] #add noise to mean precip
    P = np.ones_like(x)*(Pyr) #add that noise over whole glacier
    Tyr=Tmean+Tp[yr] #add noise to mean temp
    T_wk    = (Tyr)*np.ones_like(x) - gam*(zb+H);    # add the noisy temp to the whole glacier and adjust for lapse rate
    melt = np.maximum(0,mu*T_wk);      # convert the temeperature to melt/ablation with meltrate factor--don't "negative" melt
    b = P-melt;     # accumulation/ablation profile for this year

    
    ######
    #nyrs=50. #number of years to run simulation for. Explicit time stepping, so it takes a while if you want high-resolution space (e.g. 2 minutes for 25m grid spacing and 50-year model run on my macbook)

    ## ------------------------------------
    #define grid spacing and time
    ## ------------------------------------
    xmx = 10000.          # max of the domain size in m
    delx = 25.              # grid spacing in m
    nxs = np.round(xmx/delx) + 1  # number of grid points
    x = np.arange(0,xmx+delx,delx)                   # x array (each point)
    
    #delt = 0.0004 # time step in yrs. Stable size is ~delx^2 / R, where R is ~3-10e5. Seems like the smaller delx is, the smaller R must be.
    dfa=10.e5 
    delt = delx**2/dfa #time step in years, varies based on grid size
    ts = 0             # starting time in yrs
    tf=nyrs-1 #number of years of model run
    ntsm=np.floor((tf-ts)/delt) + 1 # number of time steps ('round' used because need nts as integer)
    #nyrs = tf-ts       # just final time - start time

    
    ## ------------------------------------
    # glacier bed geometries
    ## ------------------------------------
    bedcoor=bedprof[0:-1,0]
    bedelev=bedprof[0:-1,1]
    zb=np.interp(x,bedcoor,bedelev) #interpolate bed onto model grid
    
    ## ------------------------------------
    # mass balance input
    ## ------------------------------------    
    elam=2800.0
    slp=116. #slope of mass balance as function of position. ~110 - 120, e.g. b=(icsu-ela)/116. 
    b_mean=(zb-elam)/slp
    sigELA=100.
    ELAnoise=sigELA*np.random.randn(nyrs+1);
    ela_var=elam+ELAnoise
    #ela_var=elam*np.zeros(np.size(ELAnoise))
    b_var=np.zeros([nyrs+1,len(x)])
    for jj in range(len(ELAnoise)):
        b_var[jj,:]=(zb-ela_var[jj])/slp #b_var dimesions: rows=time (years), columns=specific mass balance (meters) for that time step
        #b_var[10,0:100]=1.01*b_var[10,0:100]
    ########
    D={}     
    D = glacier(D,b_mean,delx,nxs,x,zb,bedprof,nyrs,spin=True) #run the model
    Dspin=D
    Hinit=Dspin['Hinit']
    binit=D['binit']
    Qpinit=D['Qpinit']
    D = glacier(D,b_var,delx,nxs,x,zb,bedprof,nyrs,spin=False,Hinit=Hinit)
    if 'H' in D:
        x=D['x'] #x nodes
        H=D['H'] #thickness at x at time t
        zb=D['zb'] #bed
        L=D['L'] #glacier length through time
        t=D['t'] # times
        Pyr=D['Pyr'] #annual precip
        Tyr=D['Tyr'] #annual T
        ELA=D['ELA'] #ELA for each t
        bnet=D['bnet'] #net annual balance for t
        b=D['b'] #specific balance
        Qp=D['Qp'] #flux
        
    #fid=os.path.join(
    Fwrite = False #true or false
    if Fwrite:
        filesaver(D,savedir)
    
    profA=np.max(np.nonzero(zb[zb>2000])) #index of glacier where bed is 2000 m
    profB=np.max(np.nonzero(zb[zb>2200])) 
    profC=np.max(np.nonzero(zb[zb>2400]))
    
    ##### Plotting #####
    plotting = True #true or false
    if plotting:
        ### tons of plots here, uncomment the ones you are interested in.
    
    ##     Glacier Geometry ###
        plt.figure(1,figsize=(20,10))
        plt.axis([0, 10000, 1000, 4500])
        plt.ion()
        plt.show()
        fig = plt.gcf()
        fig.canvas.manager.window.raise_()
    
        for jj in range(len(t)):
            plt.clf()
            plt.plot(x,H[jj,:] + zb,'c')
            plt.plot(x,zb,'k')
            plt.fill_between(x,zb,H[jj,:] + zb, color='c', alpha='0.5')
            plt.grid(True)
            plt.title('Glacier Geometry')
            plt.draw()
            plt.pause(0.001)
        ######
    
        #####Initial Thickness
        plt.figure(2)
        plt.plot(x,Hinit,'c*')
        plt.title('Initial Thickness')
        #####Initial Geometry
        plt.figure()
        plt.plot(x,Hinit + zb,'c')
        plt.plot(x,zb,'k')
        plt.title('Initial Geometry')
        #####Length through time
        plt.figure()
        plt.plot(t,L)
        plt.title('Length')
        plt.show()
        ######
    
    #    ## Thickness Deviation #####     
    #    plt.figure(3,figsize=(20,10))
    #    plt.axis([0, 10000, -100, 100])
    #    plt.ion()
    #    plt.show()
    #    fig = plt.gcf()
    #    fig.canvas.manager.window.raise_()
    #
    #    for kk in range(len(t)):
    #        tstr='time = %s years' % t[kk]
    #        plt.clf()
    #        plt.plot(x,H[kk,:] - Hinit,'c')
    #        plt.axis([0, 10000, -30, 30])
    #        plt.plot(x,zb,'k')
    #        plt.fill_between(x,zb,H[kk,:] + zb, color='c', alpha='0.5')
    #        plt.grid(True)
    #        plt.title('Deviation from Hinit')
    #        plt.text(8000,20,tstr)
    #        plt.draw()
    #        plt.pause(0.7)        
    #    ######
    
    
    #    ## Thickness Deviation #####     
    #    plt.figure(4,figsize=(20,10))
    #    plt.axis([0, 10000, -100, 100])
    #    plt.ion()
    #    plt.show()
    #    fig = plt.gcf()
    #    fig.canvas.manager.window.raise_()
    #
    #    for kk in range(len(t)):
    #        tstr='time = %s years' % t[kk]
    #        plt.clf()
    #        plt.plot(x,(H[kk,:] - Hinit)/Hinit,'c')
    #        plt.axis([0, 10000, -1, 1])
    #        plt.plot(x,zb,'k')
    #        plt.fill_between(x,zb,H[kk,:] + zb, color='c', alpha='0.5')
    #        plt.grid(True)
    #        plt.title('Deviation from Hinit')
    #        plt.text(8000,20,tstr)
    #        plt.draw()
    #        plt.pause(0.7)        
    #    ######
    
    #    ### specific mass balance through time  #####     
    #    plt.figure(5,figsize=(20,10))
    #    plt.axis([0, 10000, -100, 100])
    #    plt.ion()
    #    plt.show()
    #    fig = plt.gcf()
    #    fig.canvas.manager.window.raise_()
    #
    #    for kk in range(len(t)):
    #        tstr='time = %s years' % t[kk]
    #        plt.clf()
    #        plt.plot(x,b[kk,:]-binit,'b')
    #        plt.axis([0, 10000, -0.5, 0.5])
    #        plt.grid(True)
    #        plt.title('Mass Balance')
    #        plt.text(8000,0,tstr)
    #        plt.draw()
    #        plt.pause(0.3)        
    #    #######
    
        ## Thickness at profiles, length, and ELA
        plt.figure(6,figsize=(20,10))
    
        plt.subplot(5,1,1)
        plt.title('Profile elevation,Length,ELA')
        plt.plot(t,H[:,profC]+zb[profC])
        plt.grid(True)
    
        plt.subplot(5,1,2)
        plt.plot(t,H[:,profB]+zb[profB])
        plt.grid(True)
    
        plt.subplot(5,1,3)
        plt.plot(t,H[:,profA]+zb[profA])
        plt.grid(True)
    
        plt.subplot(5,1,4)
        plt.plot(t,L)
        plt.grid(True)
    
        plt.subplot(5,1,5)
        plt.plot(t,ELA)
        plt.grid(True)
        #plt.savefig(savedir + '/SE_' + bedtype +'.png', bbox_inches='tight')
            
        plt.show()
        fig = plt.gcf()
        fig.canvas.manager.window.raise_()
        ######
    
#        ## Thickness at profiles, length, and net b_dot through time
#        plt.figure(7,figsize=(20,10))
#    
#
#        plt.subplot(5,1,1)
#        plt.title('Surface Elevation,Length,bnet')
#        plt.plot(t,H[:,profC]+zb[profC])
#        plt.grid(True)
#    
#        plt.subplot(5,1,2)
#        plt.plot(t,H[:,profB]+zb[profB])
#        plt.grid(True)
#        
#        plt.subplot(5,1,3)
#        plt.plot(t,H[:,profA]+zb[profA])
#        plt.grid(True)
#
#        plt.subplot(5,1,4)
#        plt.plot(t,L)
#        plt.grid(True)
#    
#        plt.subplot(5,1,5)
#        plt.plot(t,Tyr-Tmean,'b')
#        #plt.plot(t,Pyr-Pmean[profA],'g')
#        plt.plot(t,Pyr-Pmean,'g')
#        plt.grid(True)
#        plt.subplot(5,1,5)
#        plt.plot(t,bnet,'b')
#        plt.grid(True)
#        
#
#        plt.show()
#        fig = plt.gcf()
#        fig.canvas.manager.window.raise_()
#        ######

#        #### Thickness deviation from initial at profiles, length, and net b_dot through time
#        plt.figure(8,figsize=(20,10))
#    
#
#        plt.subplot(5,1,1)
#        plt.title('Profile deviation,Length,bnet')
#        plt.plot(t,H[:,profC]-Hinit[profC])
#        plt.grid(True)
#    
#        plt.subplot(5,1,2)
#        plt.plot(t,H[:,profB]-Hinit[profB])
#        plt.grid(True)
#
#        plt.subplot(5,1,3)
#        plt.plot(t,H[:,profA]-Hinit[profA])
#        plt.grid(True)
#
#        plt.subplot(5,1,4)
#        plt.plot(t,L)
#        plt.grid(True)
#
#        plt.subplot(5,1,5)
#        plt.plot(t,Tyr-Tmean,'b')
#        plt.plot(t,Pyr-Pmean[profA],'g')
#        plt.grid(True)
#        plt.subplot(5,1,5)
#        plt.plot(t,bnet,'b')
#        plt.grid(True)
#
#        plt.show()
#        fig = plt.gcf()
#        fig.canvas.manager.window.raise_()
#        #####
     
        ##### flux
        plt.figure(9,figsize=(20,10))

        plt.title('Flux')
        plt.plot(t,H[:,profC]+zb[profC])
        plt.plot(x,Qp[11,:]-Qpinit,'k')
        plt.plot(x,Qp[12,:]-Qpinit,'b')
        plt.plot(x,Qp[13,:]-Qpinit,'r')
        plt.plot(x,Qp[14,:]-Qpinit,'g')
        plt.plot(x,Qp[15,:]-Qpinit,'m')
        plt.plot(x,Qp[16,:]-Qpinit,'k--')
        plt.plot(x,Qp[17,:]-Qpinit,'b--')
        plt.plot(x,Qp[18,:]-Qpinit,'r--')
        plt.plot(x,Qp[19,:]-Qpinit,'g--')
        plt.grid(True)
        #plt.savefig(savedir + '/flux_' + bedtype + '.eps', bbox_inches='tight')
        plt.show()
        fig = plt.gcf()
        fig.canvas.manager.window.raise_()
        
    ######  end plots
            
    
    #Pnorm=Pyr-Pmean
    #Tnorm=Tyr-Tmean
    #plt.figure()
    #plt.plot(t,Pnorm,'b')
    #plt.plot(t,Tnorm,'g')
   

    elapsed=time.time()-tic
    elapsed_min=elapsed/60.
    mins=np.floor(elapsed_min)
    secs=(elapsed_min-mins)*60 

    print "run took %s seconds" % elapsed


