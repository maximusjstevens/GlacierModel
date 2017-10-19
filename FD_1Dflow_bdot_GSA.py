### 1-D Glacier Flow Model based on continuity: dH/dt = dQ/dx + b_dot
### Finite Difference, explicit time steps
### Python version coded by Max Stevens, 8/14/15, based on (copied from) Matlab model from Gerard Roe/Kat Huybers/Summer Rupper
### References: Oerlemans 2001: Glaciers and climate change. Page 61. 

### this version is to use mass balance perturbations.

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time
import matplotlib
import os
import sys
matplotlib.rcParams['axes.formatter.useoffset']=False


class Glacier:

	def __init__(self, spot, D, spin):

		## -------------------------------------
		# define parameters and constants
		## -------------------------------------
		self.rho=917.0              # density of ice in kg/m^3
		self.g=9.8                  # gravity in m/s^2
		self.n=3.0                    # glen's flow law constant
		self.mu = 0.6      # melt rate in m /yr /degC
		self.gam = 6.0e-3    # atmospheric lapse rate in degC/m

		self.s_per_year = 365.25*24*3600  # number of seconds in a year (s yr^-1)
		self.A_T = 6.8e-24*self.s_per_year   # softness parameter for -5 degree ice (yr^-1 Pa^-3). Units given in Oerlemans are per second, so need sPerYear conversion
		self.fd = 2*self.A_T/(self.n+2)  # constants lumped together, and put into standard units.  (Pa^-3 yr^-1)
		self.fs = 1.0e-20*self.s_per_year  # sliding parameter fromOerlemans, page 62 (Pa^-3 m^2 yr-1). Units given in Oerlemans are per second

		## ------------------------------------
		#define grid spacing 
		## ------------------------------------
		self.xmx = 10000.          # max of the domain size in m
		self.delx = 20.              # grid spacing in m
		self.nxs = np.round(self.xmx/self.delx) + 1  # number of grid points
		self.x = np.arange(0,self.xmx+self.delx,self.delx)                   # x array (each point)

		bedtype = 'nisq' #choose 'nisq','lin', or 'poly'
	
		if bedtype == 'nisq':
			##### nisqually bed, from kate allstadt's digitization of Kennard (1984?) bed from monopulse radar. Probably very wrong, but good for qualitative studies I reckon.
			DataPath = os.path.join(spot,'NisqBed.csv')
			self.bedprof=np.genfromtxt(DataPath,delimiter=',')
	
		elif bedtype == 'lin':
			###### linear bed that is Nisqually-ish
			xx=np.arange(0,10100.,100)
			yy=-0.3*xx+4000
			zz=np.vstack((xx,yy))
			self.bedprof=zz.T

		elif bedtype == 'poly':
			##### 2nd order poly fit to Nisqually bed
			xx=np.arange(0,10100.,100)
			yy=3.26188711e-05*xx**2 + -6.48240291e-01*xx +  4.33066591e+03
			zz=np.vstack((xx,yy))
			self.bedprof=zz.T

		## ------------------------------------
		# glacier bed geometries
		## ------------------------------------
		bedcoor=self.bedprof[0:-1,0]
		bedelev=self.bedprof[0:-1,1]
		self.zb=np.interp(self.x,bedcoor,bedelev) #interpolate bed onto model grid

		## ------------------------------------
		# load mass balance
		## ------------------------------------		
		self.mbr_s=np.genfromtxt(os.path.join(spot,'Nisq_MB_Bs_elevation.csv'),delimiter=',')[:,1:] #mass balance record from Al
		# mbr_w=np.genfromtxt('Nisq_MB_1917_2015.csv',delimiter=',') #mass balance record from Al
		self.mbr=np.genfromtxt(os.path.join(spot,'Nisq_MB_1917_2015.csv'),delimiter=',') #mass balance record from Al
		self.hyp = np.genfromtxt(os.path.join(spot,'Nisq_hyp.txt'))

		## ------------------------------------
		# time 
		## ------------------------------------
		self.mbr_years=self.mbr[0,:]
		self.nyrs=self.mbr_years[-1]-self.mbr_years[0] #number of years to run simulation for. Explicit time stepping, so it takes a while if you want high-resolution space (e.g. 2 minutes for 25m grid spacing and 50-year model run on my macbook)

		dfa=10.e5 
		self.delt = self.delx**2/dfa #time step in years, varies based on grid size

		if spin:
			self.ts = 0.0
			self.tf = 100.0

		else:
			self.ts = self.mbr_years[0]
			self.tf = self.mbr_years[-1] #number of years of model run

		# nts=np.floor((tf-ts)/delt) + 1 # number of time steps ('round' used because need nts as integer)
		self.ntsm=np.floor((self.tf-self.ts)/self.delt) + 1 # number of time steps ('round' used because need nts as integer)
		# self.nyrs = self.tf-self.ts       # just final time - start time
		self.tt = np.arange(self.ts,self.tf)

		ba=self.mbr[3,:]
		# ba_n=(ba-np.mean(ba))/(np.std(ba))   
		elam=2500.0 #mean ELA
		slp=116. #slope of mass balance as function of position. ~110 - 120, e.g. b=(icsu-ela)/116. 
		if spin:
			self.b_in=(self.zb-elam)/slp #mean specific mass balance profile
		# sigELA=100.
		# ELAnoise=sigELA*np.random.randn(nyrs+1);
		# ELAnoise=ba*elam/10
		ela_var=-300.0*ba+elam
		# ela_var=elam-ELAnoise #variable annual ELA
		# b_var = (zb-ELAnoise)/slp
		#ela_var=elam*np.zeros(np.size(ELAnoise))
		b_var=np.zeros([int(self.nyrs+1),len(self.x)])
		# rr=np.linspace(0,100,len(ELAnoise))
		for jj in range(len(ela_var)):
			b_var[jj,:]=(self.zb-ela_var[jj])/slp #b_var dimesions: rows=time (years), columns=specific mass balance (meters) for that time step
		#     #b_var[10,0:100]=1.01*b_var[10,0:100]
		#     b_var[jj,:]=b_mean
		#     b_var[jj,:]=(self.zb-elam-rr[jj])/slp
		if not spin:
			self.b_in = b_var
		# b_var[10,:]=(self.zb-elam+200)/slp # + is ela lowering, - is ela higher
		#b_var[10,:]=(self.zb-elam-100)/slp
		#b_var[10,0:126]=b_var[10,0:126]*1.1    
		########
		# D={}     
		# D = glacier(D,b_mean,self.delx,nxs,x,self.zb,bedprof,nyrs,spin=True) #run the model

		# D = glacier(D,b_var,self.delx,nxs,x,self.zb,bedprof,nyrs,spin=False,Hinit=Hinit)



		

	def glacier(self,Hinit=None):


		# nts=np.floor((tf-ts)/delt) + 1 # number of time steps ('round' used because need nts as integer)
		#nyrs = tf-ts       # just final time - start time

		## ------------------------------------
		# initialize arrays
		## ------------------------------------

		if Hinit is None:
			H0 = 20.0*np.ones_like(self.x)     # initial thickness of ice (in m)--do not confuse with h, ice elevation.  h = H+zb
			H0[-1]=0.0                 # impose 0 thickness at end of domain.
		else:
			H0 = Hinit

		Qp=np.zeros_like(self.x)      # Qp equals jj+1/2 flux
		Qp_slide=np.zeros_like(self.x)
		Qp_def=np.zeros_like(self.x)
		Qm=np.zeros_like(self.x)      # Qm equals jj-1/2 flux
		up = np.zeros_like(self.x)    # velocity at jj +1/2 grid 
		um = np.zeros_like(self.x)    # velcotiy at jj-1/2 grid
		dHdt=np.zeros_like(self.x)        #  slope of ice thickness -- do not confuse with dhdx, surface slope.  dhdx = dzdx+dHdx
		H = H0                    # initialize height array

		yr = 0         # for counting years
		idx_out = 0    # 
		deltout = 1    # for counting years
		nouts = np.round(self.ntsm/1)   # for counting years
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
			t_hold=np.zeros((int(nyrs)))
			L_hold=np.zeros_like(t_hold)
			H_hold=np.zeros((int(nyrs),int(nxs)))
			Qp_hold=np.zeros_like(H_hold)
			Qp_slide_hold=np.zeros_like(H_hold)
			Qp_def_hold=np.zeros_like(H_hold)
			Qm_hold=np.zeros_like(H_hold)
			up_hold=np.zeros_like(H_hold)
			um_hold=np.zeros_like(H_hold)
			dHdt_hold=np.zeros_like(H_hold)
			b_hold=np.zeros_like(H_hold)
			ELA_hold=np.zeros_like(t_hold)
			bnet_hold=np.zeros_like(t_hold)
			Tyr_hold=np.zeros_like(H_hold) 
			Pyr_hold=np.zeros_like(H_hold)  

		## ------------------------------------
		# begin loop over time
		## ------------------------------------
		yy=0 #for tracking
		
		Pmean = 10.0*np.cos(np.pi/180.0*(1/45.0*self.x-30.0))+1.0
		modelyear=1900.0
		for ii in range(int(self.ntsm)):
			seasontime=modelyear-np.floor(modelyear)
	  
			t = self.delt*(ii) # time in years

			if spin:
				#ela=3000.0
				#icsu=H+self.zb
				#b=(icsu-ela)/116.
				b=self.b_in
				
			if not spin and t==np.floor(t):
				#b = 0.004*self.zb-11.236
				b=self.b_in[yr,:]
				#ela=3000.0
				#icsu=H+self.zb
				#b=(icsu-ela)/116.
				#if seasontime>0.5 and seasontime<0.83:

				#Pyr=Pmean+Pp[yr] #add noise to mean precip
				#P = np.ones_like(self.x)*(Pyr) #add that noise over whole glacier
				#Tyr=Tmean+Tp[yr] #add noise to mean temp
				#T_wk    = (Tyr)*np.ones_like(self.x) - gam*(self.zb+H);    # add the noisy temp to the whole glacier and adjust for lapse rate
				#melt = np.maximum(0,mu*T_wk);      # convert the temeperature to melt/ablation with meltrate factor--don't "negative" melt
				#b = P-melt;     # accumulation/ablation profile for this year
				yr=yr+1
				
				
			##### calculate flux and velocity
			## ------------------------------------
			# divide/top
			## ------------------------------------
			H_ave_div =(H[0] + H[1])/2          # average glacier thickness between divide and second point
			dHdx_div = (H[1] - H[0])/self.delx      # linear slope of glacier thickness between divide and second point
			dzbdx_div = (self.zb[1] - self.zb[0])/self.delx     # linear slope of bed between divide and second point
			dhdx_div = dHdx_div+dzbdx_div             # linear surface slope of ice between divide and second point

			Qp[0] = -(dhdx_div)**3 * H_ave_div**4*(self.rho*self.g)**3*(self.fd*H_ave_div+self.fs/H_ave_div)       # flux at plus half grid point
			up[0] = Qp[0]/H_ave_div    # velocity at half grid point

			Qp_slide[0]=self.fs/H_ave_div
			Qp_def[0]=self.fd*H_ave_div

			Qm[0] = 0                      # flux at minus half grid point = 0 because we start at the divide--nothing flowing in
			um[0] = Qm[0]/H[0]       # velocity is also zero coming in

			dHdt[0] = b[0] - Qp[0]/(self.delx/2)    # change in thickness at the divide for this timestep is based on snowfall in/out, flux out

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

			Qp_slide[toe]=0
			Qp_def[toe]=0

			H_ave_toe = (0+H[toe-1])/2              # thickness of glacier is average of  last point with thickness and a thickness of 0     
			dHdx_toe = 0-H[toe-1]/self.delx              # slope is based on last point with thickness and a thickness of 0
			dzbdx_toe =(self.zb[toe]-self.zb[toe-1])/self.delx       # bed slope between toe and last point with thickness
			dhdx_toe = dHdx_toe+dzbdx_toe            # surface slope between toe and last point with thickness

			Qm[toe] = -(self.rho*self.g)**3*H_ave_toe**4*(dhdx_toe)**3*(self.fd*H_ave_toe+self.fs/H_ave_toe)    # flux coming from the last point with thickness and melting out before it reaches a new grid point
			um[toe] = Qm[toe]/H_ave_toe         # velocity at the toe

			dHdt[toe] = b[toe] -(Qp[toe]- Qm[toe])/self.delx     # change in height of the glacier at the toe based on the snowfall in/melt out and the flux in/out


			## ------------------------------------
			# beyond glacier toe: no glacier flux
			## ------------------------------------
			beyond=beyond[1:] #exclude first point of beyond, which is the toe
			dHdt[beyond] = b[beyond]  # thickness change is only from any change in accumulaiton in/melt out

			# no flux in, no flux out
			Qp[beyond] = 0
			Qp_slide[beyond]=0
			Qp_def[beyond]=0 
			Qm[beyond] = 0

			um[beyond] =0 
			up[beyond] = 0
		
			## ------------------------------------
			#### within the glacier
			## ------------------------------------
			within_tup=np.nonzero(H>0) #indices of nodes where glacier exists
			within=within_tup[0]
			if np.size(within)==np.size(self.zb):
				within=within[0:-1]
			within=within[1:] #exclude the first, which is the divide
		 
			# first for the flux going out
			H_ave_wo = (H[within+1] + H[within])/2              # the average thickness between jj and jj+1
			dHdx_wo = (H[within+1] - H[within])/self.delx           # the linear ice thickness slope between jj and jj+1
			dzbdx_wo = (self.zb[within+1] - self.zb[within])/self.delx       # the linear bed slope between jj and jj+1
			dhdx_wo = dHdx_wo + dzbdx_wo                      # the surface slope between jj and jj+1

			Qp[within] = -(self.rho*self.g)**3*H_ave_wo**4*(dhdx_wo)**3*(self.fd*H_ave_wo+self.fs/H_ave_wo)  # flux coming out
			up[within] = Qp[within]/H_ave_wo                           # velocity coming out

			Qp_slide[within]=self.fs/H_ave_wo
			Qp_def[within]=self.fd*H_ave_wo        
			# now for the flux coming in
			H_ave_wi = (H[within-1] + H[within])/2               # the average thickness between jj and jj-1
			dHdx_wi = (H[within] - H[within-1])/self.delx            # the linear ice thickness slope between jj and jj-1
			dzbdx_wi = (self.zb[within] - self.zb[within-1])/self.delx       # the linear bed slope between jj and jj-1
			dhdx_wi = dHdx_wi + dzbdx_wi                      # the surface slope between jj and jj-1

			Qm[within] = -(self.rho*self.g)**3*H_ave_wi**4*dhdx_wi**3*(self.fd*H_ave_wi+self.fs/H_ave_wi)        # flux coming in
			um[within] = Qm[within]/H_ave_wi            # velocity coming in
			###### done with flux/velocity
						 
			# change in height within glacier for this timestep
			dHdt[within] = b[within] - (Qp[within] - Qm[within])/self.delx       # based on accumulation in/melt out and flux in/out 
			dHdt[int(self.nxs-1)] = 0 # enforce no change at edge of domain

			### done with space
			
			
			ELAind=np.max(np.nonzero(b[b>0])) #index of the ELA
			L=np.min(self.x[H<=0]) #glacier length
			
			glacier_ind=np.nonzero(H[H>0]) #indices of where the glacier exists
			bnet=np.sum(b[glacier_ind])/self.delx #net mass balance for this time step
			
				
			if not spin and t==np.floor(t): #keep the model output from this time step
				t_hold[yy]=t
				L_hold[yy]=L
				H_hold[yy,:]=H 
				Qp_hold[yy,:]=Qp #units m^2 yr^-1
				Qp_slide_hold[yy,:]=Qp_slide
				Qp_def_hold[yy,:]=Qp_def
				Qm_hold[yy,:]=Qm
				up_hold[yy,:]=up #units m yr^-1
				um_hold[yy,:]=um
				dHdt_hold[yy,:]=dHdt
				b_hold[yy,:]=b
				bnet_hold[yy]=bnet
				#Tyr_hold[yy,:]=T_wk
				#Pyr_hold[yy,:]=Pyr
				ELA_hold[yy]=H[ELAind]+self.zb[ELAind]
				yy=yy+1  
				
			
			# For the next timestep, the thickness is equal to the thickness profile from the previous timstep, plus the amount that the thickness has changed at each point
			# over this time period, (dHdt*delt).  The ice thickness at each point
			# needs to be positive, so we choose the max of 0 and the summation.
			# since dHdt can be negative, we want to be sure we don't have values
			# for a negative ice thickness.
			H_old=H
			H = np.maximum(0 , (H + (dHdt*self.delt)))
			modelyear=modelyear+self.delt
			
			if spin and (H_old==H).all(): #stop spin up once glacier has reached steady state
				print('spin up done at t =',t)
				break
		  
		##### end time loop
		
		if spin:
			#Hinit = H
			D['Hinit']=H
			D['x']=self.x
			D['zb']=self.zb
			D['binit']=b
			D['Qpinit']=Qp
			D['uinit']=up
			D['Linit']=L
			print("spin")
		if not spin:
			D['t']=t_hold
			D['L']=L_hold
			D['H']=H_hold
			D['Qp']=Qp_hold
			D['Qp_slide']=Qp_slide_hold
			D['Qp_def']=Qp_def_hold
			D['Qm']=Qm_hold
			D['up']=up_hold
			D['um']=um_hold
			D['dHdt']=dHdt_hold
			D['x']=self.x
			D['zb']=self.zb
			D['b']=b_hold
			#D['Pyr']=Pyr_hold
			#D['Tyr']=Tyr_hold
			D['ELA']=ELA_hold
			D['bnet']=bnet_hold
			print("no spin")
		
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
	
	spot = os.path.dirname(sys.argv[0])
	savedir=os.path.join(spot,'resultsFolder/GSA')
	if not os.path.exists(savedir):
		os.makedirs(savedir)

	plt.close("all")
	D={}
	spin = True
	gS = Glacier(spot, D, spin)
	D = gS.glacier()

	Dspin=D
	Hinit=Dspin['Hinit']
	binit=D['binit']
	Qpinit=D['Qpinit']
	print("init length=", D['Linit'])

	spin = False
	gNS = Glacier(spot, D, spin)
	D = gNS.glacier(Hinit)







	

	
	if 'H' in D:
		x=D['x'] #x nodes
		H=D['H'] #thickness at x at time t
		zb=D['zb'] #bed
		L=D['L'] #glacier length through time
		t=D['t'] # times
		#Pyr=D['Pyr'] #annual precip
		#Tyr=D['Tyr'] #annual T
		ELA=D['ELA'] #ELA for each t
		bnet=D['bnet'] #net annual balance for t
		b=D['b'] #specific balance
		Qp=D['Qp'] #flux
		
	#fid=os.path.join(
	Fwrite = True #true or false
	if Fwrite:
		filesaver(D,savedir)
	
	profA=np.max(np.nonzero(zb[zb>1550])) #index of glacier where bed is 2000 m
	profB=np.max(np.nonzero(zb[zb>1750])) 
	profC=np.max(np.nonzero(zb[zb>2000]))

	D['profA']=profA
	D['profB']=profB
	D['profC']=profC
	D['tt'] = tt
	
	##### Plotting #####
	plotting = False #true or false
	if plotting:
		### tons of plots here, uncomment the ones you are interested in.
	
		##Glacier Geometry ###
		plt.figure(1,figsize=(20,10))
		plt.axis([0, 10000, 1000, 4500])
		plt.ion()
		plt.show()
		fig = plt.gcf()
		# fig.canvas.manager.window.raise_()
	
		for jj in range(len(t)):
			plt.clf()
			plt.plot(x,H[jj,:] + zb,'c')
			plt.plot(x,zb,'k')
			plt.fill_between(x,zb,H[jj,:] + zb)#, color='c', alpha='0.5')
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
	
	   # ## Thickness Deviation #####     
	   # plt.figure(3,figsize=(20,10))
	   # plt.axis([0, 10000, -100, 100])
	   # plt.ion()
	   # plt.show()
	   # fig = plt.gcf()
	   # fig.canvas.manager.window.raise_()
	
	   # for kk in range(len(t)):
	   #     tstr='time = %s years' % t[kk]
	   #     plt.clf()
	   #     plt.plot(x,H[kk,:] - Hinit,'c')
	   #     plt.axis([0, 10000, -30, 30])
	   #     plt.plot(x,zb,'k')
	   #     plt.fill_between(x,zb,H[kk,:] + zb, color='c', alpha='0.5')
	   #     plt.grid(True)
	   #     plt.title('Deviation from Hinit')
	   #     plt.text(8000,20,tstr)
	   #     plt.draw()
	   #     plt.pause(0.7)        
	   # ######
	
	
		## Thickness Deviation #####     
		plt.figure(4,figsize=(20,10))
		plt.axis([0, 10000, -100, 100])
		plt.ion()
		plt.show()
		fig = plt.gcf()
		# fig.canvas.manager.window.raise_()

		for kk in range(len(t)):
		   tstr='time = %s years' % t[kk]
		   plt.clf()
		   plt.plot(x,(H[kk,:] - Hinit)/Hinit,'c')
		   plt.axis([0, 10000, -1, 1])
		   plt.plot(x,zb,'k')
		   plt.fill_between(x,zb,H[kk,:] + zb)#, color='c', alpha='0.5')
		   plt.grid(True)
		   plt.title('Deviation from Hinit')
		   plt.text(8000,20,tstr)
		   plt.draw()
		   plt.pause(0.7)        
		######
	
	   # ### specific mass balance through time  #####     
	   # plt.figure(5,figsize=(20,10))
	   # plt.axis([0, 10000, -100, 100])
	   # plt.ion()
	   # plt.show()
	   # fig = plt.gcf()
	   # fig.canvas.manager.window.raise_()
	
	   # for kk in range(len(t)):
	   #     tstr='time = %s years' % t[kk]
	   #     plt.clf()
	   #     plt.plot(x,b[kk,:]-binit,'b')
	   #     plt.axis([0, 10000, -0.5, 0.5])
	   #     plt.grid(True)
	   #     plt.title('Mass Balance')
	   #     plt.text(8000,0,tstr)
	   #     plt.draw()
	   #     plt.pause(0.3)        
	   # #######
	
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
		# fig.canvas.manager.window.raise_()
		######
	
		## Thickness at profiles, length, and net b_dot through time
		plt.figure(7,figsize=(20,10))
	

		plt.subplot(5,1,1)
		plt.title('Surface Elevation,Length,bnet')
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
		plt.plot(t,bnet,'b')
		#plt.plot(t,Pyr-Pmean[profA],'g')
		#plt.plot(t,Pyr-Pmean,'g')
		plt.grid(True)
		plt.subplot(5,1,5)
		plt.plot(t,bnet,'b')
		plt.grid(True)
		

		plt.show()
		fig = plt.gcf()
		# fig.canvas.manager.window.raise_()
		######

	   # #### Thickness deviation from initial at profiles, length, and net b_dot through time
	   # plt.figure(8,figsize=(20,10))
   

	   # plt.subplot(5,1,1)
	   # plt.title('Profile deviation,Length,bnet')
	   # plt.plot(t,H[:,profC]-Hinit[profC])
	   # plt.grid(True)
   
	   # plt.subplot(5,1,2)
	   # plt.plot(t,H[:,profB]-Hinit[profB])
	   # plt.grid(True)

	   # plt.subplot(5,1,3)
	   # plt.plot(t,H[:,profA]-Hinit[profA])
	   # plt.grid(True)

	   # plt.subplot(5,1,4)
	   # plt.plot(t,L)
	   # plt.grid(True)

	   # plt.subplot(5,1,5)
	   # plt.plot(t,Tyr-Tmean,'b')
	   # plt.plot(t,Pyr-Pmean[profA],'g')
	   # plt.grid(True)
	   # plt.subplot(5,1,5)
	   # plt.plot(t,bnet,'b')
	   # plt.grid(True)

	   # plt.show()
	   # fig = plt.gcf()
	   # fig.canvas.manager.window.raise_()
	   # #####
	 
		##### flux
		plt.figure(9,figsize=(20,10))

		plt.title('Flux')
		plt.plot(t,H[:,profC]+zb[profC])
		#plt.plot(x,Qp[11,:]-Qpinit,'k',label='11')
		#plt.plot(x,Qp[12,:]-Qpinit,'b',label='12')
		#plt.plot(x,Qp[13,:]-Qpinit,'r',label='13')
		#plt.plot(x,Qp[14,:]-Qpinit,'g',label='14')
		#plt.plot(x,Qp[15,:]-Qpinit,'m',label='15')
		#plt.plot(x,Qp[16,:]-Qpinit,'k--',label='16')
		#plt.plot(x,Qp[17,:]-Qpinit,'b--',label='17')
		#plt.plot(x,Qp[18,:]-Qpinit,'r--',label='18')
		#plt.plot(x,Qp[19,:]-Qpinit,'g--',label='19')
		
		plt.plot(x,Qp[11,:]-Qp[10,:],'k',label='11')
		plt.plot(x,Qp[12,:]-Qp[11,:],'b',label='12')
		plt.plot(x,Qp[13,:]-Qp[12,:],'r',label='13')
		plt.plot(x,Qp[14,:]-Qp[13,:],'g',label='14')
		plt.plot(x,Qp[15,:]-Qp[14,:],'m',label='15')
		plt.plot(x,Qp[16,:]-Qp[15,:],'k--',label='16')
		plt.plot(x,Qp[17,:]-Qp[16,:],'b--',label='17')
		plt.plot(x,Qp[18,:]-Qp[17,:],'r--',label='18')
		plt.plot(x,Qp[19,:]-Qp[18,:],'g--',label='19')        
		
		plt.grid(True)
		plt.legend()
		#plt.savefig(savedir + '/flux_' + bedtype + '.eps', bbox_inches='tight')
		plt.show()
		fig = plt.gcf()
		# fig.canvas.manager.window.raise_()
		

	######  end plots
			
	
	#Pnorm=Pyr-Pmean
	#Tnorm=Tyr-Tmean
	#plt.figure()
	#plt.plot(t,Pnorm,'b')
	#plt.plot(t,Tnorm,'g')
	np.save('gmodel.npy',D)

	elapsed=time.time()-tic
	elapsed_min=elapsed/60.
	mins=np.floor(elapsed_min)
	secs=(elapsed_min-mins)*60 

	print("run took %s seconds" % elapsed)


