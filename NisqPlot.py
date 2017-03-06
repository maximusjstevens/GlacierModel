# NisqPlot.py
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time
import matplotlib

D2=np.load('gmodel.npy')
D=D2[()]

x=D['x'] #x nodes
H=D['H'] #thickness at x at time t
zb=D['zb'] #bed
L=D['L'] #glacier length through time
t=D['tt'] # times
#Pyr=D['Pyr'] #annual precip
#Tyr=D['Tyr'] #annual T
ELA=D['ELA'] #ELA for each t
bnet=D['bnet'] #net annual balance for t
b=D['b'] #specific balance
Qp=D['Qp'] #flux
Hinit=D['Hinit']
binit=D['binit']
Qpinit=D['Qpinit']

# profA = D['profA']
# profB = D['profB']
# profC = D['profC']

profA=np.max(np.nonzero(zb[zb>1600])) #index of glacier where bed is 2000 m
profB=np.max(np.nonzero(zb[zb>1800])) 
profC=np.max(np.nonzero(zb[zb>2000]))

# ##     Glacier Geometry ###
# plt.figure(1,figsize=(20,10))
# plt.axis([0, 10000, 1000, 4500])
# plt.ion()
# plt.show()
# fig = plt.gcf()
# # fig.canvas.manager.window.raise_()

# for jj in range(len(t)):
#     plt.clf()
#     plt.plot(x,H[jj,:] + zb,'c')
#     plt.plot(x,zb,'k')
#     plt.fill_between(x,zb,H[jj,:] + zb)#, color='c', alpha='0.5')
#     plt.grid(True)
#     plt.title('Glacier Geometry')
#     plt.draw()
#     plt.pause(0.001)
# ######

#####Initial Thickness
plt.figure(1)
plt.plot(x,H[-1,:],'c*')
plt.title('Final Thickness')
# #####Initial Geometry
# plt.figure()
# plt.plot(x,Hinit + zb,'c')
# plt.plot(x,zb,'k')
# plt.title('Initial Geometry')
# #####Length through time
# plt.figure()
# plt.plot(t,L)
# plt.title('Length')
# plt.show()
# ######


#####Initial Thickness
plt.figure(2)
plt.plot(x,Hinit,'c*')
plt.title('Initial Thickness')
#####Initial Geometry
plt.figure(51)
plt.plot(x,Hinit + zb,'c')
plt.plot(x,zb,'k')
plt.title('Initial Geometry')
#####Length through time
plt.figure(52)
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
# # fig.canvas.manager.window.raise_()

# for kk in range(len(t)):
#    tstr='time = %s years' % t[kk]
#    plt.clf()
#    plt.plot(x,H[kk,:] - Hinit,'c')
#    plt.axis([0, 10000, -30, 30])
#    plt.plot(x,zb,'k')
#    plt.fill_between(x,zb,H[kk,:] + zb)#, color='c', alpha='0.5')
#    plt.grid(True)
#    plt.title('Deviation from Hinit')
#    plt.text(8000,20,tstr)
#    plt.draw()
#    plt.pause(0.7)        
# ######


# ## Thickness Deviation #####     
# plt.figure(4,figsize=(20,10))
# plt.axis([0, 10000, -100, 100])
# plt.ion()
# plt.show()
# fig = plt.gcf()
# # fig.canvas.manager.window.raise_()

# for kk in range(len(t)):
#    tstr='time = %s years' % t[kk]
#    plt.clf()
#    plt.plot(x,(H[kk,:] - Hinit)/Hinit,'c')
#    plt.axis([0, 10000, -1, 1])
#    plt.plot(x,zb,'k')
#    plt.fill_between(x,zb,H[kk,:] + zb)#, color='c', alpha='0.5')
#    plt.grid(True)
#    plt.title('Deviation from Hinit')
#    plt.text(8000,20,tstr)
#    plt.draw()
#    plt.pause(0.7)        
# ######

# ### specific mass balance through time  #####     
# plt.figure(5,figsize=(20,10))
# plt.axis([0, 10000, -100, 100])
# plt.ion()
# plt.show()
# fig = plt.gcf()
# fig.canvas.manager.window.raise_()

# for kk in range(len(t)):
#    tstr='time = %s years' % t[kk]
#    plt.clf()
#    plt.plot(x,b[kk,:]-binit,'b')
#    plt.axis([0, 10000, -0.5, 0.5])
#    plt.grid(True)
#    plt.title('Mass Balance')
#    plt.text(8000,0,tstr)
#    plt.draw()
#    plt.pause(0.3)        
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
# fig = plt.gcf()
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
# plt.subplot(5,1,5)
# plt.plot(t,bnet,'b')
# plt.grid(True)


plt.show()
# fig = plt.gcf()
# fig.canvas.manager.window.raise_()
######

#### Thickness deviation from initial at profiles, length, and net b_dot through time
plt.figure(8,figsize=(20,10))


plt.subplot(5,1,1)
plt.title('Profile deviation,Length,bnet')
plt.plot(t,H[:,profC]-np.mean(H[:,profC]),'g')
plt.grid(True)

plt.subplot(5,1,2)
plt.plot(t,H[:,profB]-np.mean(H[:,profB]),'b')
plt.grid(True)

plt.subplot(5,1,3)
plt.plot(t,H[:,profA]-np.mean(H[:,profA]),'r')
plt.grid(True)

plt.subplot(5,1,4)
plt.plot(t,L,'k')
plt.grid(True)

plt.subplot(5,1,5)
plt.plot(t,bnet,'k')
# plt.plot(t,Tyr-Tmean,'b')
# plt.plot(t,Pyr-Pmean[profA],'g')
plt.grid(True)
# plt.subplot(5,1,5)
# plt.plot(t,bnet,'b')
# plt.grid(True)
plt.savefig('MB_AL1917.eps', bbox_inches='tight')

# plt.show()
# fig = plt.gcf()
# fig.canvas.manager.window.raise_()
####

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
# fig = plt.gcf()
# fig.canvas.manager.window.raise_()