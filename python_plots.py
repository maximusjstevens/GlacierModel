import matplotlib.pyplot as plt
## Thickness Deviation #####     
plt.figure(figsize=(20,10))
plt.axis([0, 10000, -100, 100])
plt.ion()
plt.show()
fig = plt.gcf()
fig.canvas.manager.window.raise_()

for kk in range(len(t)):
    tstr='time = %s years' % t[kk]
    plt.clf()
    plt.plot(x,Qp[kk,:]-Qp[0,:],'k')
    plt.axis([0, 10000, 0, 100])
    plt.grid(True)
    plt.title('Flux')
    plt.text(8000,20,tstr)
    plt.draw()
    plt.pause(0.7)        
######

plt.figure(figsize=(20,10))

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

