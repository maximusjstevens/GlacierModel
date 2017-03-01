% bed slopes: create the bed profile (obsolete 5/8/13)

clear all; close all;

xmx = 20000;          % max of the domain size in m
delx = 100;              % grid spacing in m
nxs = round(xmx/delx) + 1;  % number of grid points
x = 0:delx:xmx;  

% del = 15e3/log(3);
% zb =  3000.*exp(-x/del);     % bed profile in m  
%zb = (-0.1*x+1000);  

blah=erf(x-10000);

figure(1);clf;

plot(x,blah);