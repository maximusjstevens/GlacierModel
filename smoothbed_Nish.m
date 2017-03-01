bedprof=importdata('/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/model/NisqBed.csv');
xx=bedprof(:,1);
yy=bedprof(:,2);

delb = 15.e3/log(4.);
zb =  4300.*exp(-xx.^1.05/delb);

zz=3.26188711e-05*xx.^2 +  -6.48240291e-01*xx +  4.33066591e+03;


figure(1);
clf;
hold on;
plot(xx,yy,'b')
plot(xx,zb,'r')
plot(xx,zz,'g')