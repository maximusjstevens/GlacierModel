close all
set(groot,'DefaultAxesFontSize', 14)
set(groot,'DefaultTextFontSize', 14)
set(groot,'DefaultTextFontSize', 14)
% set(0,'DefaultAxesTickDir','out')
set(groot,'DefaultAxesTickDir', 'both')
set(groot,  'defaultAxesTickDirMode', 'manual');

if exist('TransA')==0
    disp('importing')
    TransA=importdata('/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/DATA/Transects/YEARLYasciisformean/TransectA_allyears_nan_SF.csv');
    TransB=importdata('/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/DATA/Transects/YEARLYasciisformean/TransectB_allyears_nan_SF.csv');
    TransC=importdata('/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/DATA/Transects/YEARLYasciisformean/TransectC_allyears_nan_SF.csv');
    
    MBin=importdata('/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/DATA/MassBalance/Al_massbalance.csv');
    MB=MBin.data;
    
    Glength_in=importdata('/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/DATA/Terminus/NisqLength151006.csv');
    Glength=Glength_in.data;
end

fig1=figure(1);
% fig1.PaperUnits = 'inches';
% fig1.PaperPosition = [0 0 12 10];
% fig1.PaperPositionMode = 'manual';
% fig1.PaperUnits = 'inches';
fig1.PaperPositionMode = 'auto';
set(fig1, 'Position', [30 30 1400 1000])
clf;
hold on;
grid on;
box on;
subplot(3,1,1)
hold on;
box on;
grid on;
plot(TransC(:,1),TransC(:,2),'b.','markersize',30)
plot(TransC(:,1),TransC(:,2),'b-','linewidth',1.5)
xlim([1930 2015])
title('Nisqually Glacier Profile Mean Elevations')
% ylabel('Elevation (m)')


subplot(3,1,2)
hold on;
grid on;
box on;
plot(TransB(:,1),TransB(:,2),'b.','markersize',30)
plot(TransB(:,1),TransB(:,2),'b-','linewidth',1.5)
xlim([1930 2015])
ylabel('Elevation (m)')

subplot(3,1,3)
hold on;
grid on;
box on;
plot(TransA(:,1),TransA(:,2),'b.','markersize',30)
plot(TransA(:,1),TransA(:,2),'b-','linewidth',1.5)
xlim([1930 2015])
% ylabel('Elevation (m)')
xlabel('Year')

print('profile_data.eps','-depsc','-r0')
