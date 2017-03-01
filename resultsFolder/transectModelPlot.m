% close all;clear all;
close all
set(groot,'DefaultAxesFontSize', 14)
set(groot,'DefaultTextFontSize', 14)
set(groot,'DefaultTextFontSize', 14)
% set(0,'DefaultAxesTickDir','out')
set(groot,'DefaultAxesTickDir', 'both')
set(groot,  'defaultAxesTickDirMode', 'manual');


datadir='/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/model/resultsFolder/nov29_2';
if exist('npy','var')==0
    files = dir( fullfile(datadir,'*.csv') );   %# list all *.xyz files
    files = {files.name}';                      %'# file names
    
    % data = cell(numel(files),1);                %# store file contents
    npy=struct();
    for ii=1:length(files)
        
        [pathstr,name,ext] = fileparts(files{ii});
        npy.(name)=importdata(strcat(datadir,'/',files{ii}));
        
    end
end

npy.surfinit=npy.bed+npy.Hinit;
pAind=find(npy.surfinit>=1600,1,'last');
pBind=find(npy.surfinit>=1800,1,'last');
pCind=find(npy.surfinit>=2000,1,'last');
[r,c]=size(npy.thickness);
frt=4;

pAx=[npy.xnodes(pAind) npy.xnodes(pAind)];
pBx=[npy.xnodes(pBind) npy.xnodes(pBind)];
pCx=[npy.xnodes(pCind) npy.xnodes(pCind)];
pAy=[npy.bed(pAind) npy.surfinit(pAind)+20];
pBy=[npy.bed(pBind) npy.surfinit(pBind)+20];
pCy=[npy.bed(pCind) npy.surfinit(pCind)+20];

npy.Aelev=npy.thickness(:,pAind)+npy.bed(pAind);
npy.Belev=npy.thickness(:,pBind)+npy.bed(pBind);
npy.Celev=npy.thickness(:,pCind)+npy.bed(pCind);

% figure(10);
% clf;
% hold on;
% grid on;
% box on;
% plot(npy.xnodes,npy.bed+npy.Hinit,'c','linewidth',1.5)
% plot(npy.xnodes,npy.bed,'k','linewidth',1.5)
% X = [npy.xnodes; flipud(npy.xnodes)];
% Y = [npy.bed; flipud(npy.bed+npy.Hinit)];
% fill(X,Y,'c');
% plot(npy.xnodes,npy.bed'+npy.thickness(20,:),'k','linewidth',1.5)
%





fig13=figure(13);
% fig13.PaperUnits = 'inches';
fig13.PaperPositionMode = 'auto';
set(fig13, 'Position', [30 30 1400 1000])
clf;
hold on;
grid on;
box on;

subplot(3,1,1)
hold on;
box on;
grid on;
plot(npy.time,npy.Celev,'k.','markersize',30)
plot(npy.time,npy.Celev,'k-','linewidth',1.5)
axis([0,100,1990,2010])
title('Modeled Profile Mean Elevations')
xlabel('Time (years)','fontsize',14)
% ylabel('Glacier Surface Elevation (meters)','fontsize',14)

subplot(3,1,2)
hold on;
box on;
grid on;
plot(npy.time,npy.Belev,'k.','markersize',30)
plot(npy.time,npy.Belev,'k-','linewidth',1.5)
axis([0,100,1780,1810])
xlabel('Time (years)','fontsize',14)
ylabel('Glacier Surface Elevation (meters)','fontsize',14)

subplot(3,1,3)
hold on;
box on;
grid on;
plot(npy.time,npy.Aelev,'k.','markersize',30)
plot(npy.time,npy.Aelev,'k-','linewidth',1.5)
axis([0,100,1570,1610])
xlabel('Time (years)','fontsize',14)
% ylabel('Glacier Surface Elevation (meters)','fontsize',14)

print('profile_model.eps','-depsc','-r0')


