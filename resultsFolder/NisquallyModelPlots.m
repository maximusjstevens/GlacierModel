close all;clear all;

set(groot,'DefaultAxesFontName', 'Arial')
set(groot,'DefaultTextFontName', 'Arial')
set(groot,'DefaultAxesFontSize', 20)
% set(groot,'DefaultAxesLineWidth', 0.8)
set(groot,'DefaultTextFontSize', 20)
set(groot,'DefaultTextFontSize', 20)
set(groot,'DefaultLineLineWidth', 2)
% set(0,'DefaultAxesTickDir','out')
set(groot,'DefaultAxesTickDir', 'both')
set(groot,  'defaultAxesTickDirMode', 'manual');



saver1=false;
saver=false;

datadir='/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/model/resultsFolder/nov29_2';
% datadir_low='/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/model/resultsFolder/elalower';
datadir_low='/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/model/resultsFolder/elalow11';
datadir_high='/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/model/resultsFolder/elahigher';
datadir_saw='/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/model/resultsFolder/elasaw';
datadir_noise='/Users/maxstev/Documents/Grad_School/Nisqually/Nisqually/model/resultsFolder/elarandom';

if exist('npy','var')==0
    files = dir( fullfile(datadir,'*.csv') );   %# list all *.xyz files
    files = {files.name}';                      %'# file names
    
    % data = cell(numel(files),1);                %# store file contents
    npy=struct();
    for ii=1:length(files)
        
        [pathstr,name,ext] = fileparts(files{ii});
        npy.(name)=importdata(strcat(datadir,'/',files{ii}));
        
    end
    npy.surfinit=npy.bed+npy.Hinit;
    pAind=find(npy.surfinit>=1600,1,'last');
    pBind=find(npy.surfinit>=1800,1,'last');
    pCind=find(npy.surfinit>=2000,1,'last');
    npy.Aelev=npy.thickness(:,pAind)+npy.bed(pAind);
    npy.Belev=npy.thickness(:,pBind)+npy.bed(pBind);
    npy.Celev=npy.thickness(:,pCind)+npy.bed(pCind);
    npy.year=(1915:2015)';
end

if exist('npylow','var')==0
    files = dir( fullfile(datadir_low,'*.csv') );   %# list all *.xyz files
    files = {files.name}';                      %'# file names
    
    % data = cell(numel(files),1);                %# store file contents
    npylow=struct();
    for ii=1:length(files)
        
        [pathstr,name,ext] = fileparts(files{ii});
        npylow.(name)=importdata(strcat(datadir_low,'/',files{ii}));
        
    end
    npylow.surfinit=npylow.bed+npylow.Hinit;
    pAindela=find(npylow.surfinit>=1600,1,'last');
    pBindela=find(npylow.surfinit>=1800,1,'last');
    pCindela=find(npylow.surfinit>=2000,1,'last');
    npylow.Aelev=npylow.thickness(:,pAindela)+npylow.bed(pAindela);
    npylow.Belev=npylow.thickness(:,pBindela)+npylow.bed(pBindela);
    npylow.Celev=npylow.thickness(:,pCindela)+npylow.bed(pCindela);
    
    
end

if exist('npyhigh','var')==0
    files = dir( fullfile(datadir_high,'*.csv') );   %# list all *.xyz files
    files = {files.name}';                      %'# file names
    
    % data = cell(numel(files),1);                %# store file contents
    npyhigh=struct();
    for ii=1:length(files)
        
        [pathstr,name,ext] = fileparts(files{ii});
        npyhigh.(name)=importdata(strcat(datadir_high,'/',files{ii}));
        
    end
    npyhigh.surfinit=npyhigh.bed+npyhigh.Hinit;
    npyhigh.Aelev=npyhigh.thickness(:,pAindela)+npyhigh.bed(pAindela);
    npyhigh.Belev=npyhigh.thickness(:,pBindela)+npyhigh.bed(pBindela);
    npyhigh.Celev=npyhigh.thickness(:,pCindela)+npyhigh.bed(pCindela);
end

if exist('npysaw','var')==0
    files = dir( fullfile(datadir_saw,'*.csv') );   %# list all *.xyz files
    files = {files.name}';                      %'# file names
    
    % data = cell(numel(files),1);                %# store file contents
    npysaw=struct();
    for ii=1:length(files)
        
        [pathstr,name,ext] = fileparts(files{ii});
        npysaw.(name)=importdata(strcat(datadir_saw,'/',files{ii}));
        
    end
    npysaw.surfinit=npysaw.bed+npysaw.Hinit;
    npysaw.Aelev=npysaw.thickness(:,pAindela)+npysaw.bed(pAindela);
    npysaw.Belev=npysaw.thickness(:,pBindela)+npysaw.bed(pBindela);
    npysaw.Celev=npysaw.thickness(:,pCindela)+npysaw.bed(pCindela);
end


if exist('npynoise','var')==0
    files = dir( fullfile(datadir_noise,'*.csv') );   %# list all *.xyz files
    files = {files.name}';                      %'# file names
    
    % data = cell(numel(files),1);                %# store file contents
    npynoise=struct();
    for ii=1:length(files)
        
        [pathstr,name,ext] = fileparts(files{ii});
        npynoise.(name)=importdata(strcat(datadir_noise,'/',files{ii}));
        
    end
    npynoise.surfinit=npynoise.bed+npynoise.Hinit;
    npynoise.Aelev=npynoise.thickness(:,pAindela)+npynoise.bed(pAindela);
    npynoise.Belev=npynoise.thickness(:,pBindela)+npynoise.bed(pBindela);
    npynoise.Celev=npynoise.thickness(:,pCindela)+npynoise.bed(pCindela);
end


% [r,c]=size(npy.thickness);


xl=[1930 2015];
yrvc=1930:5:2015;
% xtks=datenum(yrvc,ones(size(yrvc)),ones(size(yrvc)));
xtks=yrvc;
ppos=[0.1 0.1 0.9 0.7];

col1=[12 195 82] ./ 255;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig6=figure(6);
clf;
set(fig6, 'Units', 'normalized')
set(fig6, 'Position', [0.1 0.0 0.9 0.9])
fig6.PaperPositionMode = 'auto';
ax0=gca;
set(ax0,'Position',[.13 .11 0.775  0.79],'Units','normalized');
set(ax0,'color','w')
set(ax0,'box','on','LineWidth',0.8,'XColor','k','YColor','k')
set(ax0,'YGrid','off','XGrid','off','XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[])

ax1=axes;
set(ax1,'color','none')
set(ax1,'Box','off')
set(ax1,'YGrid','off','XGrid','on')
set(ax1,'Position',[.13 .3 0.775  0.6],'Units','normalized');
set(ax1,'YAxisLocation','right')
set(ax1,'XAxisLocation','top')
% set(ax1,'YLim',yls)
hold(ax1,'on')

a1=plot(npy.year,npy.Aelev-mean(npy.Aelev),'--','Color',col1,'linewidth',2);
a3=plot(npy.year,npy.Aelev-mean(npy.Aelev),'*','Color',col1,'markersize',5);

b1=plot(npy.year,npy.Belev-mean(npy.Belev),'--','Color','b','linewidth',2);
b3=plot(npy.year,npy.Belev-mean(npy.Belev),'*','Color','b','markersize',5);

c1=plot(npy.year,npy.Celev-mean(npy.Celev),'--','Color','r','linewidth',2);
c3=plot(npy.year,npy.Celev-mean(npy.Celev),'*','Color','r','markersize',5);

xlim(xl)
set(ax1,'XTick',xtks)
% datetick('x','keeplimits','keepticks')
set(ax1,'XTickLabel',[])
% set(ax1,'YTickLabel',{'-60','-40','-20','0','20','40'})
ylabel('Profile elevation deviation from mean (m)')


co3=[0.4 0.4 0.4];
ax3=axes('Position',[0.13 0.11 0.775 0.15],'Units','normalized');
set(ax3,'Box','off','Color','none')
set(ax3,'YGrid','off','XGrid','on')
hold(ax3,'on')
Tp3=bar(npy.year,npy.bnet,'EdgeColor',co3,'FaceColor',co3);
xlim(xl)
ylim([-25 25])
set(ax3,'XTick',xtks,'YColor',co3)
% set(ax3,'YTick',[],'YTickLabel',[])
set(ax3,'YAxisLocation','right')
ylabel('Mass Balance (m)','fontsize',16,'Color',co3)
xlabel('Year')
% ax3=ax0;
% set(ax3,'Box','on','linewidth',0.8)
% set(ax3,'Color','none')

% ax2=axes('Position',[0.13 0.11 0.775 0.29],'Units','normalized');
% ax2=axes('Position',[0.13 0.11 0.775 0.45],'Units','normalized');
ax2=axes('Position',[0.13 0.11 0.775 0.38],'Units','normalized');
set(ax2,'Box','off','Color','none')
set(ax2,'YGrid','off','XGrid','on')
hold(ax2,'on')
Tp1=plot(npy.year,npy.length,'k--','linewidth',2);
Tp2=plot(npy.year,npy.length,'k*','markersize',5);
xlim(xl)
ylim([6500 7600])
set(ax2,'XTick',xtks)
set(ax2,'YTick',[6800 7000 7200 7400 7600])
set(ax2,'YTickLabel',[6800 7000 7200 7400 7600])
% datetick('x','keeplimits','keepticks')
% xlabel('Year')
ylabel('             Glacier length (m)')
% leg1=legend([a1 b1 c1 Tp1],{'Profile 1/A','Profile 2/B','Profile 3/C','Terminus'});
% set(leg1,'Color','w')

% tb1=annotation('textbox',[0.18 0.16 0.3 0.1],'String','Relative mass balance' ,'FitBoxToText','on');
% tb1.EdgeColor='none';
print('ModelProfiles_ClimateForcing_MB.eps','-depsc','-r0')

if saver1
    print('ModelProfiles_ClimateForcing.eps','-depsc','-r0')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xli=[8 24];
fig1=figure(1);
clf;
set(fig1, 'Units', 'normalized')
set(fig1, 'Position', [0.0 0.0 1 1])
fig1.PaperPositionMode = 'auto';
ax0=gca;
set(ax0,'Position',[.13 .11 0.775  0.79],'Units','normalized');
set(ax0,'color','w')
set(ax0,'box','on','LineWidth',0.8)
set(ax0,'YGrid','off','XGrid','off','XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[])

ax1=axes;
set(ax1,'color','w')
set(ax1,'Box','on','LineWidth',1.5)
set(ax1,'YGrid','on','XGrid','on')
set(ax1,'Position',[.13 .74 0.775  0.16],'Units','normalized');
set(ax1,'YAxisLocation','left')
set(ax1,'XAxisLocation','bottom')
set(ax1,'XTickLabels',[])
% set(ax1,'YLim',yls)
hold(ax1,'on')
plot(npylow.time,npylow.Celev);
plot(npyhigh.time,npyhigh.Celev);
plot(npysaw.time,npysaw.Celev);
xlim(xli)
ylim([1997 2005])
leg1=legend('Experiment 1','Experiment 2','Experiment 3');
ylabel('C Elevation (m)')

ax2=axes;
set(ax2,'color','w')
set(ax2,'Box','on','LineWidth',1.5)
set(ax2,'YGrid','on','XGrid','on')
set(ax2,'Position',[.13 .58 0.775  0.16],'Units','normalized');
set(ax2,'YAxisLocation','right')
set(ax2,'XAxisLocation','bottom')
set(ax2,'XTickLabels',[])
hold(ax2,'on')
plot(npylow.time,npylow.Belev);
plot(npyhigh.time,npyhigh.Belev);
plot(npysaw.time,npysaw.Belev);
xlim(xli)
ylim([1796 1804])
ylabel('B Elevation (m)')

ax3=axes;
set(ax3,'color','w')
set(ax3,'Box','on','LineWidth',1.5)
set(ax3,'YGrid','on','XGrid','on')
set(ax3,'Position',[.13 .42 0.775  0.16],'Units','normalized');
set(ax3,'YAxisLocation','left')
set(ax3,'XAxisLocation','bottom')
set(ax3,'XTickLabels',[])
hold(ax3,'on')
plot(npylow.time,npylow.Aelev);
plot(npyhigh.time,npyhigh.Aelev);
plot(npysaw.time,npysaw.Aelev);
plot(npylow.time,npylow.Aelev+npyhigh.Aelev-npyhigh.Aelev(1),'k')
xlim(xli)
ylim([1596 1604])
ylabel('A Elevation (m)')

ax4=axes;
set(ax4,'color','w')
set(ax4,'Box','on','LineWidth',1.5)
set(ax4,'YGrid','on','XGrid','on')
set(ax4,'Position',[.13 .26 0.775  0.16],'Units','normalized');
set(ax4,'YAxisLocation','right')
set(ax4,'XAxisLocation','bottom')
set(ax4,'XTickLabels',[])
hold(ax4,'on')
plot(npylow.time,npylow.length);
plot(npyhigh.time,npyhigh.length);
plot(npysaw.time,npysaw.length);
xlim(xli)
ylabel('Length')

ax5=axes;
set(ax5,'color','w')
set(ax5,'Box','on','LineWidth',1.5)
set(ax5,'YGrid','on','XGrid','on')
set(ax5,'Position',[.13 .1 0.775  0.16],'Units','normalized');
set(ax5,'YAxisLocation','left')
set(ax5,'XAxisLocation','bottom')
hold(ax5,'on')
plot(npylow.time,npylow.bnet);
plot(npyhigh.time,npyhigh.bnet);
plot(npysaw.time,npysaw.bnet);
xlim(xli)
ylabel('Mass Balance (m)')
xlabel('Time (years)')
if saver
    print('E123_prof_2.eps','-depsc','-r0')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xli=[0 25];
colr=[0.494 0.184 0.556];
fig2=figure(2);
clf;
set(fig2, 'Units', 'normalized')
set(fig2, 'Position', [0.0 0.0 1 1])
fig2.PaperPositionMode = 'auto';
ax0=gca;
set(ax0,'Position',[.13 .11 0.775  0.79],'Units','normalized');
set(ax0,'color','w')
set(ax0,'box','on','LineWidth',0.8)
set(ax0,'YGrid','off','XGrid','off','XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[])

ax1=axes;
set(ax1,'color','w')
set(ax1,'Box','on','LineWidth',1.5)
set(ax1,'YGrid','on','XGrid','on')
set(ax1,'Position',[.13 .74 0.775  0.16],'Units','normalized');
set(ax1,'YAxisLocation','left')
set(ax1,'XAxisLocation','bottom')
set(ax1,'XTickLabels',[])
% set(ax1,'YLim',yls)
hold(ax1,'on')
plot(npynoise.time,npynoise.Celev,'color',colr);
plot(npyhigh.time,npyhigh.Celev);
plot(npysaw.time,npysaw.Celev);
xlim(xli)
ylim([1985 2015])
leg1=legend('Experiment 1','Experiment 2','Experiment 3');
ylabel('C Elevation (m)')


ax2=axes;
set(ax2,'color','w')
set(ax2,'Box','on','LineWidth',1.5)
set(ax2,'YGrid','on','XGrid','on')
set(ax2,'Position',[.13 .58 0.775  0.16],'Units','normalized');
set(ax2,'YAxisLocation','right')
set(ax2,'XAxisLocation','bottom')
set(ax2,'XTickLabels',[])
hold(ax2,'on')
plot(npynoise.time,npynoise.Belev,'color',colr);
plot(npyhigh.time,npyhigh.Belev);
plot(npysaw.time,npysaw.Belev);
xlim(xli)
ylim([1785 1815])
ylabel('B Elevation (m)')

ax3=axes;
set(ax3,'color','w')
set(ax3,'Box','on','LineWidth',1.5)
set(ax3,'YGrid','on','XGrid','on')
set(ax3,'Position',[.13 .42 0.775  0.16],'Units','normalized');
set(ax3,'YAxisLocation','left')
set(ax3,'XAxisLocation','bottom')
set(ax3,'XTickLabels',[])
hold(ax3,'on')
plot(npynoise.time,npynoise.Aelev,'color',colr);
plot(npyhigh.time,npyhigh.Aelev);
plot(npysaw.time,npysaw.Aelev);
xlim(xli)
ylim([1585 1615])
ylabel('A Elevation (m)')

ax4=axes;
set(ax4,'color','w')
set(ax4,'Box','on','LineWidth',1.5)
set(ax4,'YGrid','on','XGrid','on')
set(ax4,'Position',[.13 .26 0.775  0.16],'Units','normalized');
set(ax4,'YAxisLocation','right')
set(ax4,'XAxisLocation','bottom')
set(ax4,'XTickLabels',[])
hold(ax4,'on')
plot(npynoise.time,npynoise.length,'color',colr);
plot(npyhigh.time,npyhigh.length);
plot(npysaw.time,npysaw.length);
xlim(xli)
ylabel('Length')

ax5=axes;
set(ax5,'color','w')
set(ax5,'Box','on','LineWidth',1.5)
set(ax5,'YGrid','on','XGrid','on')
set(ax5,'Position',[.13 .1 0.775  0.16],'Units','normalized');
set(ax5,'YAxisLocation','left')
set(ax5,'XAxisLocation','bottom')
hold(ax5,'on')
plot(npynoise.time,npynoise.bnet,'color',colr);
plot(npyhigh.time,npyhigh.bnet);
plot(npysaw.time,npysaw.bnet);
xlim(xli)
ylabel('Mass Balance (m)')
xlabel('Time (years)')
if saver
    print('E4_prof_2.eps','-depsc','-r0')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xli=[8 24];
fig3=figure(3);
clf;
set(fig3, 'Units', 'normalized')
set(fig3, 'Position', [0.0 0.0 1 1])
fig3.PaperPositionMode = 'auto';
% ax0=gca;
% set(ax0,'Position',[.13 .11 0.775  0.79],'Units','normalized');
% set(ax0,'color','w')
% set(ax0,'box','on','LineWidth',0.8)
% set(ax0,'YGrid','off','XGrid','off','XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[])

ax1=axes;
set(ax1,'color','w')
set(ax1,'Box','on','LineWidth',1.5)
set(ax1,'YGrid','on','XGrid','on')
set(ax1,'Position',[.13 .74 0.775  0.16],'Units','normalized');
set(ax1,'YAxisLocation','left')
set(ax1,'XAxisLocation','bottom')
set(ax1,'XTickLabels',[])
% set(ax1,'YLim',yls)
hold(ax1,'on')
plot(npylow.time,npylow.Celev);
plot(npyhigh.time,npyhigh.Celev);
plot(npysaw.time,npysaw.Celev);
plot(npynoise.time,npynoise.Celev);
xlim(xli)
leg1=legend('Experiment 1','Experiment 2','Experiment 3','Experiment 4');
ylabel('C Elevation (m)')


ax2=axes;
set(ax2,'color','w')
set(ax2,'Box','on','LineWidth',1.5)
set(ax2,'YGrid','on','XGrid','on')
set(ax2,'Position',[.13 .58 0.775  0.16],'Units','normalized');
set(ax2,'YAxisLocation','right')
set(ax2,'XAxisLocation','bottom')
set(ax2,'XTickLabels',[])
hold(ax2,'on')
plot(npylow.time,npylow.Belev);
plot(npyhigh.time,npyhigh.Belev);
plot(npysaw.time,npysaw.Belev);
plot(npynoise.time,npynoise.Belev);
xlim(xli)
ylabel('B Elevation (m)')

ax3=axes;
set(ax3,'color','w')
set(ax3,'Box','on','LineWidth',1.5)
set(ax3,'YGrid','on','XGrid','on')
set(ax3,'Position',[.13 .42 0.775  0.16],'Units','normalized');
set(ax3,'YAxisLocation','left')
set(ax3,'XAxisLocation','bottom')
set(ax3,'XTickLabels',[])
hold(ax3,'on')
plot(npylow.time,npylow.Aelev);
plot(npyhigh.time,npyhigh.Aelev);
plot(npysaw.time,npysaw.Aelev);
plot(npynoise.time,npynoise.Aelev);
xlim(xli)
ylabel('A Elevation (m)')

ax4=axes;
set(ax4,'color','w')
set(ax4,'Box','on','LineWidth',1.5)
set(ax4,'YGrid','on','XGrid','on')
set(ax4,'Position',[.13 .26 0.775  0.16],'Units','normalized');
set(ax4,'YAxisLocation','right')
set(ax4,'XAxisLocation','bottom')
set(ax4,'XTickLabels',[])
hold(ax4,'on')
plot(npylow.time,npylow.length);
plot(npyhigh.time,npyhigh.length);
plot(npysaw.time,npysaw.length);
plot(npynoise.time,npynoise.length);
xlim(xli)
ylabel('Length')

ax5=axes;
set(ax5,'color','w')
set(ax5,'Box','on','LineWidth',1.5)
set(ax5,'YGrid','on','XGrid','on')
set(ax5,'Position',[.13 .1 0.775  0.16],'Units','normalized');
set(ax5,'YAxisLocation','left')
set(ax5,'XAxisLocation','bottom')
hold(ax5,'on')
plot(npylow.time,npylow.bnet);
plot(npyhigh.time,npyhigh.bnet);
plot(npysaw.time,npysaw.bnet);
pX=plot(npynoise.time,npynoise.bnet);
xlim(xli)
ylabel('Mass Balance (m)')
xlabel('Time (years)')
if saver
    print('E1234_prof.eps','-depsc','-r0')
end

% %%%%%%
% figure(4)
% clf
% hold on;
% plot(npylow.xnodes(1:end),npylow.thickness(11,1:end)-(npylow.Hinit(1:end))')
% plot(npylow.xnodes(1:end),npylow.thickness(12,1:end)-(npylow.Hinit(1:end))')
% plot(npylow.xnodes(1:end),npylow.thickness(13,1:end)-(npylow.Hinit(1:end))')
% plot(npylow.xnodes(1:end),npylow.thickness(14,1:end)-(npylow.Hinit(1:end))')
% plot(npylow.xnodes(1:end),npylow.thickness(15,1:end)-(npylow.Hinit(1:end))')
% plot(npylow.xnodes(1:end),npylow.thickness(16,1:end)-(npylow.Hinit(1:end))')
% plot(npylow.xnodes(1:end),npylow.thickness(17,1:end)-(npylow.Hinit(1:end))')
% legend('11','12','13','14','15','16','17')
% ylim([0 4])
% % 
% % for ii=5:20
% %     disp(ii)
% %     plot(npylow.xnodes(600:630),npylow.thickness(ii,600:630)-(npylow.Hinit(600:630))')
% %     ylim([-1 15])
% %     pause
% % end
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     