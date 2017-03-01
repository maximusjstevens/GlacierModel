close all;clear all;

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
vidrecord=false;

if vidrecord
    vidObj = VideoWriter('Nisq_term.avi');
    vidObj.FrameRate = frt;
    vidObj.Quality=100;
    open(vidObj);
end

fig11=figure(11);
clf;
box on
set(fig11,'color','white');
set(fig11, 'Units', 'normalized')
set(fig11, 'Position', [0 0 0.90 0.80])
set(gca,'nextplot','replacechildren');
grid on;
for jj=1:r
    clf;
    hold on
    box on
    tt=(npy.bed'+npy.thickness(jj,:))';
    plot(npy.xnodes,tt,'c','linewidth',1.5)
    plot(npy.xnodes,npy.bed,'k','linewidth',1.5)

    X = [npy.xnodes; flipud(npy.xnodes)];
    Y = [npy.bed; flipud(tt)];
    fill(X,Y,'c');
    plot(npy.xnodes,npy.bed+npy.Hinit,'k--','linewidth',1.5)
    plot(pAx,pAy,'k','linewidth',4)
    plot(pBx,pBy,'k','linewidth',4)
    plot(pCx,pCy,'k','linewidth',4)
    %     axis([5500,7000,1500,2000])
    xlim([5700,7700])
    %     axis([4500,5000,1800,2100])
    axis equal
    set(gca,'color','white');
    set(gca,'fontsize',14)
    set(gca,'Units','normalized')
    set(gca,'Position',[.05 .05 0.93 0.93])
    xlabel('Distance (meters)','fontsize',14)
    ylabel('Elevation (meters)','fontsize',14)
    
    st=sprintf('Time = %d years',jj);
    t = annotation('textbox');
    t.String=st;
    t.FontSize = 20;
    t.FontWeight='bold';
    t.LineStyle = 'none';
    
    if vidrecord
        currFrame = getframe(fig11);
        writeVideo(vidObj,currFrame);
    end
        %pause(0.1)
    
    %     filename = sprintf('images/image%04d.jpg',jj);
    %     export_fig(gcf,filename,'-r72', '-q95', '-transparent', '-nocrop');
    
    
end

if vidrecord
    close(vidObj);
end

if vidrecord
    vidObj = VideoWriter('Nisq_len.avi');
    vidObj.FrameRate = 4;
    vidObj.Quality=100;
    open(vidObj);
end

fig12=figure(12);
clf;
box on
set(fig12,'color','white');
set(fig12, 'Units', 'normalized')
set(fig12, 'Position', [0 0 0.90 0.50])
set(gca,'nextplot','replacechildren');
grid on;
for jj=1:r
    clf;
    hold on
    box on
    grid on
    plot(npy.time(1:jj),npy.length(1:jj),'k.','markersize',30)
    plot(npy.time(1:jj),npy.length(1:jj),'k-','linewidth',1.5)
    axis([0,101,6500,8000])
    set(gca,'color','white');
    set(gca,'fontsize',14)
    set(gca,'Units','normalized')
    set(gca,'Position',[.1 .1 0.88 0.88])
    xlabel('Time (years)','fontsize',14)
    ylabel('Glacier Length (meters)','fontsize',14)
    
    
    %     st=sprintf('Time = %d years',jj);
    %     t = annotation('textbox');
    %     t.String=st;
    %     t.FontSize = 20;
    %     t.FontWeight='bold';
    %     t.LineStyle = 'none';
    if vidrecord
        currFrame = getframe(fig12);
        writeVideo(vidObj,currFrame);
    end
       %pause(0.1)
    
    %     filename = sprintf('images/image%04d.jpg',jj);
    %     export_fig(gcf,filename,'-r72', '-q95', '-transparent', '-nocrop');
    
    
end

if vidrecord
    close(vidObj);
end

if vidrecord
    vidObj = VideoWriter('Nisq_ProfA.avi');
    vidObj.FrameRate = 4;
    vidObj.Quality=100;
    open(vidObj);
end

fig13=figure(13);
clf;
box on
set(fig13,'color','white');
set(fig13, 'Units', 'normalized')
set(fig13, 'Position', [0 0 0.90 0.50])
set(gca,'nextplot','replacechildren');
grid on;
for jj=1:r
    clf;
    hold on
    box on
    grid on
    tt=(npy.bed(pAind)+npy.thickness(jj,pAind))';
    plot(npy.time(1:jj),npy.Aelev(1:jj),'k.','markersize',30)
    plot(npy.time(1:jj),npy.Aelev(1:jj),'k-','linewidth',1.5)
    axis([0,101,1570,1610])
    set(gca,'color','white');
    set(gca,'fontsize',14)
    set(gca,'Units','normalized')
    set(gca,'Position',[.1 .1 0.88 0.88])
    xlabel('Time (years)','fontsize',14)
    ylabel('Glacier Surface Elevation (meters)','fontsize',14)
    
    
    %     st=sprintf('Time = %d years',jj);
    %     t = annotation('textbox');
    %     t.String=st;
    %     t.FontSize = 20;
    %     t.FontWeight='bold';
    %     t.LineStyle = 'none';
    if vidrecord
        currFrame = getframe(fig13);
        writeVideo(vidObj,currFrame);
    end
        %pause(0.1)

    
end

if vidrecord
    close(vidObj);
end

if vidrecord
    vidObj = VideoWriter('Nisq_profB.avi');
    vidObj.FrameRate = 4;
    vidObj.Quality=100;
    open(vidObj);
end


fig14=figure(14);
clf;
box on
set(fig14,'color','white');
set(fig14, 'Units', 'normalized')
set(fig14, 'Position', [0 0 0.90 0.50])
set(gca,'nextplot','replacechildren');
grid on;
for jj=1:r
    clf;
    hold on
    box on
    grid on
    tt=(npy.bed(pAind)+npy.thickness(jj,pAind))';
    plot(npy.time(1:jj),npy.Belev(1:jj),'k.','markersize',30)
    plot(npy.time(1:jj),npy.Belev(1:jj),'k-','linewidth',1.5)
    axis([0,101,1780,1810])
    set(gca,'color','white');
    set(gca,'fontsize',14)
    set(gca,'Units','normalized')
    set(gca,'Position',[.1 .1 0.88 0.88])
    xlabel('Time (years)','fontsize',14)
    ylabel('Glacier Surface Elevation (meters)','fontsize',14)
    
    
    %     st=sprintf('Time = %d years',jj);
    %     t = annotation('textbox');
    %     t.String=st;
    %     t.FontSize = 20;
    %     t.FontWeight='bold';
    %     t.LineStyle = 'none';
    if vidrecord
        currFrame = getframe(fig14);
        writeVideo(vidObj,currFrame);
    end
        %pause(0.1)

    
end

if vidrecord
    close(vidObj);
end

