%%% Make the movie animation

close all;
%clear all;

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

[r,c]=size(npy.thickness);
npy.surfinit=npy.bed+npy.Hinit;
pAind=find(npy.surfinit>=1600,1,'last');
pBind=find(npy.surfinit>=1800,1,'last');
pCind=find(npy.surfinit>=2000,1,'last');
[r,c]=size(npy.thickness);

frt=2;

pAx=[npy.xnodes(pAind) npy.xnodes(pAind)];
pBx=[npy.xnodes(pBind) npy.xnodes(pBind)];
pCx=[npy.xnodes(pCind) npy.xnodes(pCind)];
pAy=[npy.bed(pAind) npy.surfinit(pAind)+20];
pBy=[npy.bed(pBind) npy.surfinit(pBind)+20];
pCy=[npy.bed(pCind) npy.surfinit(pCind)+20];

npy.Aelev=npy.thickness(:,pAind)+npy.bed(pAind);
npy.Belev=npy.thickness(:,pBind)+npy.bed(pBind);
npy.Celev=npy.thickness(:,pCind)+npy.bed(pCind);

set(0,'defaultAxesFontName', 'GillSans')
set(0,'defaultTextFontName', 'GillSans')

vidrecord=false;

if vidrecord
    vidObj = VideoWriter('Nisq_all_2fps.avi');
    vidObj.FrameRate = frt;
    vidObj.Quality=100;
    open(vidObj);
end

pause(0.1)

fig11=figure(11);
clf;
% box on
set(fig11,'color','white');
set(fig11, 'Units', 'normalized')
% set(fig11, 'Position', [0 0 1 1])
set(fig11, 'Position', [0 0 1 1])
bga=gca;
set(gca,'color','none','visible','off')
ha = axes('units','normalized','position',[0 0 1 1]);
set(ha,'color','none');
uistack(ha,'bottom');
I=imread('NisquallyGE.jpg');
hi = imagesc(I);
colormap gray
set(ha,'handlevisibility','off','visible','off')
set(hi,'AlphaData',0.3)
% set(gca,'nextplot','replacechildren');

% hold on;
grid on;
box on;
fs1=20;
fs2=20;
fs3=22;
co1=[0 0.5 0];
co2=[0 0 0.75];
pause(0.1)
for jj=1:r-1
    clf;
    %% main plot
    tt=(npy.bed'+npy.thickness(jj,:))';
    toeind=find(npy.thickness(jj,:)<1,1,'first');
    p1=axes('Position',[.06 .07 0.92 0.91],'Units','normalized');
    hold(p1,'on')
    plot(p1,npy.xnodes,tt,'c','linewidth',1.5)
    plot(p1,npy.xnodes,npy.bed,'k','linewidth',1.5)
    
    X = [npy.xnodes; flipud(npy.xnodes)];
    Y = [npy.bed; flipud(tt)];
    fill(X,Y,'c');
    plot(p1,npy.xnodes,npy.bed+npy.Hinit,'k--','linewidth',1.5)
    plot(pAx,pAy,'color',co1,'linewidth',4)
    plot(pBx,pBy,'color',co2,'linewidth',4)
    %     plot(pCx,pCy,co2,'linewidth',4)
    plot(npy.xnodes(toeind),npy.bed(toeind),'k.','markersize',35)
    plot(npy.xnodes(pAind),npy.Aelev(jj),'.','markersize',35,'color',co1)
    plot(npy.xnodes(pBind),npy.Belev(jj),'.','markersize',35,'color',co2)
    
    axis([5000,8000,770,2420])
    %     p1.XLim=[5000 10000];
    %     p1.YLim=[900 2200];
    %         axis equal
    set(p1,'color','none');
    set(p1,'fontsize',fs2)
    %     set(p1,'fontweight','bold')
    %     set(gca,'Units','normalized')
    %     set(gca,'Position',[.05 .05 0.93 0.93])
    xlabel('Distance (meters)','fontsize',fs1)
    ylabel('Elevation (meters)','fontsize',fs1)
    
    t = annotation('textbox');
    st={'Glacier-Model Results',sprintf('Simulation Time = %d years',jj)};
    t.String=st;
    t.FontSize = 40;
    t.FontName='GillSans';
    %     t.FontWeight='bold';
    t.LineStyle = 'none';
    t.HorizontalAlignment='center';
    t.Position=[0.25 0.8 0.1 0.1];
    
    t2 = annotation('textbox');
    t2.String='1';
    t2.FontSize = 40;
    t2.FontName='GillSans';
    t2.Color=co1;
    %     t2.FontWeight='bold';
    t2.LineStyle = 'none';
    t2.Position=[0.473 0.51 0.1 0.1];
    
    t3 = annotation('textbox');
    t3.String='2';
    t3.FontSize = 40;
    t3.FontName='GillSans';
    t3.Color=co2;
    %     t3.FontWeight='bold';
    t3.LineStyle = 'none';
    t3.Position=[0.241 0.62 0.1 0.1];
    
    t4 = annotation('textbox');
    t4.String='1';
    t4.FontSize = 40;
    t4.FontName='GillSans';
    t4.Color=co1;
    %     t4.FontWeight='bold';
    t4.LineStyle = 'none';
    t4.Position=[.581 .57 0.01 0.01];
    t4.VerticalAlignment='middle';
    
    t5 = annotation('textbox');
    t5.String='2';
    t5.FontSize = 40;
    t5.FontName='GillSans';
    t5.FontName='GillSans';
    t5.Color=co2;
    %     t5.FontWeight='bold';
    t5.LineStyle = 'none';
    t5.Position=[.581 .79 0.01 0.01];
    t5.VerticalAlignment='middle';
    
    %% text bit
    tpos=[0.1 0.2 0.4 0.2];
    fstext=40;
    if jj<49 && jj>34
        st={'When there is more melt than accumulation, the glacier retreats.'};
        tt = annotation('textbox');
        tt.String=st;
        tt.FontSize = fstext;
        tt.FontName='GillSans';
        %         tt.FontWeight='bold';
        tt.LineStyle = 'none';
        tt.Position=tpos;
        
    elseif jj > 51 && jj < 66
        st={'When there is more accumulation than melt, the glacier advances.'};
        tt = annotation('textbox');
        tt.String=st;
        tt.FontSize = fstext;
        tt.FontName='GillSans';
        %         tt.FontWeight='bold';
        tt.LineStyle = 'none';
        tt.Position=tpos;
        
    elseif jj > 68 && jj < 83
        %         st='The Nisqually Glacier has been retreating at a rate of about 20 meters per year since 1980';
        st={'Bulges of ice, called kinematic waves, can move down the glacier and cause the surface elevation to rise and fall.'};
        tt = annotation('textbox');
        tt.String=st;
        tt.FontSize = fstext;
        tt.FontName='GillSans';
        %         tt.FontWeight='bold';
        tt.LineStyle = 'none';
        tt.Position=tpos;
        
    elseif jj < 32 && jj > 17
        %         st='The Nisqually Glacier has been retreating at a rate of about 20 meters per year since 1980';
        st={'A glacier model does not reproduce or predict a glacier''s behavior exactly, but it can help scientists understand interpret the glacier data.'};
        tt = annotation('textbox');
        tt.String=st;
        tt.FontSize = fstext;
        tt.FontName='GillSans';
        %         tt.FontWeight='bold';
        tt.LineStyle = 'none';
        tt.Position=tpos;
        
    elseif jj < 15
        %         st='The Nisqually Glacier has been retreating at a rate of about 20 meters per year since 1980';
        st={'Here, a glacier model is used to simulate the Nisqually Glacier using temperature and precipitation data from Longmire.'};
        tt = annotation('textbox');
        tt.String=st;
        tt.FontSize = fstext;
        tt.FontName='GillSans';
        %         tt.FontWeight='bold';
        tt.LineStyle = 'none';
        tt.Position=tpos;
        
    elseif jj > 85
        %         st='The Nisqually Glacier has been retreating at a rate of about 20 meters per year since 1980';
        st={'A kinematic wave can be seen after time = 60 years in the simulation: the peak of the wave passes Profile 2 and then passes Profile 1 several years later'};
        tt = annotation('textbox');
        tt.String=st;
        tt.FontSize = fstext;
        tt.FontName='GillSans';
        %         tt.FontWeight='bold';
        tt.LineStyle = 'none';
        tt.Position=tpos;
    end
    %%
    
    
    
    
    
    %% length
    p2=axes('Position',[.58 .15 0.4 0.19]);
    hold(p2,'on')
    grid(p2,'on')
    box(p2,'on')
    %     plot(p2,npy.time(1:jj),npy.length(1:jj),'b.','markersize',20)
    plot(p2,npy.time(1:jj),npy.length(1:jj),'k.','markersize',30)
    plot(p2,npy.time(1:jj),npy.length(1:jj),'k-','linewidth',1.5)
    p2.XLim=[0 101];
    p2.YLim=[6500 8000];
    set(p2,'fontsize',fs2,'LineWidth',1)
    set(p2,'color','none','TickDir','both');
    xlb=xlabel('Time (years)','fontsize',fs2);
    set(xlb,'Position',get(xlb,'Position')+[0 600 0])
    ylabel('Glacier Length (meters)','fontsize',fs2)
    title('Glacier Length','fontsize',fs3)
    %     set(p2,'grid','on')
    %     p2.YLim=[2300 3300];
    %%
    
    p3a=axes('Position',[.565 .55 0.4 0.4]);
    ylabel('Surface Elevation (meters)','fontsize',fs2)
    %     set(p3a,'color','none','visible','off')
    set(p3a,'visible','off')
    %     set(p3a,'XColor','none','YTick',[])
    %     set(p3a,'YColor','none')
    %     set(p3a,'Box','off')
    set(p3a,'fontsize',fs2)
    
    set(findall(p3a, 'type', 'text'), 'visible', 'on')
    
    
    %% transect 1
    p3=axes('Position',[.58 .55 0.4 0.19]);
    hold(p3,'on')
    grid(p3,'on')
    box(p3,'on')
    %     plot(p2,npy.time(1:jj),npy.length(1:jj),'b.','markersize',20)
    plot(p3,npy.time(1:jj),npy.Aelev(1:jj),'.','markersize',30,'color',co1)
    plot(p3,npy.time(1:jj),npy.Aelev(1:jj),'-','linewidth',1.5,'color',co1)
    p3.XLim=[0 101];
    p3.YLim=[1570 1610];
    set(p3,'fontsize',fs2,'LineWidth',1)
    set(p3,'color','none','TickDir','both');
    xlabel('Time (years)','fontsize',fs2)
    if jj>85
        xs=[70 70];
        plot(xs, p3.YLim, 'k-','linewidth',3)
    end
    
    
    %     ylabel('Surface Elevation (meters)','fontsize',fs2)
    %     title('Glacier-Surface Elevation at Profile 1','fontsize',fs3)
    %     set(p2,'grid','on')
    %     p2.YLim=[2300 3300];
    %%
    
    %% transect 2
    p4=axes('Position',[.58 .77 0.4 0.19]);
    hold(p4,'on')
    grid(p4,'on')
    box(p4,'on')
    %     plot(p2,npy.time(1:jj),npy.length(1:jj),'b.','markersize',20)
    plot(p4,npy.time(1:jj),npy.Belev(1:jj),'.','markersize',30,'color',co2)
    plot(p4,npy.time(1:jj),npy.Belev(1:jj),'-','linewidth',1.5,'color',co2)
    p4.XLim=[0 101];
    p4.YLim=[1780 1810];
    if jj>85
        xs=[64 64];
        plot(xs, p4.YLim, 'k-','linewidth',3)
    end
    set(p4,'XTickLabels',[],'TickDir','both')
    set(p4,'fontsize',fs2,'LineWidth',1)
    set(p4,'color','none');
    %     xlabel('Time (years)','fontsize',fs2)
    %     ylabel('Surface Elevation (meters)','fontsize',fs2)
    title('Glacier-Surface Elevation at Profiles','fontsize',fs3)
    %     set(p2,'grid','on')
    %     p2.YLim=[2300 3300];
    %%
    pause(0.1)
    if vidrecord
        currFrame = getframe(fig11);
        writeVideo(vidObj,currFrame);
    end
    
    
    
end

if vidrecord
    close(vidObj);
end


