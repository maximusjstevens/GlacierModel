clear all; 
close all;
tic
%% --------------------------------
%define parameters and constants
%----------------------------------
rho=917;              % density of ice in kg/m^3
g=9.8;                  % gravity in m/s^2
n=3;                    % glen's flow law constant

s_per_year = 365.35*24*3600;        % number of seconds in a year (s yr^-1)
A_T = 6.8e-24*s_per_year;       % softness parameter for -5 degree ice (yr^-1 Pa^-3)
fd = 2*A_T/(n+2);                           % constants lumped together, and put into standard units.  (Pa^-3 yr^-1)
fs = 5.7e-20*s_per_year;                % sliding parameter fromOerlemans; (Pa^-3 m^2 yr-1)

%% --------------------------------
%define grid spacing and time
%----------------------------------

xmx = 20000;          % max of the domain size in m
delx = 100;              % grid spacing in m
nxs = round(xmx/delx) + 1;  % number of grid points
x = 0:delx:xmx;                   % x array (each point)

delt = 0.0125; % time step in yrs 
ts = 0;             % starting time in yrs
tf = 500;      % final time in yrs
nts=floor((tf-ts)/delt) +1; % number of time steps ('round' used because need nts as integer)
nyrs = tf-ts;       % just final time - start time

%% -----------------
% define accumulation/ablation
%-----------------

% climate foricing for steady state
b = 3-(4/5000)*x;      % mass balance in /yr

%% ---------------------------------
% glacier bed geometries
%---------------------------------
del = 15e3/log(3);
zb =  3000.*exp(-x/del);     % bed profile in m  
%zb = (-0.1*x+1000);  

%% -----------------------
% initialize arrays
%------------------------
% load glac_init_class.mat
% H0  = interp1(x_init,h_init,x,'linear');

H0 = ones(size(x));     % initial thickness of ice (in m)--do not confuse with h, ice elevation.  h = H+zb
Qp=zeros(size(x));      % Qp equals j+1/2 flux
Qm=zeros(size(x));      % Qm equals j-1/2 flux
up = zeros(size(x));    % velocity at j +1/2 grid 
um = zeros(size(x));    % velcotiy at j-1/2 grid
dHdt=zeros(size(x));        %  slope of ice thickness -- do not confuse with dhdx, surface slope.  dhdx = dzdx+dHdx
H = H0;                    % initialize height array

yr = 0;         % for counting years
idx_out = 0;    % 
deltout = 1;    % for counting years
nouts = round(nts/1);   % for counting years
edge_out = zeros(length(nouts));    % for counting
t_out = zeros(length(nouts));           % for counting

%% -----------------------------------------
% begin loop over time
%-----------------------------------------

 for i=1:nts

     
    t = delt*(i-1); % time in years 
    if i==nts
        disp(t)
    end
    %-----------------------------------------
% begin loop over space
%-----------------------------------------
%    for j=1:nxs-1  % 
       
%        if j==1      % at the glacier divide

                H_ave_d =(H(1) + H(2))/2;          % average glacier thickness between divide and second point
                dHdx_d = (H(2) - H(1))/delx;      % linear slope of glacier thickness between divide and second point
                dzbdx_d = (zb(2) - zb(1))/delx;     % linear slope of bed between divide and second point
                dhdx_d = dHdx_d+dzbdx_d;             % linear surface slope of ice between divide and second point
            
                Qp(1) = -(dhdx_d).^3 * H_ave_d.^4*(rho*g)^3*(fd*H_ave_d+fs/H_ave_d);       % flux at plus half grid point
                up(1) = Qp(1)/H_ave_d;    % velocity at half grid point
            
                Qm(1) = 0;                      % flux at minus half grid point = 0 because we start at the divide--nothing flowing in
                um(1) = Qm(1)/H(1);       % velocity is also zero coming in
            
                dHdt(1) = b(1) - Qp(1)/(delx/2);    % change in thickness at the divide for this timestep is based on snowfall in/out, flux out
           
%        elseif H(j)==0 && H(j-1)>0 % at the glacier toe
         ind1=find(H<=0, 1 );
         if isempty(ind1)==1
             ind1=nxs;
         end
        
                Qp(ind1) = 0;   % nothing to flux out of the end
                up(ind1) = 0;   % no velocity b/c no ice
           
                H_ave_toe = (0+H(ind1-1))/2;              % thickness of glacier is average of  last point with thickness and a thickness of 0     
                dHdx_toe = 0-H(ind1-1)/delx;              % slope is based on last point with thickness and a thickness of 0
                dzbdx_toe =(zb(ind1)-zb(ind1-1))/delx;       % bed slope between toe and last point with thickness
                dhdx_toe = dHdx_toe+dzbdx_toe;            % surface slope between toe and last point with thickness
           
                Qm(ind1) = -(rho*g)^3*H_ave_toe.^4.*(dhdx_toe).^3.*(fd.*H_ave_toe+fs./H_ave_toe);    % flux coming from the last point with thickness and melting out before it reaches a new grid point
                um(ind1) = Qm(ind1)/H_ave_toe;         % velocity at the toe
           
                dHdt(ind1) = b(ind1) -(Qp(ind1)- Qm(ind1))/delx;     % change in height of the glacier at the toe based on the snowfall in/melt out and the flux in/out
                edge = ind1; 	%index of glacier toe - used for fancy plotting
         
           
%        elseif H(j)<=0 && H(j-1)<=0 % beyond glacier toe - no glacier flux
            ind2=find(H<=0);
            ind2=ind2(2:end);
                dHdt(ind2) = b(ind2);  % thickness change is only from any change in accumulaiton in/melt out
           
                % no flux in, no flux out
                Qp(ind2) = 0; 
                Qm(ind2) = 0;
           
                um(ind2) =0; 
                up(ind2) = 0;
            
%        else  % within the glacier
            ind3=2:ind1-1;
           % first for the flux going out
                H_ave_o = (H(ind3+1) + H(ind3))/2;              % the average thickness between j and j+1
                dHdx_o = (H(ind3+1) - H(ind3))/delx;           % the linear ice thickness slope between j and j+1
                dzbdx_o = (zb(ind3+1) - zb(ind3))/delx;       % the linear bed slope between j and j+1
                dhdx_o = dHdx_o + dzbdx_o;                      % the surface slope between j and j+1
           
                Qp(ind3) = -(rho*g)^3*H_ave_o.^4.*(dhdx_o).^3.*(fd.*H_ave_o+fs./H_ave_o);  % flux coming out
                up(ind3) = Qp(ind3)./H_ave_o;                           % velocity coming out
           
           % now for the flux coming in
                H_ave_i = (H(ind3-1) + H(ind3))/2;               % the average thickness between j and j-1
                dHdx_i = (H(ind3) - H(ind3-1))/delx;            % the linear ice thickness slope between j and j-1
                dzbdx_i = (zb(ind3) - zb(ind3-1))/delx;       % the linear bed slope between j and j-1
                dhdx_i = dHdx_i + dzbdx_i;                      % the surface slope between j and j-1

                Qm(ind3) = -(rho*g)^3.*H_ave_i.^4.*dhdx_i.^3.*(fd.*H_ave_i+fs./H_ave_i);        % flux coming in
                um(ind3) = Qm(ind3)/H_ave_i;            % velocity coming in
           
           % change in height at point j for this timestep
           dHdt(ind3) = b(ind3) - (Qp(ind3) - Qm(ind3))/delx;       % based on accumulation in/melt out and flux in/out
           
       
%        end  
       
       dHdt(nxs) = 0; % enforce no change at edge of domian


%    end  % done with each point in space
    
   % For the next timestep, the thickness is equal to the thickness profile from the previous timstep, plus the amount that the thickness has changed at each point
   % over this time period, (dHdt*delt).  The ice thickness at each point
   % needs to be positive, so we choose the max of 0 and the summation.
   % since dHdt can be negative, we want to be sure we don't have values
   % for a negative ice thickness.
   
   H = max(0 , (H + (dHdt*delt)));
   H(end)=0;
    
 end  % done with each time step
 
 Hinit = H;
 save Hinit.mat Hinit

%%
figure;        

subplot(221);  plot(x,Qm,'co'); ylabel('flux');title('western cv'); 
subplot(222);  plot(x,Qp,'co'); title('eastern cv');
subplot(223);  plot(x,um,'co'); ylabel('velocity'); xlabel('distance (m)'); 
subplot(224);  plot(x,up,'co'); xlabel('distance (m)'); 


figure;
plot(x,H,'c*');

figure;
plot(x,H + zb,'c*'); hold on; plot(x,zb,'k','linewidth',1)

toc



