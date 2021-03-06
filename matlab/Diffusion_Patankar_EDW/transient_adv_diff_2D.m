%%  ESS 590D: HEAT AND MASS FLOW MODELLING  SPRING 2006
%%
%%  Finite-Volume solution for transient phi(x,z)
%%  diffusion with advection
%%
%%    d(rho phi)/dt + div( rho u phi ) = div( Gamma grad phi) + S 
%%
%% For example,
%%  phi   = temperature                         (deg)
%%  Gamma = k/c
%%             k = thermal conductivity         (W m^-1 deg^-1)
%%             c = specific heat capacity       (J kg^-1 deg^-1)
%%  S     = sources                             (W m^-3)
%%        + S_C + S_P phi
%%
%%  This code uses PATANKAR POWER-LAW SCHEME for advective fluxes
%%
%%  Externally pecified velocity field (u,w) must astisfy 
%%   mass conservation
%%
%%   d(rho)/dt + div( rho u) = 0
%%
%%  Boundary Conditions:
%%   All boundaries are finite-volume interfaces
%%   (Patankar Practice B, page 69)
%%   However, there is an additional "point" X_B located on boundary.
%%   Its value and coefficients are determined by phi value or gradient
%%   boundary conditions there.
%%   Adjacent interior point is X_I
%%   (X_I - X_B) = dX_B
%%
%%    phi_B = phi_0
%%        a_B = 1,  all other a_J = 0
%%        b   = phi_0
%%
%%    d phi/dx = dphi_dx_0
%%        use gradient between X_B and X_I (interior point)
%%        d phi/dx|B = Gamma_I (phi_I - phi_B)/dx_B
%%    or
%%        phi_B = dphi_dx_0 * dx_B + phi_I
%%
%%
%%  vector bc_u = phi_u       specified value of phi on upper boundary
%%  vector bc_d = phi_d       specified value of phi on lower boundary
%%  vector bc_w = dphi_dx_w   specified gradient on west boundary
%%  vector bc_e = dphi_dx_e   specified gradient on east boundary
%%
%%  THIS CODE APPLIES GRADIENT BOUNDARY CONDITIONS ON VOLUME INTERFACES
%%  on East and West
%%  APPLIES FIXED-VALUE BOUNDARY CONDITIONS ON NODES
%%  in Up and Down directions
%%
%% Material Properties
%%  K can vary with position - K(x,z)
%%
%%
%%  In the spirit of matlab's vector capabilities, the program
%%   avoids loops over spatial indices.
%%  Instead, matrices are formed for all spatially varying parameters
%%   and variables.
%%
%%                                   Ed Waddington      May 9 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   global c_x  c_z  c_prime
%
% set strain rate for velocity calculations
%   u = c_x + c_prime * x
%   w = c_z - c_1 * z
   c_x = -1;
   c_z = 0.5;
   c_prime = 0;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%  set up finite volumes
%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This code sets up uniform-sized volumes (could be generalized)
%   There are points on all boundaries, which are also
%   volume interfaces  (Patankar Practice B, page 69)
%
   Lx    = 10;               %  x-extent of domain
   nx_fv = 10;                %  number of finite volumes in x
   nx_P  = nx_fv + 2;        %  number of nodes in x
%
%
   Lz    = 1;               %  z-extent of domain
   nz_fv = 1;                %  number of finite volumes in z
   nz_P  = nz_fv + 2;        %  number of nodes in z
%
%
%  set up matrices X_edges, Z_edges - defining (x,z) coordinates 
%  of finite-volume edges
%
%    set x-widths of volumes 
%    this bit of code sets up uniform edge spacing
%   ----------------------------------------------
     dx = Lx/nx_fv;              %  z-widths of volumes
     x_edges_vec = dx * (0:nx_fv);
%
  %%   X_edges = repmat(x_edges_vec, nz_fv+1, 1);
%
     dz = Lz/nz_fv;               %  z-widths of volumes
     z_edges_vec = dz * (0:nz_fv)';
%
     Z_edges = repmat(z_edges_vec, 1, nx_fv+1);

%
% set nodes at finite-volume centers and on all boundaries
     x_P_vec = [ x_edges_vec(1) ...
                (x_edges_vec(1:end-1) + diff(x_edges_vec)/2) ...
                 x_edges_vec(end) ]; 
%
  %%   X_P = repmat(x_P_vec, nz_P, 1);
%
%
     z_P_vec = [ z_edges_vec(1) ...
                (z_edges_vec(1:end-1) + diff(z_edges_vec)/2)' ...
                 z_edges_vec(end) ]'; 
%
   %%  Z_P = repmat(z_P_vec, 1, nx_P);
%
%
%
%%  Boundary Conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%  Boundary conditions along each boundary are carried in an array
%  labelled bc_x.  For each boundary point, 
%        first entry is numerical value of BC,
%        second entry indicates type of BC
%               1 for specified value
%               2 for specified gradient
%  ---------------------------------------------------
%
%  phi_e at east boundary
     bc_e_0 = 0;
     bc_type   = 2;
     bc_e = repmat( [ bc_e_0  bc_type ], nz_P, 1 );
%
%%  phi_w at west boundary
     bc_w_0 = 1;
     bc_type   = 0;
     bc_w = repmat( [ bc_w_0  bc_type ], nz_P, 1 );
%
%
%%  phi gradient on upper boundary
     bc_u_0 = 0;
     bc_type = 2;
     bc_u   = repmat( [ bc_u_0  bc_type ]', 1, nx_P );
%
%
%%  phi gradient on lower boundary
     bc_d_0 = 0;
     bc_type = 2;
     bc_d   = repmat( [ bc_d_0  bc_type ]', 1, nx_P );
%

%
%
%  set time-step control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      dt = 1;
      nt = 25;
%
%  set up initial condition
%  ------------------------
      phi_0 = zeros( nz_P, nx_P );
%
%
%%  material parameters
%%%%%%%%%%%%%%%%%%%%%%%
% set up thermal-property matrix - uniform conductivity
     Gamma_0 = -1e6;                     %   W m^-1 deg^-1    
     Gamma_P = repmat( Gamma_0, nz_P, nx_P );
%
%
% set up density matrix
     rho_0 = 1;                       %   kg m^-3
     rho = repmat( rho_0, nz_P, nx_P );
%

%
%%  source term   
%%%%%%%%%%%%%%%%%%%%%
%   S = S_C + S_P * phi_P 
%   ----------------------
     S_C_0  = 0;                   %  W m^-3 s^-1
     S_P_0  = 0;                   %  W m^-3 s^-1 deg^-1
%
     S_C  = repmat( S_C_0, nz_P, nx_P );  %  constant source term
     S_P  = repmat( S_P_0, nz_P, nx_P );  %  phi_P-dependent source
%
%

%
% This is it! solve for phi using Finite Volumes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phi = fvm_transient_adv_diff_2D( x_edges_vec, z_edges_vec, ...
                  x_P_vec, z_P_vec, nt, dt, Gamma_P, rho, S_C, S_P, ...
                  bc_u, bc_d, bc_e, bc_w, phi_0 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%
%%  Plot results
%%%%%%%%%%%%%%%%%

%  recover node positions (volume centers and BC points)
%  -----------------------------------------------------
     X_P = repmat(x_P_vec, nz_P, 1);
     Z_P = repmat(z_P_vec, 1, nx_P);
%

%       
%  plot phi(x,z) solution as perspective surface
%  ----------------------------------------------
  %%  nt = size(phi,2)/nx_P;
  %%  j_stop = 0;
  %%  for i_time = 1:nt
  %%      j_start = j_stop + 1;
  %%      j_stop  = j_start -1 + nx_P; 
  %%      phi_plot = phi(:,j_start:j_stop);
  %%      surf(X_P, Z_P, phi_plot )
  %%      colorbar
  %%      grid on; box on
  %%      xlabel('Distance  (m)')
  %%      ylabel('Depth  (m)')
  %%      zlabel('Temperature  (deg C)')
  %%      pause
  %%  end
    
%
%
%
%
%

