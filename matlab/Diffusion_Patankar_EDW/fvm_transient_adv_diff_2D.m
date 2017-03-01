function phi = fvm_transient_adv_diff_2D( x_edges_vec, z_edges_vec, ...
                     x_P_vec, z_P_vec, nt, dt, Gamma_P, rho, S_C, S_P, ...
                      bc_u, bc_d, bc_e, bc_w, phi_0 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  find transient phi(x,z), for combined advection and diffusion
%%  with given heat sources S and conductivity Gamma_P.
%%
%%    d(rho phi)/dt + div( rho u phi) = div( K grad(phi) ) + S
%%
%%  Uses Patankar Finite Volume method.
%%  Advective flux included through PATANKAR POWER-LAW SCHEME
%%
%%
%%  EXTERNALLY SPECIFIED VELOCITY MUST SATISFY CONTINUITY
%%   i.e. 
%%     div(rho u) = 0
%%
%%  Boundary Conditions:
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
%%  For each boundary there is a matrix bc_x
%%     with either 2 rows (Up,Down) or 2 columns (East,West).
%%  e.g. East  
%%   bc_e is [nz_P x 2] where nz_P is number of nodes on boundary 
%%     - first column contains value of BC at that node on boundary
%%     - second column specifies type of BC
%%            - 1 -> specified phi
%%            - 2 -> specified phi gradient across boundary
%%
%%   Similarly bc_w, bc_u, bc_d on other boundaries
%%
%%
%%  a_P = coefficient of phi_P at point P = (i,j)
%%  a_E = coefficient of phi_E at point east of (i,j)
%%  a_E, a_W, a_U, a_D for points East, West, Up, Down
%%
%%  System of equations to be solved for phi
%%
%%       big_A * phi = rhs
%%
%%  big_A = matrix with coefficients a_E etc, modified to include BC
%%  phi   = unknown values in Control-volumes
%%  rhs   = right-hand-side column vector containing sources and BC terms
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   global c_x  c_z  c_prime
%
% get number of Finite-Volume nodes
% ----------------------------------
    nx_P = length(x_P_vec);
    nz_P = length(z_P_vec);
%
% get number of finite-volume faces
     nz_fv = nz_P - 2;
     nx_fv = nx_P - 2;
%
% set up arrays defining nodes locations
     X_P = repmat(x_P_vec, nz_P, 1);
     Z_P = repmat(z_P_vec, 1, nx_P);


%  get dX, dZ - spacings between volume edges
%  -------------------------------------------
   %%    dx_vec = [ 0  diff(x_edges_vec)  0];
   %%    dz_vec = [ 0  diff(z_edges_vec)'  0]';
   %% This leads to division by zero for Peclet numbers
       dx_vec = [ 1  diff(x_edges_vec)  1 ];
       dz_vec = [ 1  diff(z_edges_vec)'  1 ]';
%
      dX = repmat( dx_vec, nz_P, 1 );
      dZ = repmat( dz_vec, 1, nx_P );
%

%
%  get dX_E, dx_W, dZ_U, dZ_D - 
%          - distances to adjacent volume centers
%  -----------------------------------------------
%
       dX_e = diff(X_P')';
    % add dummy row at eastern edge to maintain size [nx_P, nz_P]
       dX_e = [ dX_e dX_e(:,end) ];
%
       dX_w = diff(X_P')';
    % add dummy row at western edge to maintain size [nx_P, nz_P]
       dX_w = [ dX_w(:,1) dX_w ];
%
       dZ_u = diff(Z_P);
    % add dummy row at top edge to maintain size [nx_P, nz_P]
       dZ_u = [ dZ_u(1,:)' dZ_u' ]';
%
       dZ_d = diff(Z_P);
    %  add dummy row at bottom edge to maintain size [nx_P, nz_P]
       dZ_d = [ dZ_d' dZ_d(end,:)' ]';

%

%  get f_e, f_w, f_u, f_d 
%          - fractional distances to interface in adjacent volume
%  --------------------------------------------------------------
       f_e = [(1 -(x_edges_vec -x_P_vec(1:end-1))./dX_e(1,1:end-1)) 0];
       f_e = repmat( f_e, nz_P, 1 );
       
       f_w = [0 (1-(x_P_vec(2:end) -x_edges_vec )./dX_w(1,2:end) )];
       f_w = repmat( f_w, nz_P, 1 );
       
       f_u = [0 (1 -(z_P_vec(2:end) -z_edges_vec)./dZ_u(2:end,1) )']';
       f_u = repmat( f_u, 1, nx_P );
%
     f_d = [(1 -(z_edges_vec -z_P_vec(1:end-1) )./dZ_d(1:end-1,1))' 0]';
     f_d = repmat( f_d, 1, nx_P );
%



% set up b matrix - includes constant part of source matrix S 
%         S = S_C + S_P * phi_P    (units are W m^-3)
% ----------------------------------------------------------
    b_0 = S_C.* dX.* dZ;
%
%
% set up matrices for conductivities in adjacent volumes
% ------------------------------------------------------
    Gamma_E   =  [ Gamma_P(:,2:end) Gamma_P(:, end) ];
    Gamma_W   =  [ Gamma_P(:,1) Gamma_P(:,1:end-1) ];
    Gamma_U   =  [ Gamma_P(1,:)' Gamma_P(1:end-1,:)' ]';
    Gamma_D   =  [ Gamma_P(2:end,:)' Gamma_P(end,:)' ]';
%
%
% set up interface conductivity matrices
% --------------------------------------
    Gamma_e =  1./ ( (1 - f_e)./Gamma_P + f_e./Gamma_E );
    Gamma_w =  1./ ( (1 - f_w)./Gamma_P + f_w./Gamma_W );
    Gamma_u =  1./ ( (1 - f_u)./Gamma_P + f_u./Gamma_U );
    Gamma_d =  1./ ( (1 - f_d)./Gamma_P + f_d./Gamma_D );
%
%
% Advection
%%%%%%%%%%%%%%%%%
%  set up positions X and Z where velocities u and w
%  are calculated 
%  - u on East and West interfaces
%  - w on Up and Down interfaces
     X_u_velo = repmat( x_edges_vec, nz_P, 1 );
     Z_u_velo = repmat( z_P_vec, 1, nx_P-1 );
     X_w_velo = repmat( x_P_vec, nz_P-1, 1 );
     Z_w_velo = repmat( z_edges_vec, 1, nx_P );
%
%  find the velocities
     u_edges = u(X_u_velo, Z_u_velo);
     w_edges = w(X_w_velo, Z_w_velo);
%
%  split them out into East and West, Up and Down
     u_e = [ u_edges  u_edges(:, end) ];
     u_w = [ u_edges(:, 1)  u_edges ];
%
      w_u = [ w_edges(1,:)'  w_edges' ]';
      w_d = [ w_edges'  w_edges(end,:)' ]';
      
   %%   clear u_edges w_edges
%
%
% set up Flow and Diffusion coefficients
% ---------------------------------------
    D_e = (Gamma_e ./ dX_e).* dZ;
    D_w = (Gamma_w ./ dX_w).* dZ; 
    D_u = (Gamma_u ./ dZ_u).* dX;
    D_d = (Gamma_d ./ dZ_d).* dX;
%
    F_e = ( rho.* u_e ).* dZ;
    F_w = ( rho.* u_w ).* dZ; 
    F_u = ( rho.* w_u ).* dX;
    F_d = ( rho.* w_d ).* dX;
%

% set up Peclet numbers
% ----------------------
    P_e = F_e./ D_e;
    P_w = F_w./ D_w;
    P_u = F_u./ D_u;
    P_d = F_d./ D_d;

%
% set up coefficients at adjacent points
% ---------------------------------------
     a_E = D_e .* A( P_e ) + F_upwind( -F_e );
     a_W = D_w .* A( P_w ) + F_upwind(  F_w );
     a_U = D_u .* A( P_u ) + F_upwind(  F_u );
     a_D = D_d .* A( P_d ) + F_upwind( -F_d );
%
     a_P_0 = rho.* dX.*dZ./ dt;
%
%         
% set up coefficient at central point P (incl source dependent on phi)
% --------------------------------------------------------------------
    a_P = a_E + a_W + a_U + a_D + a_P_0 - S_P.* dX.* dZ; 
%
%
%       
%  plot phi_0 initial condition
%  -----------------------------
        figure(1)
        surf(X_P, Z_P, phi_0 )
        colorbar
        grid on; box on
        xlabel('Distance  (m)')
        ylabel('Depth  (m)')
        zlabel('Temperature  (deg C)')
        pause(1)

%

%%%%%%%%%%%%%%%%%%%%%%%%%
%   time step loop here
%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Boundary conditions are reset each time step because b 
%  depends on phi_t (previous solution)
%  While this may be redundant for modifications to the
%  coefficients a_E etc, it may be conceptually clearer to 
%  set the coefficients at the boundary at the same time
%  as we set the values for b, which _does_ depend on time
%  for transient boundary conditions, and when old solution
%  _phi_t changes with time 
%
%
%  phi stores _all_ solutions through time.
%  This will get cumbersome with very many time steps!
%  You may want to replace this with a different way 
%  to store results
% ----------------------------------------------------------
%
%
%  put initial condition phi_0 into big storage matrix phi
    phi = [ phi_0 ];
%
%  put initial condition into phi_t, which carries solution 
%  forward from previous time step
    phi_t = phi_0;

%
    for i_time = 1:nt
%
%   include old solution (initial condition) in rhs vector
        b = b_0 + a_P_0 .* phi_t;
%
%
  %  set Boundary Conditions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %   specified values phi_P = phi_value at x_P
  %     set a_P = 1, all other a_X = 0
  %           b = phi_value = bc_x(index, 1);
%
  %  specified gradients
  %     gradient = (phi_I - phi_P)/dX_I  where dX_I = (X_I - X_P)
  %  or
  %     phi_P = phi_I + dX_I * gradient 
  %  so
  %    a_P = 1
  %    a_I = 1
  %    b   = dX_I * gradient
  %  -------------------------------------
%
  %  East Boundary
  %  --------------
       bc_type = bc_e(:,2);
%
       a_P(:,end) = 1; 
       a_U(:,end) = 0;
       a_D(:,end) = 0;
       a_W(:,end) = 0;
       a_E(:,end) = 0;
%
   % index identifies boundary nodes where phi is specified
       index = find(bc_type == 1);
       if( ~isempty(index) )
          b(index, end) = bc_e(index, 1);
       end
%
   % index identifies boundary nodes where gradient is specified
       index = find( bc_type == 2 );
       if( ~isempty(index) )
          a_W(index, end) = 1;
          b(index, end)   = dX_e(index, end).* bc_e(index, 1);
       end
%
%
  %  West Boundary
  %  --------------
       bc_type = bc_w(:,2);

       a_P(:,1) = 1; 
       a_U(:,1) = 0;
       a_D(:,1) = 0;
       a_W(:,1) = 0;
       a_E(:,1) = 0;
%
    % index identifies boundary nodes where phi is specified
       index = find( bc_type == 1 );
       if( ~isempty(index) )
          b(index, 1) = bc_w(index, 1);
       end
%
    % index identifies boundary nodes where gradient is specified
       index = find( bc_type == 2 );
       if( ~isempty(index) )
          a_E(index, 1) = 1;
          b(index, 1)   = - dX_w(index,1).* bc_w(index, 1);
       end
%
%
%
  %  Up Boundary
  %  --------------
       bc_type = bc_u(2,:);

       a_P(1,:) = 1; 
       a_U(1,:) = 0;
       a_D(1,:) = 0;
       a_W(1,:) = 0;
       a_E(1,:) = 0;
%
    % index identifies boundary nodes where phi is specified
       index = find( bc_type == 1 );
       if( ~isempty(index) )
          b(1, index) = bc_u(1, index);
       end
%
    % index identifies boundary nodes where gradient is specified
       index = find( bc_type == 2 );
       if( ~isempty(index) )
          a_D(1, index) = 1;
          b(1, index)   = - dZ_d(1, index).* bc_u(1, index);
       end
%
%
  %  Down Boundary
  %  --------------
       bc_type = bc_d(2,:);

       a_P(end,:) = 1; 
       a_U(end,:) = 0;
       a_D(end,:) = 0;
       a_W(end,:) = 0;
       a_E(end,:) = 0;
%
    % index identifies boundary nodes where phi is specified
       index = find( bc_type == 1 );
       if( ~isempty(index) )
          b(end, index) = bc_d(1, index);
       end
%
    % index identifies boundary nodes where gradient is specified
       index = find(bc_type == 2);
       if( ~isempty(index) )
          a_U(end, index) = 1;
          b(end, index)   = dZ_u(end, index).* bc_d(1, index);
       end
%%
%
%
%
%
   %   Now call matrix solver
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        phi_t = solver( a_E, a_W, a_U, a_D, a_P, b );
%
   %   Stash solution for plotting later
        phi = [ phi phi_t ];
%

%       
  %  plot phi(x,z) solution as perspective surface
  %  ----------------------------------------------
        figure(1)
        surf(X_P, Z_P, phi_t )
        colorbar
        grid on; box on
        xlabel('Distance  (m)')
        ylabel('Depth  (m)')
        zlabel('Temperature  (deg C)')
   %  set limits differently for different problems
        set(gca,'zlim',[0,50] ) 
        pause(.1)

    end
%
%
