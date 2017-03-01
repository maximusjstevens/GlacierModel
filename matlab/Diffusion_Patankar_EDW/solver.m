      function phi = solver( a_E, a_W, a_U, a_D, a_P, b )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  matrix solution for phi_P at points P_ij
%%  a_P = coefficient of phi_P at point P_ij
%%  a_E = coefficient of phi_E  at point east of P_ij,   i.e.  (i,j+1)
%%  a_W = coefficient of phi_W  at point west of P_ij,   i.e.  (i,j-1)
%%  a_U = coefficient of phi_U  at point up from P_ij,   i.e.  (i-1,j)
%%  a_D = coefficient of phi_D  at point down from P_ij, i.e.  (i+1,j)
%%
%%  system of equations to be solved for phi
%%
%%       big_A * phi = rhs
%%
%%  big_A = sparse matrix with diagonals -a_P, a_E etc, 
%%  phi   = unknown solution at grid points P
%%  rhs   = column vector of constant terms  derived from b
%             (includes constant source, BC, previous step, etc)
%%
%% On entry:
%%  a_E, a_W, a_U, a_D, a_P and rhs_stuff are matrices with  
%%  the same dimensions as solution matrix.
%%  They have already been modified to include Boundary Conditions
%%  Each entry X_ij in one of these matrices "X" contains terms 
%%  associated with control-volume center point P_ij.
%%
%%  a_P phi_P = a_E phi_E + a_W phi_W + a_U phi_U + a_D phi_D + b
%%
%% or:
%%
%%  - a_P phi_P + a_E phi_E + a_W phi_W + a_U phi_U + a_D phi_D = - b
%%
%%
%%  Each matrix a_E, a_W, a_U, a_D, a_P is reshaped into a
%%  column vector, which becomes a diagonal of the
%%  sparse matrix big_A
%%    a_P         -> main diagonal
%%    a_U and a_D -> adjacent diagonals
%%    a_E and a_W -> diagonals nz+1 above and below main diagonal 
%%
%%  matrix b contains rhs stuff 
%%     - source contributions S_C independent of phi
%%     - B.C. contributions
%%     - values of solution at previous step (if transient problem)
%%  b is reshaped into a column vector rhs
%%         rhs  = right-hand-side vector
%%
%%-----------------------------------------------------------------
%  find nz, nx (number of rows and columns)
         [nz, nx] = size( b );
         n_volumes = nx*nz;
%
% reshape a_X matrices into diagonals of big_A matrix
%   column 3 is center diagonal of big_A
           Diags = zeros( n_volumes, 5 );
       %%    cols = [ -(nz+1) -1  0  1  (nz+1) ];
           cols = [ -nz -1  0  1  nz ];
%
           Diags(:, 5) = reshape(a_W, n_volumes, 1 );
           Diags(:, 4) = reshape(a_U, n_volumes, 1 );
           Diags(:, 3) = reshape(-a_P, n_volumes, 1 ); % note minus sign
           Diags(:, 2) = reshape(a_D, n_volumes, 1 );
           Diags(:, 1) = reshape(a_E, n_volumes, 1 );
%
%   note that spdiags stuffs the transpose of big_A, because
%   all diagonal elements in column 1 refer to grid point (1,1)
%   and I want column 1 to transpose to the row equation for phi(1,1)
           big_A = spdiags( Diags, cols, n_volumes, n_volumes)';
%
%                
% reshape rhs matrix into column vector
         rhs = reshape(-b, n_volumes, 1);   % note minus sign on b
%
%
%        Solve for phi
          phi = big_A \ rhs;
%
%        Reshape column vector into matrix
          phi = reshape( phi, nz, nx );
%    
%
%
