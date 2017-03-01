   function F_upwind = F_upwind( F )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   F_upwind is a copy of F with negative values replaced by zero
%
%
     F_upwind = max( F, 0 );
%
%
%
%
