   function A = A( P )
% =======================
% calculates the upwinding weight as a function of Peclet number
%   This routine uses Patankar's Power Law scheme (page 95)
%      max[ 0, (1 - 0.1 |P|)^5 ]
%
%
     A = max( (1 - 0.1 * abs( P ) ).^5, zeros( size(P) ) );
%
%
%
