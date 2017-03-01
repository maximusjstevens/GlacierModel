  function w = w( xpos, zpos ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finds the vertical velocity w 
%% at positions x = xpos, z = zpos 
%                  
%  
%----------------------------------------
 
    global c_x  c_z  c_prime



    w = c_z - c_prime * zpos;
%
%
%
