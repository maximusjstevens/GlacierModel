  function u = u( xpos, zpos ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finds the horizontal velocity u 
%% at positions x = xpos, z = zpos
%%
%%
%%---------------------------------------

   global c_x  c_z  c_prime
   
    
    u =  c_x + c_prime * xpos;
%
%
%
%
