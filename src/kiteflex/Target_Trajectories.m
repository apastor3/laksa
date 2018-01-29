function [xs_target xs_p_target xs_pp_target ]   = Target_Trajectories(t,PND)
         

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute Figure-Of-Eight trajectory                             %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  s      -> arc-length parameter (0<s<2pi)                %%
%          L_CM   -> distance fro SE to the center of mass         %%
%          PND_CL -> dimensionless parameters of the close loop    %%
%                                                                  %%
% Outputs: rE   -> Vector position (SE components)                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators
NT = 4+NG+3;       % Contrl variables

rE        = zeros(3,1);
xs_target = zeros(2*NR+3,1);


if PND.CL.Target.Type == 0 % Equilibrium State
 
   Theta0       =  0.137906823298247;

   xs_target    = [zeros(2*NR,1); Theta0; 0; 0];
   xs_p_target  = [zeros(2*NR+3,1)];
   xs_pp_target = [zeros(2*NR+3,1)];
    
end


 
end