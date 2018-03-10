function DF = Fun_ODE_KS(t,X)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Time derivative of the extended state vector                   %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%                                                                    %%
% Inputs:  t               -> dimensionless time                     %%
%          X = [xs xs_dot] -> extended state vector                  %%
% Outputs: DF              -> d X/dt                                 %%
%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global PND 

[rg vk ak omega_BE alfa_BE T_Bp T_Bm m_Bp m_Bm fa ma alfa beta ...
          Rp Rm Elong_p Elong_m DF RHS] = Fun_ODE_Full_Output_KS(t,X,PND);


end