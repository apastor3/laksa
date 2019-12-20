function DF = Fun_ODE_KE(t,X)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute Right-Hand Side of the equations of motion             %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:                                                                  %
%           t  -> Time                                                       %
%           X  -> Extended state vector                                      %
%   Outputs:                                                                 %
%           DF -> Time derivative of the extended state vector               % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global PND

[R_KE rK vK eulerK omegaK F_Tether M_Tether FA MA alfa beta  rM DF T_U T_D] = Fun_ODE_Full_Output_KE(t,X,PND);


end