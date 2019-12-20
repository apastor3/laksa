function DF =  Fun_ODE_Full_KT(t,xs_amp)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and  Jose A. Serrano-Iglesia           %
% Language  : Matlab                                                         %
% Synopsis  : RHS with Full Outputs                                          %
% Copyright : Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%     t                       -> normalized time                          %%
%     Extended state vector   -> xs_amp = [xs xs_p]                       %%
% Outputs                                                                 %%
%     DF     -> Righthand side of the first order Equations               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global PND

[xc RBE rK vK omegaK FA MA alfa beta Ups Ups_xs Phi Phi_xs Q DF RHS]=  Fun_ODE_Full_Output_KT(t,xs_amp,PND);





end