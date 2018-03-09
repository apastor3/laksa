function [Q f m alfa beta] = Fun_Q_KA(vw,vg,omega,Ups_s,Phi,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Generalized forces                                             %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  vw    -> components of the dimensionless wind velocity  %%
%                   vector in the body frame                       %%
%          vg    -> components of the dimensionless velocity vector%% 
%                   of the kite in the body frame                  %%
%          omega -> kite angular velocity components in body frame %% 
%          Ups_s -> see Fun_Matrix_Upsilon                         %%
%          Phi   -> see Fun_Matrix_Omega                           %%
%          PND   -> dimensionless parameters of the system         %%
%                                                                  %%
% Outputs: Q     -> generalized forces                             %%
%          f     -> components of the aerodynamic force of the     %%
%                   kite in the body frame                         %%
%          m     -> components of the aerodynamic torque about the %%
%                   center of mass of the kite  in the body frame  %%
%          alfa  -> angle of attack  (rad)                         %% 
%          beta  -> sideslip angle   (rad)                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Components of the aerodinamic velocity vector in the body frame
VA     = vg-vw;
%% Components in the body frame of the aerodynamic force and torque 
[f m alfa beta] = Aerokite(VA,omega,PND,0,0,0);
%% Generalized forces
Q        = f'*Ups_s+m'*Phi;
Q        = Q';

end