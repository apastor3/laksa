function [U U_xs] = Fun_Potential_KA(xs,xc,PND)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute potential and its gradient                             %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  xs    -> state vector                                   %%
%          xc    -> control vector                                 %%
%          PND   -> dimensionless parameters of the system         %%
%                                                                  %%
% Outputs: U    -> Dimensionless potential                         %%
%          U_xs -> partial derivatives of U with respect to        %%
%                  the state vector components                     %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xA = PND.Tether.XA ;     % XA Attachment Point Coordinates
zA = PND.Tether.ZA ;     % ZA Attachment Point Coordinates

% State Vector
varphi = xs(1,1);
gamma  = xs(2,1);
eta    = xs(3,1);
theta  = xs(4,1);
% Control Vector
l      = xc(1,1);
delta  = xc(2,1);
% Potential
U = l*cos(gamma)*cos(eta+delta)-xA*(sin(gamma)*cos(theta)+cos(gamma)*sin(theta)*cos(eta))-zA*(sin(gamma)*sin(theta)-cos(gamma)*cos(theta)*cos(eta));
%%%%%%%%%%%%%%%%%%%%%%%%
% Potential Gradient  %%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Derivative with respect to varphi
U_xs(1,1) = 0;
%% Derivative with respect to gamma
U_xs(2,1) = -l*sin(gamma)*cos(eta+delta)-xA*(cos(gamma)*cos(theta)-sin(gamma)*sin(theta)*cos(eta))-zA*(cos(gamma)*sin(theta)+sin(gamma)*cos(theta)*cos(eta));
%% Derivative with respect to eta
U_xs(3,1) = -l*cos(gamma)*sin(eta+delta)+xA*cos(gamma)*sin(theta)*sin(eta)-zA*cos(gamma)*cos(theta)*sin(eta);
%% Derivative with respectto theta
U_xs(4,1) = xA*(sin(gamma)*sin(theta)-cos(gamma)*cos(theta)*cos(eta))-zA*(sin(gamma)*cos(theta)+cos(gamma)*sin(theta)*cos(eta));
 

end