function [U U_xs] = Fun_Potential_KS(xs,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors    : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,          %
% Language  : Matlab                                                         %
% Synopsis  : Potential energy and its gradient                              %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %%
% Inputs:  xs              -> state vector of the kite               %%
%          xc = [l delta]  -> control vector of the kite             %%
% Outputs: U               -> Dimensionless potential energy         %%
%          U_xs            -> Gradient of U with respect to xs       %%
%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xA = PND.Tether.XA;     % Attachment Point Coordinates
zA = PND.Tether.ZA;     % Attachment Point Coordinates

% State Vector
varphi = xs(1,1);
gamma  = xs(2,1);
eta    = xs(3,1);
theta  = xs(4,1);
chi    = xs(5,1);
% Control Vector
l0     = PND.Tether.l;
l1     = PND.Bar.Ls;
% Potential
U =   cos(gamma)*( l0*cos(eta) + l1*cos(eta+chi) )...
     -xA*(sin(gamma)*cos(theta)+cos(gamma)*sin(theta)*cos(eta))-zA*(sin(gamma)*sin(theta)-cos(gamma)*cos(theta)*cos(eta));
%%%%%%%%%%%%%%%%%%%%%%%%
% Potential Gradient  %%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Derivative with respect to varphi
U_xs(1,1) = 0;
%% Derivative with respect to gamma
U_xs(2,1) = -sin(gamma)*( l0*cos(eta) + l1*cos(eta+chi) )...
            -xA*(cos(gamma)*cos(theta)-sin(gamma)*sin(theta)*cos(eta))-zA*(cos(gamma)*sin(theta)+sin(gamma)*cos(theta)*cos(eta));
%% Derivative with respect to eta
U_xs(3,1) = -cos(gamma)*( l0*sin(eta) + l1*sin(eta+chi) )...
            +xA*cos(gamma)*sin(theta)*sin(eta)-zA*cos(gamma)*cos(theta)*sin(eta);
%% Derivative with respectto theta
U_xs(4,1) = xA*(sin(gamma)*sin(theta)-cos(gamma)*cos(theta)*cos(eta))-zA*(sin(gamma)*cos(theta)+cos(gamma)*sin(theta)*cos(eta));
%% Derivative with respect to chi
U_xs(5,1) = -l1*cos(gamma)*sin(eta+chi);
    


end