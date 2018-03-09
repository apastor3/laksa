function [Ups_s Ups_c Ups_s_xs Ups_s_xc Ups_c_xs Ups_c_xc] = Fun_Matrix_Upsilon_KA(xs,xc,PND)
         
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Matrix Upsilon and its gradients                               %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  xs                 -> state vector                      %%
%          xc                 -> control vector                    %%
%          PND                -> dimensionless parameters          %%       
% Outputs: Ups_s, Ups_c       -> matrices for the computation of   %%
%          the components in SB of the dimesionless velocity of the%%
%          center of mass of the kite with respect to SE:          %%
%          vg = Ups_s*xs_p + Ups_c*xc_p.                           %%
%                                                                  %%
%          Ups_s_xs, Ups_s_xc -> tensors with the partial          %%
%          derivatives of Ups_s with respect to to  the components %%
%          of xs and xc.                                           %%
%                                                                  %%
%          Ups_c_xs, Ups_c_xc -> tensors with the partial          %%
%          derivatives of Ups_c with respect to to  the components %%
%          of xs and xc.                                           %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   


% Recover dimensionless Parameters
xA = PND.Tether.XA;     % Attachment Point Coordinates
zA = PND.Tether.ZA;     % Attachment Point Coordinates

% Recover state vector variables
varphi = xs(1,1);
gamma  = xs(2,1);
eta    = xs(3,1);
theta  = xs(4,1);
% Recover control vector variables
l      = xc(1,1);
delta  = xc(2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Compute matrices Ups_s and Ups_c                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize Matrices
Ups_s = zeros(3,4);
Ups_c = zeros(3,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upsilon_S 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ups_s(1,1) =   l*(sin(delta)*sin(gamma)*sin(theta)-cos(gamma)*sin(eta+delta)*cos(theta))-zA*cos(gamma)*sin(eta);
Ups_s(2,1) =  sin(gamma)*(xA*sin(theta)-zA*cos(theta)-l*cos(delta))-cos(gamma)*cos(eta)*(xA*cos(theta)+zA*sin(theta));
Ups_s(3,1) =  -l*(sin(delta)*sin(gamma)*cos(theta)+cos(gamma)*sin(eta+delta)*sin(theta))+xA*cos(gamma)*sin(eta);

Ups_s(1,2) = -l*cos(eta+delta)*cos(theta)-zA*cos(eta);
Ups_s(2,2) =  (xA*cos(theta)+zA*sin(theta))*sin(eta);
Ups_s(3,2) = -l*cos(eta+delta)*sin(theta)+xA*cos(eta);

Ups_s(1,3) = -l*sin(delta)*sin(theta);
Ups_s(2,3) =  l*cos(delta)-xA*sin(theta)+zA*cos(theta); 
Ups_s(3,3) =  l*sin(delta)*cos(theta);

Ups_s(1,4) = -zA;
Ups_s(2,4) = 0;
Ups_s(3,4) = xA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient of Upsilon_S with respect to the state vector %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ups_s_xs  = zeros(3,4,4); % partial Ups_s/partial x_s
%% Derivative with respect to varphi is equal to zero
%% Derivative with respect to gammma is:
Ups_s_xs(1,1,2) =   l*(sin(delta)*cos(gamma)*sin(theta)+sin(gamma)*sin(eta+delta)*cos(theta))+zA*sin(gamma)*sin(eta);
Ups_s_xs(2,1,2) =  cos(gamma)*(xA*sin(theta)-zA*cos(theta)-l*cos(delta))+sin(gamma)*cos(eta)*(xA*cos(theta)+zA*sin(theta));
Ups_s_xs(3,1,2) =  -l*(sin(delta)*cos(gamma)*cos(theta)-sin(gamma)*sin(eta+delta)*sin(theta))-xA*sin(gamma)*sin(eta);
%% Derivative with respect to eta is :
Ups_s_xs(1,1,3) =  -l*cos(gamma)*cos(eta+delta)*cos(theta)-zA*cos(gamma)*cos(eta);
Ups_s_xs(2,1,3) =  cos(gamma)*sin(eta)*(xA*cos(theta)+zA*sin(theta));
Ups_s_xs(3,1,3) =  -l*cos(gamma)*cos(eta+delta)*sin(theta)+xA*cos(gamma)*cos(eta);

Ups_s_xs(1,2,3) =  l*sin(eta+delta)*cos(theta)+zA*sin(eta);
Ups_s_xs(2,2,3) =  (xA*cos(theta)+zA*sin(theta))*cos(eta);
Ups_s_xs(3,2,3) =  l*sin(eta+delta)*sin(theta)-xA*sin(eta);
%% Derivative with respect to theta is :
Ups_s_xs(1,1,4) =   l*(sin(delta)*sin(gamma)*cos(theta)+cos(gamma)*sin(eta+delta)*sin(theta));
Ups_s_xs(2,1,4) =  sin(gamma)*(xA*cos(theta)+zA*sin(theta))-cos(gamma)*cos(eta)*(zA*cos(theta)-xA*sin(theta));
Ups_s_xs(3,1,4) =   l*(sin(delta)*sin(gamma)*sin(theta)-cos(gamma)*sin(eta+delta)*cos(theta));

Ups_s_xs(1,2,4) =  l*cos(eta+delta)*sin(theta);
Ups_s_xs(2,2,4) = -(xA*sin(theta)-zA*cos(theta))*sin(eta);
Ups_s_xs(3,2,4) = -l*cos(eta+delta)*cos(theta);

Ups_s_xs(1,3,4) =  -l*sin(delta)*cos(theta);
Ups_s_xs(2,3,4) =  -xA*cos(theta)-zA*sin(theta); 
Ups_s_xs(3,3,4) =  -l*sin(delta)*sin(theta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient of Upsilon_S with respect to the control vector %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ups_s_xc  = zeros(3,4,2); % partial Ups_s/partial x_s
%% Derivative with respect to l is:
Ups_s_xc(1,1,1) =   sin(delta)*sin(gamma)*sin(theta)-cos(gamma)*sin(eta+delta)*cos(theta);
Ups_s_xc(2,1,1) =  -sin(gamma)*cos(delta);
Ups_s_xc(3,1,1) =  -(sin(delta)*sin(gamma)*cos(theta)+cos(gamma)*sin(eta+delta)*sin(theta));

Ups_s_xc(1,2,1) = -cos(eta+delta)*cos(theta);
Ups_s_xc(3,2,1) = -cos(eta+delta)*sin(theta);

Ups_s_xc(1,3,1) = -sin(delta)*sin(theta);
Ups_s_xc(2,3,1) =  cos(delta); 
Ups_s_xc(3,3,1) =  sin(delta)*cos(theta);

%% Derivative with respect to delta is:
Ups_s_xc(1,1,2) =   l*(cos(delta)*sin(gamma)*sin(theta)-cos(gamma)*cos(eta+delta)*cos(theta));
Ups_s_xc(2,1,2) =   l*sin(gamma)*sin(delta);
Ups_s_xc(3,1,2) =  -l*(cos(delta)*sin(gamma)*cos(theta)+cos(gamma)*cos(eta+delta)*sin(theta));

Ups_s_xc(1,2,2) =   l*sin(eta+delta)*cos(theta);
Ups_s_xc(3,2,2) =   l*sin(eta+delta)*sin(theta);

Ups_s_xc(1,3,2) =  -l*cos(delta)*sin(theta);
Ups_s_xc(2,3,2) =  -l*sin(delta); 
Ups_s_xc(3,3,2) =   l*cos(delta)*cos(theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upsilon_C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ups_c(1,1) =  cos(delta)*sin(theta);        Ups_c(1,2) = -l*sin(delta)*sin(theta);
Ups_c(2,1) =  sin(delta);                   Ups_c(2,2) =  l*cos(delta);
Ups_c(3,1) = -cos(delta)*cos(theta);        Ups_c(3,2) =  l*sin(delta)*cos(theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient of Upsilon_C with respect to the state vector %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ups_c_xs  = zeros(3,2,4); % partial Ups_c/partial x_s
%% Derivative with respect to varphi is equal to zero
%% Derivative with respect to gamma is equal to zero
%% Derivative with respect to eta   is equal to zero
%% Derivative with respect to theta is:
Ups_c_xs(1,1,4) =  cos(delta)*cos(theta);      Ups_c_xs(1,2,4) =  -l*sin(delta)*cos(theta);
Ups_c_xs(3,1,4) =  cos(delta)*sin(theta);      Ups_c_xs(3,2,4) =  -l*sin(delta)*sin(theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient of Upsilon_C with respect to the control vector %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ups_c_xc  = zeros(3,2,2); % partial Ups_c/partial x_c
% Derivative with respect to l
Ups_c_xc(1,2,1) = -sin(delta)*sin(theta);
Ups_c_xc(2,2,1) =  cos(delta);
Ups_c_xc(3,2,1) =  sin(delta)*cos(theta);
% Derivative with respect to delta
Ups_c_xc(1,1,2) = -sin(delta)*sin(theta);        Ups_c_xc(1,2,2) = -l*cos(delta)*sin(theta);
Ups_c_xc(2,1,2) =  cos(delta);                   Ups_c_xc(2,2,2) = -l*sin(delta);
Ups_c_xc(3,1,2) =  sin(delta)*cos(theta);        Ups_c_xc(3,2,2) =  l*cos(delta)*cos(theta);


end