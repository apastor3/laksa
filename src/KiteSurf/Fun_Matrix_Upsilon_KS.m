function [Ups_s Ups_s_xs ] = Fun_Matrix_Upsilon_KS(xs,PND)
         
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Velocity matrix and its gradient                               %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %%
% Inputs:  xs                 -> state vector of the kite            %%
%          PND                -> dimensionless parameters            %%       
% Outputs: Ups_s              -> matrices for the computation of     %%
%          the components in SB of the dimesionless velocity of the  %%
%          center of mass of the kite with respect to SE:            %%
%          vg = Ups_s*xs_p                                           %%
%                                                                    %%
%          Ups_s_xs           -> tensors with the partial            %%
%          derivatives of Ups_s with respect to to  the components   %%
%          of xs .                                                   %%
%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% Recover dimensionless Parameters
xA = PND.Tether.XA;     % Attachment Point Coordinates
zA = PND.Tether.ZA;     % Attachment Point Coordinates

% Recover state vector variables
varphi = xs(1,1);
gamma  = xs(2,1);
eta    = xs(3,1);
theta  = xs(4,1);
chi    = xs(5,1);
% Recover control vector variables
l0     = PND.Tether.l;
l1     = PND.Bar.Ls;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Compute matrices Ups_s and Ups_c                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize Matrices
Ups_s0     = zeros(3,5);
Ups_s1     = zeros(3,5);
Ups_s0_xs  = zeros(3,5,5); % partial Ups_s/partial x_s
Ups_s1_xs  = zeros(3,5,5); % partial Ups_s/partial x_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upsilon_S0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ups_s0(1,1) =  -cos(gamma)*sin(eta)*(zA+l0*cos(theta));
Ups_s0(2,1) =  sin(gamma)*(xA*sin(theta)-zA*cos(theta)-l0)-cos(gamma)*cos(eta)*(xA*cos(theta)+zA*sin(theta));
Ups_s0(3,1) =   cos(gamma)*sin(eta)*(xA-l0*sin(theta));

Ups_s0(1,2) =   -cos(eta)*(zA+l0*cos(theta));
Ups_s0(2,2) =  (xA*cos(theta)+zA*sin(theta))*sin(eta);
Ups_s0(3,2) =    cos(eta)*(xA-l0*sin(theta));

Ups_s0(2,3) =  l0-xA*sin(theta)+zA*cos(theta); 

Ups_s0(1,4) = -zA;
Ups_s0(3,4) =  xA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient of Upsilon_S0 with respect to the state vector %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Derivative with respect to varphi is equal to zero
%% Derivative with respect to gammma is:
Ups_s0_xs(1,1,2) =   sin(gamma)*sin(eta)*(zA+l0*cos(theta));
Ups_s0_xs(2,1,2) =   cos(gamma)*(xA*sin(theta)-zA*cos(theta)-l0)+sin(gamma)*cos(eta)*(xA*cos(theta)+zA*sin(theta));
Ups_s0_xs(3,1,2) =  -sin(gamma)*sin(eta)*(xA-l0*sin(theta));
%% Derivative with respect to eta is :
Ups_s0_xs(1,1,3) = -cos(gamma)*cos(eta)*(zA+l0*cos(theta));
Ups_s0_xs(2,1,3) =  cos(gamma)*sin(eta)*(xA*cos(theta)+zA*sin(theta));
Ups_s0_xs(3,1,3) =  cos(gamma)*cos(eta)*(xA-l0*sin(theta));

Ups_s0_xs(1,2,3) =  sin(eta)*(zA+l0*cos(theta));
Ups_s0_xs(2,2,3) =  (xA*cos(theta)+zA*sin(theta))*cos(eta);
Ups_s0_xs(3,2,3) =  -sin(eta)*(xA-l0*sin(theta));
%% Derivative with respect to theta is :
Ups_s0_xs(1,1,4) =    l0*cos(gamma)*sin(eta)*sin(theta);
Ups_s0_xs(2,1,4) =  sin(gamma)*(xA*cos(theta)+zA*sin(theta))-cos(gamma)*cos(eta)*(zA*cos(theta)-xA*sin(theta));
Ups_s0_xs(3,1,4) =   -l0*cos(gamma)*sin(eta)*cos(theta);

Ups_s0_xs(1,2,4) =  l0*cos(eta)*sin(theta);
Ups_s0_xs(2,2,4) = -(xA*sin(theta)-zA*cos(theta))*sin(eta);
Ups_s0_xs(3,2,4) = -l0*cos(eta)*cos(theta);

Ups_s0_xs(2,3,4) =  -xA*cos(theta)-zA*sin(theta); 
%% Derivative with respect to chi is equal to zero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upsilon_S1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ups_s1(1,1) =  cos(theta)*cos(gamma)*sin(eta+chi)-sin(theta)*sin(gamma)*sin(chi);
Ups_s1(2,1) =  sin(gamma)*cos(chi);
Ups_s1(3,1) =  sin(theta)*cos(gamma)*sin(eta+chi)+cos(theta)*sin(gamma)*sin(chi);

Ups_s1(1,2) =   cos(theta)*cos(eta+chi);
Ups_s1(3,2) =   sin(theta)*cos(eta+chi);

Ups_s1(1,3) =   sin(theta)*sin(chi);
Ups_s1(2,3) =  -cos(chi); 
Ups_s1(3,3) =  -cos(theta)*sin(chi);

Ups_s1(1,5) =   sin(theta)*sin(chi);
Ups_s1(2,5) =  -cos(chi);
Ups_s1(3,5) =  -cos(theta)*sin(chi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient of Upsilon_S1 with respect to the state vector %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Derivative with respect to varphi is equal to zero
%% Derivative with respect to gammma is:
Ups_s1_xs(1,1,2) =  -cos(theta)*sin(gamma)*sin(eta+chi)-sin(theta)*cos(gamma)*sin(chi);
Ups_s1_xs(2,1,2) =   cos(gamma)*cos(chi);
Ups_s1_xs(3,1,2) =  -sin(theta)*sin(gamma)*sin(eta+chi)+cos(theta)*cos(gamma)*sin(chi);
%% Derivative with respect to eta is :
Ups_s1_xs(1,1,3) =  cos(theta)*cos(gamma)*cos(eta+chi);
Ups_s1_xs(3,1,3) =  sin(theta)*cos(gamma)*cos(eta+chi);

Ups_s1_xs(1,2,3) =  -cos(theta)*sin(eta+chi);
Ups_s1_xs(3,2,3) =  -sin(theta)*sin(eta+chi);
%% Derivative with respect to theta is :
Ups_s1_xs(1,1,4) = -sin(theta)*cos(gamma)*sin(eta+chi)-cos(theta)*sin(gamma)*sin(chi);
Ups_s1_xs(3,1,4) =  cos(theta)*cos(gamma)*sin(eta+chi)-sin(theta)*sin(gamma)*sin(chi);

Ups_s1_xs(1,2,4) =  -sin(theta)*cos(eta+chi);
Ups_s1_xs(3,2,4) =   cos(theta)*cos(eta+chi);

Ups_s1_xs(1,3,4) =   cos(theta)*sin(chi); 
Ups_s1_xs(3,3,4) =   sin(theta)*sin(chi);

Ups_s1_xs(1,5,4) =   cos(theta)*sin(chi);
Ups_s1_xs(3,5,4) =   sin(theta)*sin(chi);
%% Derivative with respect to chi is :
Ups_s1_xs(1,1,5) =  cos(theta)*cos(gamma)*cos(eta+chi)-sin(theta)*sin(gamma)*cos(chi);
Ups_s1_xs(2,1,5) = -sin(gamma)*sin(chi);
Ups_s1_xs(3,1,5) =  sin(theta)*cos(gamma)*cos(eta+chi)+cos(theta)*sin(gamma)*cos(chi);

Ups_s1_xs(1,2,5) =  -cos(theta)*sin(eta+chi);
Ups_s1_xs(3,2,5) =  -sin(theta)*sin(eta+chi);

Ups_s1_xs(1,3,5) =   sin(theta)*cos(chi);
Ups_s1_xs(2,3,5) =   sin(chi); 
Ups_s1_xs(3,3,5) =  -cos(theta)*cos(chi);

Ups_s1_xs(1,5,5) =   sin(theta)*cos(chi);
Ups_s1_xs(2,5,5) =   sin(chi);
Ups_s1_xs(3,5,5) =  -cos(theta)*cos(chi);

%% Add Both Contribution
Ups_s     = Ups_s0     -l1*Ups_s1; 
Ups_s_xs  = Ups_s0_xs  -l1*Ups_s1_xs;

end