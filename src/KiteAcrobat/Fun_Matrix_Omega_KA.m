function [Phi Phi_xs] = Fun_Matrix_Omega_KA(xs)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Matrix Phi and its gradient                                    %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  xs     -> state vector                                  %%
%                                                                  %%
% Outputs: Phi    -> matrix for the computation of the components  %%
%                  in SB of the dimensionless angular velocity of  %%
%                  of the kite with respect to SE:                 %%
%                  omega = Phi*xs_p                                %%
%          Phi_xs -> tensor with the partial derivatives of Phi    %%
%                  with respect to the state vector components     %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recover the variables of the state Vector
varphi = xs(1,1);
gamma  = xs(2,1);
eta    = xs(3,1);
theta  = xs(4,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            Compute matrix phi                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi(1,1) = -cos(gamma)*cos(eta)*sin(theta)-sin(gamma)*cos(theta);   Phi(1,2) =  sin(eta)*sin(theta);  Phi(1,3) = cos(theta); Phi(1,4) = 0;
Phi(2,1) =  cos(gamma)*sin(eta);                                    Phi(2,2) =  cos(eta);             Phi(2,3) = 0;          Phi(2,4) = 1;
Phi(3,1) =  cos(gamma)*cos(eta)*cos(theta)-sin(gamma)*sin(theta);   Phi(3,2) = -sin(eta)*cos(theta);  Phi(3,3) = sin(theta); Phi(3,4) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient of Omega with respect to the state vector     %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi_xs  = zeros(3,4,4); % partial Omega/partial x_s
%% Derivative with respect to varphi is equal to zero
%% Derivative with respect to gamma is:
Phi_xs(1,1,2) =  sin(gamma)*cos(eta)*sin(theta)-cos(gamma)*cos(theta);   
Phi_xs(2,1,2) = -sin(gamma)*sin(eta);                                    
Phi_xs(3,1,2) = -sin(gamma)*cos(eta)*cos(theta)-cos(gamma)*sin(theta);   
%% Derivative with respect to eta is:
Phi_xs(1,1,3) =  cos(gamma)*sin(eta)*sin(theta);   Phi_xs(1,2,3) =  cos(eta)*sin(theta);  
Phi_xs(2,1,3) =  cos(gamma)*cos(eta);              Phi_xs(2,2,3) = -sin(eta);             
Phi_xs(3,1,3) = -cos(gamma)*sin(eta)*cos(theta);   Phi_xs(3,2,3) = -cos(eta)*cos(theta);  
%% Derivative with respect to theta is:
Phi_xs(1,1,4) = -cos(gamma)*cos(eta)*cos(theta)+sin(gamma)*sin(theta);   Phi_xs(1,2,4) =  sin(eta)*cos(theta);  Phi_xs(1,3,4) = -sin(theta);           
Phi_xs(3,1,4) = -cos(gamma)*cos(eta)*sin(theta)-sin(gamma)*cos(theta);   Phi_xs(3,2,4) =  sin(eta)*sin(theta);  Phi_xs(3,3,4) =  cos(theta); 



end