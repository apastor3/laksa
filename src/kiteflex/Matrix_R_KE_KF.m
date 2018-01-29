function [R_KE Grad_R_KE  ] = Matrix_R_KE_KF(xs)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : SK-SE Rotation Matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Inputs    xs -> state vector                                           %% 
%%   Output    R_BE       -> Rotation matrix Kite (Body) -Earth             %%
%              Grad_R_BE  -> derivative of R_BE with respect to Euler angles%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Theta = xs.theta;
Psi   = xs.psi;
Phi   = xs.phi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%      Kite - Earth Rotation                   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Row
R_KE(1,1) = cos(Psi)*cos(Theta);
R_KE(1,2) = sin(Psi)*cos(Theta);
R_KE(1,3) = -sin(Theta);
% Second Row
R_KE(2,1) = cos(Psi)*sin(Theta)*sin(Phi)-sin(Psi)*cos(Phi);
R_KE(2,2) = sin(Psi)*sin(Theta)*sin(Phi)+cos(Psi)*cos(Phi);
R_KE(2,3) = cos(Theta)*sin(Phi);
% Third Row
R_KE(3,1) = cos(Psi)*sin(Theta)*cos(Phi)+sin(Psi)*sin(Phi);
R_KE(3,2) = sin(Psi)*sin(Theta)*cos(Phi)-cos(Psi)*sin(Phi);
R_KE(3,3) = cos(Theta)*cos(Phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Grad_R_KE = zeros(3,3,3);

% Theta Derivative
Grad_R_KE(1,1,1) = -cos(Psi)*sin(Theta);
Grad_R_KE(1,2,1) = -sin(Psi)*sin(Theta);
Grad_R_KE(1,3,1) = -cos(Theta);
Grad_R_KE(2,1,1) = cos(Psi)*cos(Theta)*sin(Phi);
Grad_R_KE(2,2,1) = sin(Psi)*cos(Theta)*sin(Phi);
Grad_R_KE(2,3,1) =-sin(Theta)*sin(Phi);
Grad_R_KE(3,1,1) = cos(Psi)*cos(Theta)*cos(Phi);
Grad_R_KE(3,2,1) = sin(Psi)*cos(Theta)*cos(Phi);
Grad_R_KE(3,3,1) = -sin(Theta)*cos(Phi);
% Psi Derivative
Grad_R_KE(1,1,2) = -sin(Psi)*cos(Theta);
Grad_R_KE(1,2,2) =  cos(Psi)*cos(Theta);
Grad_R_KE(2,1,2) = -sin(Psi)*sin(Theta)*sin(Phi)-cos(Psi)*cos(Phi);
Grad_R_KE(2,2,2) =  cos(Psi)*sin(Theta)*sin(Phi)-sin(Psi)*cos(Phi);
Grad_R_KE(3,1,2) = -sin(Psi)*sin(Theta)*cos(Phi)+cos(Psi)*sin(Phi);
Grad_R_KE(3,2,2) =  cos(Psi)*sin(Theta)*cos(Phi)+sin(Psi)*sin(Phi);
% Phi Derivative
Grad_R_KE(2,1,3) = cos(Psi)*sin(Theta)*cos(Phi)+sin(Psi)*sin(Phi);
Grad_R_KE(2,2,3) = sin(Psi)*sin(Theta)*cos(Phi)-cos(Psi)*sin(Phi);
Grad_R_KE(2,3,3) = cos(Theta)*cos(Phi);
Grad_R_KE(3,1,3) = -cos(Psi)*sin(Theta)*sin(Phi)+sin(Psi)*cos(Phi);
Grad_R_KE(3,2,3) = -sin(Psi)*sin(Theta)*sin(Phi)-cos(Psi)*cos(Phi);
Grad_R_KE(3,3,3) = -cos(Theta)*sin(Phi);

