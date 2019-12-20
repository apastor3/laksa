function R_KE = Matrix_R_KE_KE(eulerK)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute Kite-Earth Rotation matrix                             %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:                                                                  %
%           Euler angles -> [Phi Theta Psi]                                  %
%   Outputs:                                                                 %
%           R_KE  -> Kite-Earth Rotation matrix                              % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%      Kite - Earth Rotation                   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover Euler angles
Phi   = eulerK(1,1);
Theta = eulerK(2,1);
Psi   = eulerK(3,1);

% Compute rotation matrix
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
