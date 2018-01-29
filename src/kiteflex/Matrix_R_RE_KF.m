function [R_RE ] = Matrix_R_RE_KF(xs,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : SR-SE rotation matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Inputs    xs -> state vector                                           %% 
%%   Output    R_RE       -> Rotation matrix SR-SE                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recover some important dimensions
NR          = PND.Num.N;                  % Number of Rods
NG          = PND.Gen.Num;                % Number of Generators
Nvar        = 2*PND.Num.N+3;              % Number of variables 
% Initialize matrices
R_RE        = zeros(3,3,NR);
% Compute rotation matrix 
for i=1:1:NR
    gamma  = xs.gamma(i);
    varphi = xs.varphi(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%      R_RE Rotation Matrix                    %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First Row
    R_RE(1,1,i) = cos(gamma)*cos(varphi);
    R_RE(1,2,i) = cos(gamma)*sin(varphi);
    R_RE(1,3,i) = sin(gamma);
    % Second Row
    R_RE(2,1,i) = -sin(varphi);
    R_RE(2,2,i) =  cos(varphi);
    R_RE(2,3,i) =  0;
    % Third Row
    R_RE(3,1,i) = -sin(gamma)*cos(varphi);
    R_RE(3,2,i) = -sin(gamma)*sin(varphi);
    R_RE(3,3,i) = cos(gamma);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

