function R21 = Fun_R21_KA(xs)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : S2-S1 rotation matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  xs    -> state vector                                   %%
%                                                                  %%
% Outputs: R21   -> S2-S1 rotation matrix                          %%
%                   Example V2 = R21*V1, with V2 and V1 the        %%
%                   components of a vector in the S2 and S1 frames %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recover the state variable
eta    = xs(3,1);
% Compute the rotation matrix
R21(1,1) = 1;   R21(1,2) =  0;           R21(1,3) = 0;
R21(2,1) = 0;   R21(2,2) =  cos(eta);    R21(2,3) = sin(eta);
R21(3,1) = 0;   R21(3,2) = -sin(eta);    R21(3,3) = cos(eta);


end