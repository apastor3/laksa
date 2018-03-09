function RB2 = Fun_RB2_KA(xs)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : SB-S2 rotation matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  xs    -> state vector                                   %%
%                                                                  %%
% Outputs: RB2   -> SB-S2 rotation matrix                          %%
%                   Example VB = RB2*V2, with VB and V2 the        %%
%                   components of a vector in the SB and S2 frames %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recover the state variable
theta  = xs(4,1);
% Compute the rotation matrix
RB2(1,1) = cos(theta);   RB2(1,2) = 0;         RB2(1,3) = -sin(theta);
RB2(2,1) = 0;            RB2(2,2) = 1;         RB2(2,3) = 0;
RB2(3,1) = sin(theta);   RB2(3,2) = 0;         RB2(3,3) =  cos(theta);

end