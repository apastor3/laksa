function R32 = Fun_R32_KS(xk)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Author    : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : S3-S2 Rotation matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%                                                                    %%
% Inputs:  xk    -> kite state vector                                %%
%                                                                    %%
% Outputs: R32   -> S3-S2 rotation matrix                            %%
%                   Example V3 = R32*V2, with V3 and V2 the          %%
%                   components of a vector in the S3 and S2 frames   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variable
chi      = xk(5,1);
% Initialize
R32      = zeros(3,3);

R32(1,1) = 1;   R32(1,2) =      0;        R32(1,3) = 0;
R32(2,1) = 0;   R32(2,2) =  cos(chi);     R32(2,3) = sin(chi);
R32(3,1) = 0;   R32(3,2) = -sin(chi);     R32(3,3) = cos(chi);

end