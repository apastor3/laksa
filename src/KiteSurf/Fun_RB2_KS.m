function RB2 = Fun_RB2_KS(xk)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Author    : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : SB-S2 Rotaation matrix                                         %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%                                                                    %%
% Inputs:  xk    -> kite state vector                                %%
%                                                                    %%
% Outputs: RB2   -> SB-S2 rotation matrix                            %%
%                   Example VB = RB2*V2, with VB and V2 the          %%
%                   components of a vector in the SB and S2 frames   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
theta    = xk(4,1);
%
RB2      = zeros(3,3);
%
RB2(1,1) = cos(theta);   RB2(1,2) = 0;         RB2(1,3) = -sin(theta);
RB2(2,1) = 0;            RB2(2,2) = 1;         RB2(2,3) = 0;
RB2(3,1) = sin(theta);   RB2(3,2) = 0;         RB2(3,3) =  cos(theta);

end