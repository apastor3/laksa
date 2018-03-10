function R1E = Fun_R1E_KS(xk)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : S1-SE Rotation matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %%
% Inputs:  xk    -> kite state vector                                %%
%                                                                    %%
% Outputs: R1E   -> S1-SE rotation matrix                            %%
%                   Example V1 = R1E*VE, with V1 and VE the          %%
%                   components of a vector in the S1 and SE frames   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varphi = xk(1,1);
gamma  = xk(2,1);

R1E      = zeros(3,3);

R1E(1,1) = cos(gamma)*cos(varphi); R1E(1,2) = cos(gamma)*sin(varphi);  R1E(1,3) = -sin(gamma);
R1E(2,1) = -sin(varphi);           R1E(2,2) = cos(varphi);             R1E(2,3) = 0;
R1E(3,1) = sin(gamma)*cos(varphi); R1E(3,2) = sin(gamma)*sin(varphi);  R1E(3,3) = cos(gamma);

end