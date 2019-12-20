function [R21 R21_xs R21_xs2] = Fun_R21_KT(xs,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : S2-S1 rotation matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %%
% Inputs:  xs     -> state vector                                   %%
%          PND    -> Dimensionless parameters                       %%
% Outputs: R21    -> S2-S1 rotation matrix                          %%
%                   Example V2 = R21*V1, with V2 and V1 the         %%
%                   components of a vector in the S2 and S1 frames  %%
%          R21_xs -> Derivative of R21 with respect to xs           %%
%          R21_xs2 -> Derivative of R21_xs with respect to xs       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NK     = PND.Kite.N; % Number of Kites

% Initialize the Matrix
R21     = zeros(3,3,NK);
R21_xs  = zeros(3,3,4,NK);
R21_xs2 = zeros(3,3,4,4,NK);
for i=1:1:NK
    % Recover the state variable
    eta    = xs(4*(i-1)+3,1);
    % Compute the rotation matrix
    R21(1,1,i) = 1;  
                      R21(2,2,i) =  cos(eta);    R21(2,3,i) = sin(eta);
                      R21(3,2,i) = -sin(eta);    R21(3,3,i) = cos(eta);
    % Compute the derivative of R21 with respect to xs (only derivatives with respect to eta are non-zero)  
    R21_xs(2,2,3,i) = -sin(eta);    R21_xs(2,3,3,i) =  cos(eta);
    R21_xs(3,2,3,i) = -cos(eta);    R21_xs(3,3,3,i) = -sin(eta);
 
    % Compute the derivative of  R21_xs with respect to xs (only derivatives with respect to eta are non-zero)  
    R21_xs2(2,2,3,3,i) = -cos(eta);    R21_xs2(2,3,3,3,i) = -sin(eta);
    R21_xs2(3,2,3,3,i) =  sin(eta);    R21_xs2(3,3,3,3,i) = -cos(eta);
 
end

end