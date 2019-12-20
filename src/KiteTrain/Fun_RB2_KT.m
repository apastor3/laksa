function [RB2 RB2_xs RB2_xs2] = Fun_RB2_KT(xs,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : SB-S2 rotation matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %%
% Inputs:  xs      -> state vector                                   %%
%          PND     -> Dimensionless parameters                       %%
% Outputs: RB2     -> SB-S2 rotation matrix                          %%
%                     Example VB = RB2*V2, with VB and V2 the        %%
%                     components of a vector in the SB and S2 frames %%
%          RB2_xs  -> Derivative of RB2 with respect to xs           %%
%          RB2_xs2 -> Derivative of RB2_xs with respect to xs        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NK         = PND.Kite.N; % Number of Kites
% Initialize the Matrix
RB2        = zeros(3,3,NK);
RB2_xs     = zeros(3,3,4,NK);
RB2_xs2    = zeros(3,3,4,4,NK);
for i=1:1:NK 
    % Recover the state variable
    theta  = xs(4*(i-1)+4,1);
    % Compute the rotation matrix
    RB2(1,1,i) = cos(theta);                         RB2(1,3,i) = -sin(theta);
                                 RB2(2,2,i) = 1;         
    RB2(3,1,i) = sin(theta);                         RB2(3,3,i) =  cos(theta);
    % Compute the derivative with respect to xs (only derivatives with respect to theta are non-zero)
    RB2_xs(1,1,4,i) = -sin(theta);                  RB2_xs(1,3,4,i) = -cos(theta);         
    RB2_xs(3,1,4,i) =  cos(theta);                  RB2_xs(3,3,4,i) = -sin(theta);    

    % Compute the derivative of RB"_xs with respect to xs (only derivatives with respect to theta are non-zero)
    RB2_xs2(1,1,4,4,i) = -cos(theta);               RB2_xs2(1,3,4,4,i) =  sin(theta);         
    RB2_xs2(3,1,4,4,i) = -sin(theta);               RB2_xs2(3,3,4,4,i) = -cos(theta);       
end

end