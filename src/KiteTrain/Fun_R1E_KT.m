function [R1E R1E_xs R1E_xs2]= Fun_R1E_KT(xs,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga  and Jose A. Serrano-Iglesia           %                                      
% Language  : Matlab                                                         %
% Synopsis  : S1-SE rotation matrix                                          %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  xs    -> state vector                                   %%
%          PND   -> Dimensionless parameters                       %%
% Outputs: R1E   -> S1-SE rotation matrix                          %%
%                   Example V1 = R1E*VE, with V1 and VE the        %%
%                   components of a vector in the S1 and SE frames %%
%          R1E_xs     -> Derivative of R1E with respect to xs      %%
%          R1E_xs2    -> Derivative of R1E_xs with respect to xs   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NK         = PND.Kite.N; % Number of Kites
% Initialize the Matrix
R1E        = zeros(3,3,NK);
R1E_xs     = zeros(3,3,4,NK);
R1E_xs2    = zeros(3,3,4,4,NK);
for i=1:1:NK
    % Recover the state variable
    varphi = xs(4*(i-1)+1,1);
    gamma  = xs(4*(i-1)+2,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute the rotation matrix %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R1E(1,1,i) = cos(gamma)*cos(varphi); R1E(1,2,i) = cos(gamma)*sin(varphi);  R1E(1,3,i) = -sin(gamma);
    R1E(2,1,i) = -sin(varphi);           R1E(2,2,i) = cos(varphi);             R1E(2,3,i) = 0;
    R1E(3,1,i) = sin(gamma)*cos(varphi); R1E(3,2,i) = sin(gamma)*sin(varphi);  R1E(3,3,i) = cos(gamma);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute the derivative of R1E with respect to xs   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Derivative with respect to varphi
    R1E_xs(1,1,1,i) = -cos(gamma)*sin(varphi); R1E_xs(1,2,1,i) =  cos(gamma)*cos(varphi); 
    R1E_xs(2,1,1,i) = -cos(varphi);            R1E_xs(2,2,1,i) = -sin(varphi);            
    R1E_xs(3,1,1,i) = -sin(gamma)*sin(varphi); R1E_xs(3,2,1,i) =  sin(gamma)*cos(varphi); 
    % Compute the rotation matrix
    R1E_xs(1,1,2,i) = -sin(gamma)*cos(varphi); R1E_xs(1,2,2,i) = -sin(gamma)*sin(varphi);  R1E_xs(1,3,2,i) = -cos(gamma);
    R1E_xs(3,1,2,i) =  cos(gamma)*cos(varphi); R1E_xs(3,2,2,i) =  cos(gamma)*sin(varphi);  R1E_xs(3,3,2,i) = -sin(gamma);
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute the derivative of R1E_xs with respect to xs   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Derivative with respect to varphi
    R1E_xs2(1,1,1,1,i) = -cos(gamma)*cos(varphi); R1E_xs2(1,2,1,1,i) = -cos(gamma)*sin(varphi); 
    R1E_xs2(2,1,1,1,i) =  sin(varphi);            R1E_xs2(2,2,1,1,i) = -cos(varphi);            
    R1E_xs2(3,1,1,1,i) = -sin(gamma)*cos(varphi); R1E_xs2(3,2,1,1,i) = -sin(gamma)*sin(varphi); 
  
    R1E_xs2(1,1,2,1,i) =  sin(gamma)*sin(varphi); R1E_xs2(1,2,2,1,i) = -sin(gamma)*cos(varphi);  
    R1E_xs2(3,1,2,1,i) = -cos(gamma)*sin(varphi); R1E_xs2(3,2,2,1,i) =  cos(gamma)*cos(varphi);  
  
    % Derivative with respect to gamma
    R1E_xs2(1,1,1,2,i) =  sin(gamma)*sin(varphi); R1E_xs2(1,2,1,2,i) = -sin(gamma)*cos(varphi);             
    R1E_xs2(3,1,1,2,i) = -cos(gamma)*sin(varphi); R1E_xs2(3,2,1,2,i) =  cos(gamma)*cos(varphi); 
   
    R1E_xs2(1,1,2,2,i) = -cos(gamma)*cos(varphi); R1E_xs2(1,2,2,2,i) = -cos(gamma)*sin(varphi);  R1E_xs2(1,3,2,2,i) = sin(gamma);
    R1E_xs2(3,1,2,2,i) = -sin(gamma)*cos(varphi); R1E_xs2(3,2,2,2,i) = -sin(gamma)*sin(varphi);  R1E_xs2(3,3,2,2,i) = -cos(gamma);

end