function [SR_xs SR_xc CR_xs OmR_xs] = Grad_Rod_KF(xs,xc,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Gradients of Rod Matrices                                      %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                                    %%
%     State vectors   xs -> [xs.gamma xs.varxs.phi thxc.eta xs.psi xs.phi ] %%
%     Control vector  xc -> [xc.lt xc.lb dexc.lta xc.eta]                   %%
%     Dimensionless Parameters -> PND                                       %%
% Output:                                                                   %%
% Gradients                                                                 %%
%       SR_xs  = partial Si/partial x_s                                     %%
%       SR_xc  = partial Si/partial x_c                                     %%
%       CR_xs  = partial Ci/partial x_s                                     %%
%       OmR_xs = partial Ci/partial x_c                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Number of bar
NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators

Nvar   = 2*NR+3;
Nvar_p = 2*NR+3+NG;
Nc0    = 4;

SR_xs  = zeros(3,Nvar_p,NR,Nvar); % Initialize the matrix
SR_xc  = zeros(3,Nvar_p,NR,Nc0);    % Initialize the matrix

CR_xs  = zeros(3,Nc0,NR,Nvar); % Initialize the matrix

OmR_xs = zeros(3,Nvar_p,NR,Nvar);

for i=1:1:NR
    
    % Gradient of Omi
    OmR_xs(1,NR+i,i,i) =  cos(xs.gamma(i));
    OmR_xs(3,NR+i,i,i) = -sin(xs.gamma(i));
    
    for j=1:1:i
        if j==i
            g = 0.5;
        else
            g = 1;
        end
        
        %%% Gradients of Si
        
        % Derivative with respect to gamma
        SR_xs(1,j,i,j)     =  xc.lt*g*cos(xs.gamma(j))*cos(xs.varphi(j));
        SR_xs(1,NR+j,i,j)  = -xc.lt*g*sin(xs.gamma(j))*sin(xs.varphi(j));
        SR_xs(2,j,i,j)     =  xc.lt*g*cos(xs.gamma(j))*sin(xs.varphi(j));
        SR_xs(2,NR+j,i,j)  =  xc.lt*g*sin(xs.gamma(j))*cos(xs.varphi(j));
        SR_xs(3,j,i,j)     =  xc.lt*g*sin(xs.gamma(j));
        
        % Derivative with respect to varphi
        SR_xs(1,j,i,NR+j)    = -xc.lt*g*sin(xs.gamma(j))*sin(xs.varphi(j));
        SR_xs(1,NR+j,i,NR+j) =  xc.lt*g*cos(xs.gamma(j))*cos(xs.varphi(j));
        SR_xs(2,j,i,NR+j)    =  xc.lt*g*sin(xs.gamma(j))*cos(xs.varphi(j));
        SR_xs(2,NR+j,i,NR+j) =  xc.lt*g*cos(xs.gamma(j))*sin(xs.varphi(j)); 
        
        % Derivative with respect to lt
        SR_xc(1,j,i,1)    =  g*sin(xs.gamma(j))*cos(xs.varphi(j));
        SR_xc(1,NR+j,i,1) =  g*cos(xs.gamma(j))*sin(xs.varphi(j));
        SR_xc(2,j,i,1)    =  g*sin(xs.gamma(j))*sin(xs.varphi(j));
        SR_xc(2,NR+j,i,1) = -g*cos(xs.gamma(j))*cos(xs.varphi(j));
        SR_xc(3,j,i,1)    = -g*cos(xs.gamma(j));
        
        %%% Gradients of Sc
        
        % Derivative with respect to gamma
        CR_xs(1,1,i,j)   =   + g*sin(xs.gamma(j))*cos(xs.varphi(j));  
        CR_xs(2,1,i,j)   =   + g*sin(xs.gamma(j))*sin(xs.varphi(j));  
        CR_xs(3,1,i,j)   =   - g*cos(xs.gamma(j));
       
        % Derivative with respect to varphi 
        CR_xs(1,1,i,NR+j) =  + g*cos(xs.gamma(j))*sin(xs.varphi(j));  
        CR_xs(2,1,i,NR+j) =  - g*cos(xs.gamma(j))*cos(xs.varphi(j));  
    end
end


end

