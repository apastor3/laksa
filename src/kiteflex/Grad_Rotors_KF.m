function  [SG_xs SG_xc CG_xs CG_xc ] = Grad_Rotors_KF(xs,xc,R_KE,Grad_R_KE,SK_xs,SK_xc,CK_xs,CK_xc,PND)
  
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Gradient of Rotor-Matrices                                     %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:                                                                   %%
%         State vectors xs  -> [xs.gamma xs.varphi xs.theta xs.psi xs.phi]  %% 
%         Control vector xc -> [xc.lt xc.lb dexc.lta xc.eta]                %%
%         Rotation matrix R_KE, and its gradient Grad_R_KE                  %%
%         Gradients of Kite matrices SK_xs,SK_xc,CK_xs,CK_xc                %%
%         PND               -> Dimensionless parameters                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:                                                                   %%
% Gradients                                                                 %%
%       SG_xs  = partial SG/partial x_s                                     %%
%       SG_xc  = partial SG/partial x_c                                     %%
%       CG_xs  = partial CG/partial x_s                                     %%
%       CG_xc  = partial CG/partial x_c                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Number of bars and generators

NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators

% State and contrl vector dimensions
Nvar   = 2*NR+3;
Nvar_p = 2*NR+3+NG;
Nc0    = 4;

SG_xs  = zeros(3,Nvar_p,NG,Nvar);  % Initialize the matrix
SG_xc  = zeros(3,Nvar_p,NG,Nc0);    % Initialize the matrix
CG_xs  = zeros(3,Nc0,NG,Nvar);    % Initialize the matrix
CG_xc  = zeros(3,Nc0,NG,Nc0);      % Initialize the matrix

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matrix SG
%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:NG
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Matrix SG0 = R_EK*AUX  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SG0_xs   = zeros(3,Nvar,Nvar);
    AUX      = zeros(3,3);
    Grad_AUX = zeros(3,3,3);
      
    AUX(1,1)  = PND.Gen.z(i)*cos(xs.phi)+PND.Gen.y(i)*sin(xs.phi);
    AUX(1,2)  = PND.Gen.z(i)*cos(xs.theta)*sin(xs.phi)-PND.Gen.y(i)*cos(xs.theta)*cos(xs.phi);
    
    AUX(2,1)  = -PND.Gen.x(i)*sin(xs.phi);
    AUX(2,2)  =  PND.Gen.x(i)*cos(xs.theta)*cos(xs.phi)+PND.Gen.z(i)*sin(xs.theta);
    AUX(2,3)  = -PND.Gen.z(i);
  
    AUX(3,1)  = -PND.Gen.x(i)*cos(xs.phi);
    AUX(3,2)  = -PND.Gen.y(i)*sin(xs.theta)-PND.Gen.x(i)*cos(xs.theta)*sin(xs.phi);
    AUX(3,3)  =  PND.Gen.y(i);
    
    
    % Derivative of AUX with respect to Theta
    Grad_AUX(1,2,1)  = -PND.Gen.z(i)*sin(xs.theta)*sin(xs.phi)+PND.Gen.y(i)*sin(xs.theta)*cos(xs.phi);
    Grad_AUX(2,2,1)  = -PND.Gen.x(i)*sin(xs.theta)*cos(xs.phi)+PND.Gen.z(i)*cos(xs.theta);
    Grad_AUX(3,2,1)  = -PND.Gen.y(i)*cos(xs.theta)+PND.Gen.x(i)*sin(xs.theta)*sin(xs.phi);
  
    % Derivative of AUX with respect to Phi
    Grad_AUX(1,1,3)  = -PND.Gen.z(i)*sin(xs.phi)+PND.Gen.y(i)*cos(xs.phi);
    Grad_AUX(1,2,3)  = PND.Gen.z(i)*cos(xs.theta)*cos(xs.phi)+PND.Gen.y(i)*cos(xs.theta)*sin(xs.phi);
    
    Grad_AUX(2,1,3)  = -PND.Gen.x(i)*cos(xs.phi);
    Grad_AUX(2,2,3)  = -PND.Gen.x(i)*cos(xs.theta)*sin(xs.phi);
  
    Grad_AUX(3,1,3)  =  PND.Gen.x(i)*sin(xs.phi);
    Grad_AUX(3,2,3)  = -PND.Gen.x(i)*cos(xs.theta)*cos(xs.phi);
  
    
    % Compute Matrices
    for j=1:1:3 % Derivatives with respect to theta, psi, and phi   
        SG_xs(:,2*NR+1:2*NR+3,i,2*NR+j)  = Grad_R_KE(:,:,j)'*AUX   + R_KE'* Grad_AUX(:,:,j);   
    end
 
    SG_xs(:,:,i,:) = SK_xs(:,:,:) + squeeze(SG_xs(:,:,i,:));  
    SG_xc(:,:,i,:) = SK_xc(:,:,:);    
    
    CG_xs(:,:,i,:) = CK_xs(:,:,:);
    CG_xc(:,:,i,:) = CK_xc(:,:,:);
end






end

