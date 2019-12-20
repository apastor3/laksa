function [ZX Grad_ZX Grad2_ZX] = Fun_ZX_Full_KT(R1E,R21,RB2,R1E_xs,R21_xs,RB2_xs,R1E_xs2,R21_xs2,RB2_xs2,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : Distances d_j and their gradients                              %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                            %
% Inputs:                                                                    %
%          R1E,R21,RB2              -> S1-SE rotation matrix                 % 
%          R1E_xs,R21_xs,RB2_xs,    -> Derivative of R with respect to xs    %
%          R1E_xs2,R21_xs2,RB2_xs2, -> Derivative of R_xs with respect to xs %
%          PND                      -> Dimensionless parameters              %
% Outputs: d                        -> Coordinates zeta, xi                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Initialize
ZX         = zeros(2,PND.Kite.N);
Grad_ZX    = zeros(2,4*PND.Kite.N,PND.Kite.N);
Grad2_ZX   = zeros(2,4*PND.Kite.N,4*PND.Kite.N,PND.Kite.N);

% Compute first kite
ZX(1,1)    = 0; 
ZX(2,1)    = sqrt(PND.Tether.L(1)^2-PND.Tether.YA(1)^2);
for i=2:1:PND.Kite.N
   % Initialize 
   C        = zeros(2,3); % (Plus/Minus, x/y/z,Kite)
   Grad_C   = zeros(2,3,4*PND.Kite.N); 
   Grad2_C  = zeros(2,3,4*PND.Kite.N,4*PND.Kite.N);
   
   % Auxiliary Calculation   
   R2E_i  = squeeze(R21(:,:,i))*squeeze(R1E(:,:,i));
   REB_iM = squeeze(R1E(:,:,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))';
  
   ATT_C(:,1)   = [ PND.Tether.XC(i-1)  PND.Tether.YC(i-1) PND.Tether.ZC(i-1)]';
   ATT_C(:,2)   = [ PND.Tether.XC(i-1) -PND.Tether.YC(i-1) PND.Tether.ZC(i-1)]';
  
   ATT_A(:,1)   = [0                    PND.Tether.YA(i)    0]'; 
   ATT_A(:,2)   = [0                   -PND.Tether.YA(i)    0]'; 
 
   
   for k=1:1:2
     C(k,:)     = ATT_A(:,k)-R2E_i*REB_iM*ATT_C(:,k);
     Kappa(k)   = sqrt( PND.Tether.L(i)^2 - C(k,1)^2);
   end
   
   [ZX(:,i) R2 Aux_S1 Aux_S2 Aux_Sq] = Circle(C(1,2),C(1,3),Kappa(1),C(2,2),C(2,3),Kappa(2));
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%  Compute Gradients of C                  %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for k=1:1:2
       % Derivative with respect to varphi^(n-1)
       Ind           = 4*(i-2)+1;
       Grad_C(k,:,Ind) = -R2E_i*squeeze(R1E_xs(:,:,1,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);

       % Derivative with respect to gamma^(n-1)
       Ind           = 4*(i-2)+2;
       Grad_C(k,:,Ind) = -R2E_i*squeeze(R1E_xs(:,:,2,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);

       % Derivative with respect to eta^(n-1)
       Ind           = 4*(i-2)+3;
       Grad_C(k,:,Ind) = -R2E_i*squeeze(R1E(:,:,i-1))'*squeeze(R21_xs(:,:,3,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);

       % Derivative with respect to theta^(n-1)
       Ind           = 4*(i-2)+4;
       Grad_C(k,:,Ind) = -R2E_i*squeeze(R1E(:,:,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2_xs(:,:,4,i-1))'*ATT_C(:,k);

       % Derivative with respect to varphi^n
       Ind           = 4*(i-1)+1;
       Grad_C(k,:,Ind) = -R21(:,:,i)*R1E_xs(:,:,1,i)*REB_iM*ATT_C(:,k);

       % Derivative with respect to gamma^n
       Ind           = 4*(i-1)+2;
       Grad_C(k,:,Ind) = -R21(:,:,i)*R1E_xs(:,:,2,i)*REB_iM*ATT_C(:,k);

       % Derivative with respect to eta^n
       Ind           = 4*(i-1)+3;
       Grad_C(k,:,Ind) = -R21_xs(:,:,3,i)*R1E(:,:,i)*REB_iM*ATT_C(:,k);
   end
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%           Gradient of ZX           %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   for k=1:1:2
      Grad_Kappa(k,:) = -C(k,1)*Grad_C(k,1,:)/Kappa(k); 
   end
   
    Grad_R2(1,:)   = 2*(C(1,2)-C(2,2))*(Grad_C(1,2,:)-Grad_C(2,2,:))...
                       + 2*(C(1,3)-C(2,3))*(Grad_C(1,3,:)-Grad_C(2,3,:));
    Grad_S1        =  2*(Kappa(1)+Kappa(2))*(Grad_Kappa(1,:) + Grad_Kappa(2,:)) - Grad_R2(1,:);
    Grad_S2        = -2*(Kappa(1)-Kappa(2))*(Grad_Kappa(1,:) - Grad_Kappa(2,:)) + Grad_R2(1,:);
    Grad_Sq        = (Aux_S2*Grad_S1 + Aux_S1*Grad_S2)/(2*Aux_Sq);

    Grad_ZX(1,:,i) =                  (Grad_C(1,2,:)+Grad_C(2,2,:))/2;

    Grad_ZX(1,:,i) = Grad_ZX(1,:,i) + (Kappa(1)^2-Kappa(2)^2)/(2*R2)*squeeze((Grad_C(2,2,:) - Grad_C(1,2,:)))';
    Grad_ZX(1,:,i) = Grad_ZX(1,:,i) + (Kappa(1)*Grad_Kappa(1,:) - Kappa(2)*Grad_Kappa(2,:))*(C(2,2)-C(1,2))/R2;
    Grad_ZX(1,:,i) = Grad_ZX(1,:,i) - (Kappa(1)^2-Kappa(2)^2)/(2*R2^2)*(C(2,2)-C(1,2))*Grad_R2(1,:);

    Grad_ZX(1,:,i) = Grad_ZX(1,:,i) +  Grad_Sq(1,:)/(2*R2)*(C(2,3)-C(1,3));
    Grad_ZX(1,:,i) = Grad_ZX(1,:,i) +  Aux_Sq/(2*R2)*squeeze((Grad_C(2,3,:)-Grad_C(1,3,:)))';
    Grad_ZX(1,:,i) = Grad_ZX(1,:,i) -  Aux_Sq*(C(2,3) - C(1,3))/(2*R2^2)*Grad_R2(1,:);

    Grad_ZX(2,:,i) =                  (Grad_C(1,3,:)+Grad_C(2,3,:))/2;
    Grad_ZX(2,:,i) = Grad_ZX(2,:,i) + (Kappa(1)^2-Kappa(2)^2)/(2*R2)*squeeze((Grad_C(2,3,:) - Grad_C(1,3,:)))';
    Grad_ZX(2,:,i) = Grad_ZX(2,:,i) + (Kappa(1)*Grad_Kappa(1,:) - Kappa(2)*Grad_Kappa(2,:))*(C(2,3)-C(1,3))/R2;
    Grad_ZX(2,:,i) = Grad_ZX(2,:,i) - (Kappa(1)^2-Kappa(2)^2)/(2*R2^2)*(C(2,3)-C(1,3))*Grad_R2(1,:);
    Grad_ZX(2,:,i) = Grad_ZX(2,:,i) +  Grad_Sq(1,:)/(2*R2)*(C(1,2)-C(2,2));
    Grad_ZX(2,:,i) = Grad_ZX(2,:,i) +  Aux_Sq/(2*R2)*squeeze((Grad_C(1,2,:)-Grad_C(2,2,:)))';
    Grad_ZX(2,:,i) = Grad_ZX(2,:,i) -  Aux_Sq*(C(1,2) - C(2,2))/(2*R2^2)*Grad_R2(1,:);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%  Compute the Derivatives of the Gradients of the Distances  with respect to xs     %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   R2E_i_varphi  = R21(:,:,i)*R1E_xs(:,:,1,i);
   R2E_i_gamma   = R21(:,:,i)*R1E_xs(:,:,2,i);
   R2E_i_eta     = R21_xs(:,:,3,i)*R1E(:,:,i);
   
   REB_iM_varphi = squeeze(R1E_xs(:,:,1,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))';
   REB_iM_gamma  = squeeze(R1E_xs(:,:,2,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))';
   REB_iM_eta    = squeeze(R1E(:,:,i-1))'*squeeze(R21_xs(:,:,3,i-1))'*squeeze(RB2(:,:,i-1))';
   REB_iM_theta  = squeeze(R1E(:,:,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2_xs(:,:,4,i-1))';
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Derivative with respect to varphi^(n-1) %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for k=1:1:2
       Aux0                = squeeze(R1E_xs(:,:,1,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind                 =  4*(i-2)+1;

       Ind2                =  4*(i-2)+1;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs2(:,:,1,1,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+2;  %% Derivative with respect to gamma^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs2(:,:,1,2,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+3;  %% Derivative with respect to eta^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs(:,:,1,i-1))'*squeeze(R21_xs(:,:,3,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+4;  %% Derivative with respect to eta^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs(:,:,1,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2_xs(:,:,4,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-1)+1;  %% Derivative with respect to varphi^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_varphi*Aux0;
       Ind2                =  4*(i-1)+2;  %% Derivative with respect to gamma^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_gamma*Aux0;
       Ind2                =  4*(i-1)+3;  %% Derivative with respect to eta^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_eta*Aux0;

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       % Derivative with respect to gamma^(n-1)   %%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Aux0                = squeeze(R1E_xs(:,:,2,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind                 =  4*(i-2)+2;

       Ind2                =  4*(i-2)+1;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs2(:,:,2,1,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+2;  %% Derivative with respect to gamma^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs2(:,:,2,2,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+3;  %% Derivative with respect to eta^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs(:,:,2,i-1))'*squeeze(R21_xs(:,:,3,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+4;  %% Derivative with respect to theta^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs(:,:,2,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2_xs(:,:,4,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-1)+1;  %% Derivative with respect to varphi^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_varphi*Aux0;
       Ind2                =  4*(i-1)+2;  %% Derivative with respect to gamma^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_gamma*Aux0;
       Ind2                =  4*(i-1)+3;  %% Derivative with respect to eta^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_eta*Aux0;
   
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Derivative with respect to eta^(n-1)     %%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Aux0                = squeeze(R1E(:,:,i-1))'*squeeze(R21_xs(:,:,3,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind                 = 4*(i-2)+3;

       Ind2                =  4*(i-2)+1;  %% Derivative with respect to varphi^(n-1) 
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs(:,:,1,i-1))'*squeeze(R21_xs(:,:,3,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+2;  %% Derivative with respect to gamma^(n-1) 
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs(:,:,2,i-1))'*squeeze(R21_xs(:,:,3,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+3;  %% Derivative with respect to eta^(n-1) 
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E(:,:,i-1))'*squeeze(R21_xs2(:,:,3,3,i-1))'*squeeze(RB2(:,:,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+4;  %% Derivative with respect to theta^(n-1) 
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E(:,:,i-1))'*squeeze(R21_xs(:,:,3,i-1))'*squeeze(RB2_xs(:,:,4,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-1)+1;  %% Derivative with respect to varphi^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_varphi*Aux0;
       Ind2                =  4*(i-1)+2;  %% Derivative with respect to gamma^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_gamma*Aux0;
       Ind2                =  4*(i-1)+3;  %% Derivative with respect to eta^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_eta*Aux0;
   
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Derivative with respect to theta^(n-1)     %%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Aux0                = squeeze(R1E(:,:,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2_xs(:,:,4,i-1))'*ATT_C(:,k);
       Ind                 = 4*(i-2)+4;

       Ind2                =  4*(i-2)+1;  %% Derivative with respect to varphi^(n-1)   
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs(:,:,1,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2_xs(:,:,4,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+2;  %% Derivative with respect to gamma^(n-1)   
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E_xs(:,:,2,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2_xs(:,:,4,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+3;  %% Derivative with respect to eta^(n-1)   
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E(:,:,i-1))'*squeeze(R21_xs(:,:,3,i-1))'*squeeze(RB2_xs(:,:,4,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-2)+4;  %% Derivative with respect to theta^(n-1)   
       Grad2_C(k,:,Ind,Ind2) = -R2E_i*squeeze(R1E(:,:,i-1))'*squeeze(R21(:,:,i-1))'*squeeze(RB2_xs2(:,:,4,4,i-1))'*ATT_C(:,k);
       Ind2                =  4*(i-1)+1;  %% Derivative with respect to varphi^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_varphi*Aux0;
       Ind2                =  4*(i-1)+2;  %% Derivative with respect to gamma^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_gamma*Aux0;
       Ind2                =  4*(i-1)+3;  %% Derivative with respect to eta^(n)
       Grad2_C(k,:,Ind,Ind2) = -R2E_i_eta*Aux0;
 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Derivative with respect to varphi^n        %%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Aux0                = R21(:,:,i)*R1E_xs(:,:,1,i);
       Ind                 = 4*(i-1)+1;

       Ind2                =  4*(i-2)+1;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_varphi*ATT_C(:,k);
       Ind2                =  4*(i-2)+2;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_gamma*ATT_C(:,k);
       Ind2                =  4*(i-2)+3;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_eta*ATT_C(:,k);
       Ind2                =  4*(i-2)+4;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_theta*ATT_C(:,k);
       Ind2                =  4*(i-1)+1;  %% Derivative with respect to varphi^n  
       Grad2_C(k,:,Ind,Ind2) = -R21(:,:,i)*R1E_xs2(:,:,1,1,i)*REB_iM*ATT_C(:,k);
       Ind2                =  4*(i-1)+2;  %% Derivative with respect to gamma^n  
       Grad2_C(k,:,Ind,Ind2) = -R21(:,:,i)*R1E_xs2(:,:,1,2,i)*REB_iM*ATT_C(:,k);
       Ind2                =  4*(i-1)+3;  %% Derivative with respect to eta^n  
       Grad2_C(k,:,Ind,Ind2) = -R21_xs(:,:,3,i)*R1E_xs(:,:,1,i)*REB_iM*ATT_C(:,k);

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Derivative with respect to gamma^n           %%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Aux0                = R21(:,:,i)*R1E_xs(:,:,2,i);  
       Ind                 = 4*(i-1)+2;

       Ind2                =  4*(i-2)+1;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_varphi*ATT_C(:,k);
       Ind2                =  4*(i-2)+2;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_gamma*ATT_C(:,k);
       Ind2                =  4*(i-2)+3;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_eta*ATT_C(:,k);
       Ind2                =  4*(i-2)+4;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_theta*ATT_C(:,k);
       Ind2                =  4*(i-1)+1;  %% Derivative with respect to varphi^n    
       Grad2_C(k,:,Ind,Ind2) = -R21(:,:,i)*R1E_xs2(:,:,2,1,i)*REB_iM*ATT_C(:,k);
       Ind2                =  4*(i-1)+2;  %% Derivative with respect to gamma^n    
       Grad2_C(k,:,Ind,Ind2) = -R21(:,:,i)*R1E_xs2(:,:,2,2,i)*REB_iM*ATT_C(:,k);
       Ind2                =  4*(i-1)+3;  %% Derivative with respect to eta^n    
       Grad2_C(k,:,Ind,Ind2) = -R21_xs(:,:,3,i)*R1E_xs(:,:,2,i)*REB_iM*ATT_C(:,k);

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Derivative with respect to eta^n                %%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Aux0                = R21_xs(:,:,3,i)*R1E(:,:,i);  
       Ind                 = 4*(i-1)+3;

       Ind2                =  4*(i-2)+1;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_varphi*ATT_C(:,k);
       Ind2                =  4*(i-2)+2;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_gamma*ATT_C(:,k);
       Ind2                =  4*(i-2)+3;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_eta*ATT_C(:,k);
       Ind2                =  4*(i-2)+4;  %% Derivative with respect to varphi^(n-1)
       Grad2_C(k,:,Ind,Ind2) = -Aux0*REB_iM_theta*ATT_C(:,k);

       Ind2                =  4*(i-1)+1;  %% Derivative with respect to varphi^n    
       Grad2_C(k,:,Ind,Ind2) = -R21_xs(:,:,3,i)*R1E_xs(:,:,1,i)*REB_iM*ATT_C(:,k);
       Ind2                =  4*(i-1)+2;  %% Derivative with respect to gamma^n    
       Grad2_C(k,:,Ind,Ind2) = -R21_xs(:,:,3,i)*R1E_xs(:,:,2,i)*REB_iM*ATT_C(:,k);
       Ind2                =  4*(i-1)+3;  %% Derivative with respect to eta^n    
       Grad2_C(k,:,Ind,Ind2) = -R21_xs2(:,:,3,3,i)*R1E(:,:,i)*REB_iM*ATT_C(:,k);

   
       end % Loop for Plus and minus
   
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Index          = [4*(i-2)+1 4*(i-2)+2 4*(i-2)+3 4*(i-2)+4 4*(i-1)+1 4*(i-1)+2 4*(i-1)+3 ];
       
       Grad2_R2(:,:)  =  2*(C(1,2)-C(2,2))*(Grad2_C(1,2,:,:)-Grad2_C(2,2,:,:))...
                       + 2*(C(1,3)-C(2,3))*(Grad2_C(1,3,:,:)-Grad2_C(2,3,:,:));
                   
       Grad2_Kappa    = zeros(2,4*PND.Kite.N,4*PND.Kite.N);            
       for k=1:1:2
           Grad2_Kappa(k,:,:) = -C(k,1)*Grad2_C(k,1,:,:)/Kappa(k); 
       end
                   
       for m=1:1:length(Index)
           Ind2       = Index(m);
           Grad2_R2(:,Ind2)  = Grad2_R2(:,Ind2) + 2*(Grad_C(1,2,Ind2)-Grad_C(2,2,Ind2))*squeeze(Grad_C(1,2,:)-Grad_C(2,2,:))...
                                                + 2*(Grad_C(1,3,Ind2)-Grad_C(2,3,Ind2))*squeeze(Grad_C(1,3,:)-Grad_C(2,3,:));
           for k=1:1:2
               Grad2_Kappa(k,:,Ind2) = Grad2_Kappa(k,:,Ind2) - Grad_C(k,1,Ind2)*squeeze(Grad_C(k,1,:))'/Kappa(k)+ ...
                                                               C(k,1)*squeeze(Grad_C(k,1,:))'*Grad_Kappa(k,Ind2)/Kappa(k)^2; 
           end                                        
       end          
          
      
       Grad2_S1(:,:)   =  2*(Kappa(1)+Kappa(2))*squeeze(Grad2_Kappa(1,:,:) + Grad2_Kappa(2,:,:)) - Grad2_R2(:,:);
       Grad2_S2(:,:)   = -2*(Kappa(1)-Kappa(2))*squeeze(Grad2_Kappa(1,:,:) - Grad2_Kappa(2,:,:)) + Grad2_R2(:,:);

       for m=1:1:length(Index)
               Ind2       = Index(m);
               Grad2_S1(:,Ind2) = Grad2_S1(:,Ind2) + 2*(Grad_Kappa(1,Ind2)+Grad_Kappa(2,Ind2))*squeeze(Grad_Kappa(1,:) + Grad_Kappa(2,:))';
               Grad2_S2(:,Ind2) = Grad2_S2(:,Ind2) - 2*(Grad_Kappa(1,Ind2)-Grad_Kappa(2,Ind2))*squeeze(Grad_Kappa(1,:) - Grad_Kappa(2,:))';
       end 

       Grad2_Sq(:,:)   = (Aux_S2*Grad2_S1 + Aux_S1*Grad2_S2)/(2*Aux_Sq);
       for m=1:1:length(Index)
               Ind2                 = Index(m);
               Grad2_Sq(:,Ind2)   =  Grad2_Sq(:,Ind2)+ (Grad_S2(Ind2)*Grad_S1 + Grad_S1(Ind2)*Grad_S2)'/(2*Aux_Sq)...
                                                     - (Aux_S2*Grad_S1 +  Aux_S1*Grad_S2)'/(2*Aux_Sq^2)*Grad_Sq(Ind2);
       end

       % Do First part of the Gradient of the gradient of zeta
       Aux0 = zeros(4*PND.Kite.N,4*PND.Kite.N);

       Aux0(:,:)  =        squeeze((Grad2_C(1,2,:,:)+Grad2_C(2,2,:,:)))/2;

       Aux0(:,:)  = Aux0(:,:) + (Kappa(1)^2-Kappa(2)^2)/(2*R2)*squeeze((Grad2_C(2,2,:,:) - Grad2_C(1,2,:,:)));
       Aux0(:,:)  = Aux0(:,:) + squeeze(Kappa(1)*Grad2_Kappa(1,:,:) - Kappa(2)*Grad2_Kappa(2,:,:))*(C(2,2)-C(1,2))/R2;
       Aux0(:,:)  = Aux0(:,:) - (Kappa(1)^2-Kappa(2)^2)/(2*R2^2)*(C(2,2)-C(1,2))*squeeze(Grad2_R2(:,:));

       Aux0(:,:)  = Aux0(:,:) +  Grad2_Sq(:,:)/(2*R2)*(C(2,3)-C(1,3));
       Aux0(:,:)  = Aux0(:,:) +  Aux_Sq/(2*R2)*squeeze((Grad2_C(2,3,:,:)-Grad2_C(1,3,:,:)))';
       Aux0(:,:)  = Aux0(:,:) -  Aux_Sq*(C(2,3) - C(1,3))/(2*R2^2)*squeeze(Grad2_R2(:,:));

       % Do First part of the Gradient of the gradient of xi
       Aux1 = zeros(4*PND.Kite.N,4*PND.Kite.N);

       Aux1(:,:)  =        squeeze((Grad2_C(1,3,:,:)+Grad2_C(2,3,:,:)))/2;

       Aux1(:,:)  = Aux1(:,:) + (Kappa(1)^2-Kappa(2)^2)/(2*R2)*squeeze((Grad2_C(2,3,:,:) - Grad2_C(1,3,:,:)));
       Aux1(:,:)  = Aux1(:,:) + squeeze(Kappa(1)*Grad2_Kappa(1,:,:) - Kappa(2)*Grad2_Kappa(2,:,:))*(C(2,3)-C(1,3))/R2;
       Aux1(:,:)  = Aux1(:,:) - (Kappa(1)^2-Kappa(2)^2)/(2*R2^2)*(C(2,3)-C(1,3))*squeeze(Grad2_R2(:,:));

       Aux1(:,:)  = Aux1(:,:) +  Grad2_Sq(:,:)/(2*R2)*(C(1,2)-C(2,2));
       Aux1(:,:)  = Aux1(:,:) +  Aux_Sq/(2*R2)*squeeze((Grad2_C(1,2,:,:)-Grad2_C(2,2,:,:)))';
       Aux1(:,:)  = Aux1(:,:) -  Aux_Sq*(C(1,2) - C(2,2))/(2*R2^2)*squeeze(Grad2_R2(:,:));

       for m=1:1:length(Index)

                Ind2        = Index(m);

                % Do Second part of the Gradient of the gradient of zeta

                Aux0(:,Ind2)   = Aux0(:,Ind2)  + 2*(Kappa(1)*Grad_Kappa(1,Ind2)-Kappa(2)*Grad_Kappa(2,Ind2))/(2*R2)*squeeze((Grad_C(2,2,:) - Grad_C(1,2,:)));
                Aux0(:,Ind2)   = Aux0(:,Ind2)  - (Kappa(1)^2-Kappa(2)^2)/(2*R2^2)*Grad_R2(Ind2)*squeeze((Grad_C(2,2,:) - Grad_C(1,2,:)));

                Aux0(:,Ind2)   = Aux0(:,Ind2)  + (Grad_Kappa(1,Ind2)*Grad_Kappa(1,:) - Grad_Kappa(2,Ind2)*Grad_Kappa(2,:))'*(C(2,2)-C(1,2))/R2;
                Aux0(:,Ind2)   = Aux0(:,Ind2)  + (Kappa(1)*Grad_Kappa(1,:) - Kappa(2)*Grad_Kappa(2,:))'*(Grad_C(2,2,Ind2)-Grad_C(1,2,Ind2))/R2;
                Aux0(:,Ind2)   = Aux0(:,Ind2)  - (Kappa(1)*Grad_Kappa(1,:) - Kappa(2)*Grad_Kappa(2,:))'*(C(2,2)-C(1,2))/(R2^2)*Grad_R2(Ind2);

                Aux0(:,Ind2)   = Aux0(:,Ind2)  - 2*(Kappa(1)*Grad_Kappa(1,Ind2)-Kappa(2)*Grad_Kappa(2,Ind2))/(2*R2^2)*(C(2,2)-C(1,2))*Grad_R2(1,:)';
                Aux0(:,Ind2)   = Aux0(:,Ind2)  - (Kappa(1)^2-Kappa(2)^2)/(2*R2^2)*(Grad_C(2,2,Ind2)-Grad_C(1,2,Ind2))*Grad_R2(1,:)';
                Aux0(:,Ind2)   = Aux0(:,Ind2)  + (Kappa(1)^2-Kappa(2)^2)/(R2^3)*(C(2,2)-C(1,2))*Grad_R2(1,:)'*Grad_R2(1,Ind2);

                Aux0(:,Ind2)   = Aux0(:,Ind2) +  Grad_Sq(1,:)'/(2*R2)*(Grad_C(2,3,Ind2)-Grad_C(1,3,Ind2));
                Aux0(:,Ind2)   = Aux0(:,Ind2) -  Grad_Sq(1,:)'/(2*R2^2)*(C(2,3)-C(1,3))*Grad_R2(Ind2);

                Aux0(:,Ind2)   = Aux0(:,Ind2) +  Grad_Sq(Ind2)/(2*R2)*squeeze((Grad_C(2,3,:)-Grad_C(1,3,:)));
                Aux0(:,Ind2)   = Aux0(:,Ind2) -  Aux_Sq/(2*R2^2)*Grad_R2(Ind2)*squeeze((Grad_C(2,3,:)-Grad_C(1,3,:)));

                Aux0(:,Ind2)   = Aux0(:,Ind2) -  Grad_Sq(Ind2)*(C(2,3) - C(1,3))/(2*R2^2)*Grad_R2(1,:)';
                Aux0(:,Ind2)   = Aux0(:,Ind2) -  Aux_Sq*(Grad_C(2,3,Ind2) - Grad_C(1,3,Ind2))/(2*R2^2)*Grad_R2(1,:)';
                Aux0(:,Ind2)   = Aux0(:,Ind2) +  Aux_Sq*(C(2,3) - C(1,3))/(R2^3)*Grad_R2(1,:)'*Grad_R2(1,Ind2);



                % Do Second part of the Gradient of the gradient of xi

                Aux1(:,Ind2)   = Aux1(:,Ind2)  + 2*(Kappa(1)*Grad_Kappa(1,Ind2)-Kappa(2)*Grad_Kappa(2,Ind2))/(2*R2)*squeeze((Grad_C(2,3,:) - Grad_C(1,3,:)));
                Aux1(:,Ind2)   = Aux1(:,Ind2)  - (Kappa(1)^2-Kappa(2)^2)/(2*R2^2)*Grad_R2(Ind2)*squeeze((Grad_C(2,3,:) - Grad_C(1,3,:)));

                Aux1(:,Ind2)   = Aux1(:,Ind2)  + (Grad_Kappa(1,Ind2)*Grad_Kappa(1,:) - Grad_Kappa(2,Ind2)*Grad_Kappa(2,:))'*(C(2,3)-C(1,3))/R2;
                Aux1(:,Ind2)   = Aux1(:,Ind2)  + (Kappa(1)*Grad_Kappa(1,:) - Kappa(2)*Grad_Kappa(2,:))'*(Grad_C(2,3,Ind2)-Grad_C(1,3,Ind2))/R2;
                Aux1(:,Ind2)   = Aux1(:,Ind2)  - (Kappa(1)*Grad_Kappa(1,:) - Kappa(2)*Grad_Kappa(2,:))'*(C(2,3)-C(1,3))/(R2^2)*Grad_R2(Ind2);

                Aux1(:,Ind2)   = Aux1(:,Ind2)  - 2*(Kappa(1)*Grad_Kappa(1,Ind2)-Kappa(2)*Grad_Kappa(2,Ind2))/(2*R2^2)*(C(2,3)-C(1,3))*Grad_R2(1,:)';
                Aux1(:,Ind2)   = Aux1(:,Ind2)  - (Kappa(1)^2-Kappa(2)^2)/(2*R2^2)*(Grad_C(2,3,Ind2)-Grad_C(1,3,Ind2))*Grad_R2(1,:)';
                Aux1(:,Ind2)   = Aux1(:,Ind2)  + (Kappa(1)^2-Kappa(2)^2)/(R2^3)*(C(2,3)-C(1,3))*Grad_R2(1,:)'*Grad_R2(1,Ind2);

                Aux1(:,Ind2)   = Aux1(:,Ind2) +  Grad_Sq(1,:)'/(2*R2)*(Grad_C(1,2,Ind2)-Grad_C(2,2,Ind2));
                Aux1(:,Ind2)   = Aux1(:,Ind2) -  Grad_Sq(1,:)'/(2*R2^2)*(C(1,2)-C(2,2))*Grad_R2(Ind2);

                Aux1(:,Ind2)   = Aux1(:,Ind2) +  Grad_Sq(Ind2)/(2*R2)*squeeze((Grad_C(1,2,:)-Grad_C(2,2,:)));
                Aux1(:,Ind2)   = Aux1(:,Ind2) -  Aux_Sq/(2*R2^2)*Grad_R2(Ind2)*squeeze((Grad_C(1,2,:)-Grad_C(2,2,:)));

                Aux1(:,Ind2)   = Aux1(:,Ind2) -  Grad_Sq(Ind2)*(C(1,2) - C(2,2))/(2*R2^2)*Grad_R2(1,:)';
                Aux1(:,Ind2)   = Aux1(:,Ind2) -  Aux_Sq*(Grad_C(1,2,Ind2) - Grad_C(2,2,Ind2))/(2*R2^2)*Grad_R2(1,:)';
                Aux1(:,Ind2)   = Aux1(:,Ind2) +  Aux_Sq*(C(1,2) - C(2,2))/(R2^3)*Grad_R2(1,:)'*Grad_R2(1,Ind2);

       end

       Grad2_ZX(1,:,:,i) = Aux0;
       Grad2_ZX(2,:,:,i) = Aux1;
      
       
     
       
         
end % Kite Loop



end