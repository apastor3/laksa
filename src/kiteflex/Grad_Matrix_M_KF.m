function [Ms_xs Ms_xc Msc_xs Msc_xc Mc_xs Mc_xc] = Grad_Matrix_M_KF(R_KE,Grad_R_KE,SR,CR,OmR,SK,CK,OmK,SG,CG,xs,xc,Flag_xc,PND)
                                                                
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Gradients of Matrix M                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:                                                                   %%
%        R_KE and Grad_R_KE -> Rotation matrix and its gradient             %%
%        SR,CR,OmR          -> Rod matrices                                 %%
%        SK,CK,OmK          -> Kite matrices                                %%
%        SG,CG              -> Generator matrices                           %%
%        xs,xc              -> State anc ontrol vector                      %%
%        Flag_xc            -> = 0 Mc_xc is not computed (Mc_xc is not      %%
%                               needed for the Lagrange equations of motion %%  
%                              = 1 Mc_xc is  computed  (Mc_xc is required   %% 
%                                  to compute Lt and check the energy)      %%
%                              = 2 Ms_xc, Msc_xc and Mc_xc are not computed %% 
%                                  (Hamilton equation)                      %% 
%         PND               -> Dimensionless parameter                      %%
% Outputs                                                                   %%
%         Ms_xs             -> Gradient of Ms with respect to xs            %%
%         Ms_xc             -> Gradient of Ms with respect to xc            %%
%         Msc_xs            -> Gradient of Msc with respect to xs           %%
%         Msc_xc            -> Gradient of Msc with respect to xc           %%
%         Mc_xs             -> Gradient of Mc with respect to xs            %%
%         Mc_xc             -> Gradient of Mc with respect to xc            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensions
NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators

Nvar   = 2*NR+3;
Nvar_p = 2*NR+3+NG;
Nc     = 4+NG+3;       % Total number of control variables 
Nc0    = 4;            % Control variables affecting the kinematics

% Kite Parameters
iK        = PND.Kite.ik;
% Tether Parameter
sigmaR    = PND.Tether.Sigma;
iR        = PND.Tether.Ups;

if NG>0
    % Generator Parameters
    iG      = PND.Gen.iota; % Generator tensor of inertia (divided by Sigma*lG^2)
    sigmaG  = PND.Gen.Sigma;
    lG      = PND.Gen.l;
    nu      = PND.Gen.nu;
end

% Compute the gradients
[SR_xs SR_xc CR_xs OmR_xs]       = Grad_Rod_KF(xs,xc,PND); 
[SK_xs SK_xc CK_xs CK_xc OmK_xs] = Grad_Kite_KF(xs,xc,R_KE,Grad_R_KE,PND); 
[SG_xs SG_xc CG_xs CG_xc ]       = Grad_Rotors_KF(xs,xc,R_KE,Grad_R_KE,SK_xs,SK_xc,CK_xs,CK_xc,PND);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gradients with respect to xs  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ms_xs  = zeros(Nvar_p,Nvar_p,Nvar);
Msc_xs = zeros(Nvar_p,Nc0,Nvar);
Mc_xs  = zeros(Nc0,Nc0,Nvar);

% Compute the contribution from omega_G * iG * omega_G
theta = xs.theta;
psi   = xs.psi;
phi   = xs.phi;
Grad_MSG_xs = zeros(Nvar_p,Nvar_p,NG,Nvar);
B           = zeros(3,3,NG);
for j=1:1:NG
     
  % Terms lambda_dot by (theta_dot, psi_dot, and phi_dot )
  % Derivative with respect to theta
  Grad_MSG_xs(2*NR+3+j,2*NR+1:2*NR+3,j,2*NR+1) =  iG(1,1)*[0 -cos(nu(j))*cos(theta)+sin(nu(j))*sin(theta)*cos(phi) 0];
  Grad_MSG_xs(2*NR+1:2*NR+3,2*NR+3+j,j,2*NR+1) =  squeeze(Grad_MSG_xs(2*NR+3+j,2*NR+1:2*NR+3,j,2*NR+1))';
  % Derivative with respect to phi
  Grad_MSG_xs(2*NR+3+j,2*NR+1:2*NR+3,j,2*NR+3) =  iG(1,1)*[sin(nu(j))*cos(phi)  sin(nu(j))*cos(theta)*sin(phi) 0];
  Grad_MSG_xs(2*NR+1:2*NR+3,2*NR+3+j,j,2*NR+3) =  squeeze(Grad_MSG_xs(2*NR+3+j,2*NR+1:2*NR+3,j,2*NR+3))';

  
  % Terms (theta_dot, psi_dot, and phi_dot )*(theta_dot, psi_dot, and phi_dot )
  B(1,1:3,j)  = [ iG(1,1)*(cos(nu(j)))^2+iG(2,2)*(sin(nu(j)))^2                  0      -iG(1,1)*sin(nu(j))*cos(nu(j))+iG(2,2)*sin(nu(j))*cos(nu(j)) ];
  B(2,1:3,j)  = [0                                                            iG(2,2)               0 ];
  B(3,1:3,j)  = [-iG(1,1)*sin(nu(j))*cos(nu(j))+iG(2,2)*sin(nu(j))*cos(nu(j))    0       iG(1,1)*(sin(nu(j)))^2+iG(2,2)*(cos(nu(j)))^2];
  % Ther gradient Grad_MSG_xs is not finishes (see below)
end

for i=1:1:Nvar
    % Kite contributions
    Ms_xs(:,:,i)   = squeeze(SK_xs(:,:,i))'*SK + SK'*squeeze(SK_xs(:,:,i));
    Ms_xs(:,:,i)   = Ms_xs(:,:,i) + squeeze(OmK_xs(:,:,i))'*iK*OmK + OmK'*iK*squeeze(OmK_xs(:,:,i));

    Msc_xs(:,:,i)  = squeeze(SK_xs(:,:,i))'*CK + SK'*CK_xs(:,:,i);
    
    Mc_xs(:,:,i)  = squeeze(CK_xs(:,:,i))'*CK + CK'*CK_xs(:,:,i);
    
    % Bars contributions
    for j=1:1:NR
        Aux1 = (squeeze(SR_xs(:,:,j,i)))'*squeeze(SR(:,:,j)) + (squeeze(SR(:,:,j)))'*squeeze(SR_xs(:,:,j,i));
        Aux2 = (squeeze(OmR_xs(:,:,j,i)))'*iR*squeeze(OmR(:,:,j))+ (squeeze(OmR(:,:,j)))'*iR*squeeze(OmR_xs(:,:,j,i));
        Ms_xs(:,:,i)   = Ms_xs(:,:,i) + sigmaR*xc.lt*(Aux1+xc.lt^2*Aux2);

        Aux1 = (squeeze(SR_xs(:,:,j,i)))'*(squeeze(CR(:,:,j)))+(squeeze(SR(:,:,j)))'*(squeeze(CR_xs(:,:,j,i)));
        Msc_xs(:,:,i)  = Msc_xs(:,:,i) + sigmaR*xc.lt*Aux1;

        Aux1 = (squeeze(CR_xs(:,:,j,i)))'*(squeeze(CR(:,:,j))) + (squeeze(CR(:,:,j)))'*(squeeze(CR_xs(:,:,j,i)));
        Mc_xs(:,:,i)   = Mc_xs(:,:,i)  + sigmaR*xc.lt*Aux1;  
    end   
    
    % Generator contributions
    for j=1:1:NG
        Aux1 = (squeeze(SG_xs(:,:,j,i)))'*squeeze(SG(:,:,j)) + (squeeze(SG(:,:,j)))'*squeeze(SG_xs(:,:,j,i));
            
        %% The gradient is completed here    
        Grad_MSG_xs(2*NR+1:2*NR+3,2*NR+1:2*NR+3,j,i) = squeeze(OmK_xs(:,2*NR+1:2*NR+3,i))'*squeeze(B(:,:,j))*OmK(:,2*NR+1:2*NR+3) + OmK(:,2*NR+1:2*NR+3)'*squeeze(B(:,:,j))*squeeze(OmK_xs(:,2*NR+1:2*NR+3,i));
        
        Ms_xs(:,:,i)   = Ms_xs(:,:,i) + sigmaG*(Aux1+lG^2* Grad_MSG_xs(:,:,j,i));

        Aux1 = (squeeze(SG_xs(:,:,j,i)))'*(squeeze(CG(:,:,j)))+(squeeze(SG(:,:,j)))'*(squeeze(CG_xs(:,:,j,i)));
        Msc_xs(:,:,i)  = Msc_xs(:,:,i) + sigmaG*Aux1;

        Aux1 = (squeeze(CG_xs(:,:,j,i)))'*(squeeze(CG(:,:,j))) + (squeeze(CG(:,:,j)))'*(squeeze(CG_xs(:,:,j,i)));
        Mc_xs(:,:,i)   = Mc_xs(:,:,i)  + sigmaG*Aux1;  
    end   
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gradients with respect to xc  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ms_xc        = zeros(Nvar_p,Nvar_p,Nc0);
Msc_xc       = zeros(Nvar_p,Nc0,Nc0);
Mc_xc        = zeros(Nc0,Nc0,Nc0);

if Flag_xc<2
    for i=1:1:Nc0
        % Kite contributions
        Ms_xc(:,:,i)   = squeeze(SK_xc(:,:,i))'*SK + SK'*squeeze(SK_xc(:,:,i));
        Msc_xc(:,:,i)  = squeeze(SK_xc(:,:,i))'*CK + SK'*CK_xc(:,:,i);
        
        if Flag_xc==1
            Mc_xc(:,:,i)   = squeeze(CK_xc(:,:,i))'*CK + CK'*CK_xc(:,:,i);
        end
         % Bars contributions
        for j=1:1:NR
            Aux1 = (squeeze(SR_xc(:,:,j,i)))'*squeeze(SR(:,:,j)) + (squeeze(SR(:,:,j)))'*squeeze(SR_xc(:,:,j,i));
            Ms_xc(:,:,i)   = Ms_xc(:,:,i) + sigmaR*xc.lt*(Aux1);

            if i==1
                Aux1           = (squeeze(SR(:,:,j)))'*squeeze(SR(:,:,j));
                Aux2           = (squeeze(OmR(:,:,j)))'*iR*squeeze(OmR(:,:,j));
                Ms_xc(:,:,i)   = Ms_xc(:,:,i) + sigmaR*(Aux1+3*xc.lt^2*Aux2);
            end

            Aux1           = (squeeze(SR_xc(:,:,j,i)))'*(squeeze(CR(:,:,j)));
            Msc_xc(:,:,i)  = Msc_xc(:,:,i) + sigmaR*xc.lt*Aux1;

            if i==1
                Aux1           = (squeeze(SR(:,:,j)))'*squeeze(CR(:,:,j));
                Msc_xc(:,:,i)  = Msc_xc(:,:,i) + sigmaR*Aux1;
            end

            if Flag_xc==1 &&  i==1    
                    Aux1           = (squeeze(CR(:,:,j)))'*squeeze(CR(:,:,j));
                    Mc_xc(:,:,i)   =  Mc_xc(:,:,i) + sigmaR*Aux1;
            end

        end    
        
        % Generator contributions
        for j=1:1:NG
            Aux1 = (squeeze(SG_xc(:,:,j,i)))'*squeeze(SG(:,:,j)) + (squeeze(SG(:,:,j)))'*squeeze(SG_xc(:,:,j,i));
            Ms_xc(:,:,i)   = Ms_xc(:,:,i) + sigmaG*(Aux1);

            
            Aux1 = (squeeze(SG_xc(:,:,j,i)))'*squeeze(CG(:,:,j)) + (squeeze(SG(:,:,j)))'*squeeze(CG_xc(:,:,j,i));   
            Msc_xc(:,:,i)  = Msc_xc(:,:,i) + sigmaG*Aux1;

            if Flag_xc==1
                Aux1 = (squeeze(CG_xc(:,:,j,i)))'*squeeze(CG(:,:,j)) + (squeeze(CG(:,:,j)))'*squeeze(CG_xc(:,:,j,i));   
                Mc_xc(:,:,i)   =  Mc_xc(:,:,i) + sigmaG*Aux1;
            end
        end    
      
    end
end

end