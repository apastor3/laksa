function [xc RBE rK vK omegaK FA_K MA_K alfa beta Ups Ups_xs Phi Phi_xs Q DF RHS]=  Fun_ODE_Full_Output_KT(t,xs_amp,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and  Jose A. Serrano-Iglesia           %
% Language  : Matlab                                                         %
% Synopsis  : RHS with Full Outputs                                          %
% Copyright : Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%     t                       -> normalized time                          %%
%     Extended state vector   -> xs_amp = [xs xs_p]                       %%
%     PND                     -> Dimensionless parameters                 %%
% Outputs                                                                 %%
%     R_BE   -> SB-SE rotation matrices                                   %%
%     rK     -> Kite position vectors (SE components)                     %%
%     vK     -> Kite velocity vectors (SE components)                     %%
%     omegaK -> Kite angular velocity (SB components)                     %%
%     FA_K   -> Aerodynamic forces upon the kites (SB components)         %%
%     MA_K   -> Aerodynamic moments upon the kites(SB components)         %%
%     alfa   -> angles of attack                                          %%
%     beta   -> sideslip angles                                           %%
%     Q      -> Generalized forces                                        %%
%     DF     -> Righthand side of the first order Equations               %%
%     RHS    -> Righthand side of the second order equations              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Important Dimensions
NK                   = PND.Kite.N;
Nvar                 = 4*NK;
% Recover the state vector
xs                   = xs_amp(1:Nvar,1);
xs_p                 = xs_amp(Nvar+1:2*Nvar,1);
% Recover the control vector
xc                   = Fun_Control_KT(t,PND);

% Rotation matrices
[R1E R1E_xs R1E_xs2] = Fun_R1E_KT(xs,PND);
[R21 R21_xs R21_xs2] = Fun_R21_KT(xs,PND);
[RB2 RB2_xs RB2_xs2] = Fun_RB2_KT(xs,PND);


% Position vector
[ZX Grad_ZX Grad2_ZX] = Fun_ZX_KT(R1E,R21,RB2,R1E_xs,R21_xs,RB2_xs,R1E_xs2,R21_xs2,RB2_xs2,PND);
rK                    = Position_KT(R1E,R21,RB2,ZX,PND); 
% Kinematic matrices
[Phi Phi_xs]  = Fun_Omega_KT(xs,PND);
[Ups Ups_xs]  = Fun_Upsilon_KT(xs,ZX,Grad_ZX,Grad2_ZX,PND);


for i=1:1:PND.Kite.N 
   RBE(:,:,i)     = squeeze(RB2(:,:,i))*squeeze(R21(:,:,i))*squeeze(R1E(:,:,i));
   omegaK(:,i)    = squeeze(Phi(:,:,i))*xs_p;     % Kite angular velocity components in SB_j
   vK(:,i)        = squeeze(Ups(:,:,i))*xs_p;     % Kite velocity components in SE
   vw(:,i)        = Fun_Wind(t,rK(:,i),PND);      % Wind velocity components in SE
   
   vw_b(:,i)      = squeeze(RBE(:,:,i))*vw(:,i);  % Wind velocity components in SB
   vK_b(:,i)      = squeeze(RBE(:,:,i))*vK(:,i);  % Kite velocity components in SB
end  
% Compute matrices of the kinematics terms
M    = Matrix_M_KT(Ups,Phi,PND);
M_xs = Grad_Matrix_M_KT(Ups,Ups_xs,Phi,Phi_xs,PND);
    
% Compute the gradient of the potential
[U U_xs] = Fun_Potential_KT(xs,ZX,Grad_ZX,PND);

% Generalized Forces 
[Q FA_K MA_K alfa beta] = Fun_Q_KT(xc,vw_b,vK_b,omegaK,Ups,Phi,RBE,PND);

% Compute the RHS
AUX1 = zeros(Nvar,Nvar);
RHS1 = zeros(Nvar,1);
RHS2 = zeros(Nvar,1);

for i=1:1:Nvar
     AUX1 = AUX1+xs_p(i)*M_xs(:,:,i);
     RHS2(i,1) =  0.5*xs_p'*M_xs(:,:,i)*xs_p;
end
RHS1  = -AUX1*xs_p;
RHS   = Q+RHS1+RHS2-U_xs;

% Ordinary diferential Equations
DF(1:Nvar,1)        = xs_p(1:Nvar,1);
DF(Nvar+1:2*Nvar,1) = M\RHS;



end