function [rQ rR_Edge R_KE rR vR omegaR FA_R QR rK vK omegaK FA_K MA_K MMC_K alfa_K beta_K QK rG vG omegaG FA_G MA_G MMC_G QG DF RHS]=  Fun_ODE_Lag_Full_Output_KF(t,xs_amp,xc_amp,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Right-Hand-Side and full outputs of Lagrangian Equations       %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%     t                       -> normalized time                          %%
%     Extended state vector   -> xs_amp = [xs xs_p]                       %%
%     Extended control vector -> xc_amp                                   %%
%     PND                     -> Dimensionless parameters                 %%
% Outputs                                                                 %%
%     rQ     -> Point Q vector position  (SE components)                  %%
%     R_KE   -> SK-SE rotation matrix                                     %%
%                                                                         %%
%     rR     -> Rod position vectors (SE components)                      %%
%     vR     -> Rod velocity vectors (SE components)                      %%
%     omegaR -> Rod angular velocity (SR components)                      %%
%     FA_R   -> Aerodynamic force upon the rods (SE components)           %%
%     Q_R    -> Rod Generalized forces                                    %%
%                                                                         %%
%     rK     -> Kite position vectors (SE components)                     %%
%     vK     -> Kite velocity vectors (SE components)                     %%
%     omegaK -> Kite angular velocity (SK components)                     %%
%     FA_K   -> Aerodynamic force upon the kite (SK components)           %%
%     MA_K   -> Aerodynamic moment upon the rods (SK components)          %%
%     MMC_K  -> Motor controller torque  (SK components)                  %%
%     alfa_K -> angle of attack                                           %%
%     beta_K -> sideslip angle                                            %%
%     Q_K    -> Kite Generalized forces                                   %%
%                                                                         %%
%     rG     -> Rotor position vectors (SE components)                    %%
%     vG     -> Rotor velocity vectors (SE components)                    %%
%     omegaG -> Rotor angular velocity (SK components)                    %%
%     FA_G   -> Aerodynamic force upon the kite (SK components)           %%
%     MA_G   -> Aerodynamic moment upon the rods (SK components)          %%
%     MMC_G  -> Motor controller torque  (SK components)                  %%
%     Q_G    -> Rotor Generalized forces                                  %%
%                                                                         %%
%     DF     -> Righthand side of the first order Equations               %%
%     RHS    -> Righthand side of the second order equations              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Important Dimensions
Nvar                          = 2*PND.Num.N+3;              % Number of variables 
Nvar_p                        = 2*PND.Num.N+3+PND.Gen.Num;  % Number of variables (dot) 
Nc                            = 4+PND.Gen.Num+3;            % Number of control parameters
Nc0                           = 4;           % Number of control variables affecting the kinematics
% Recover State and control vector in structured form
[xs_str xs_p_str]             = From_xs2Var_KF(xs_amp,PND);
[xc_str xc_p_str xc_pp_str]   = From_xc2Var_KF(xc_amp,PND);
% Recover the State and control vector in column form
xs_p                          = xs_amp(Nvar+1:Nvar+Nvar_p,1);
xc                            = xc_amp(0*Nc+1:1*Nc,1);
xc_p                          = xc_amp(1*Nc+1:2*Nc,1);
xc_pp                         = xc_amp(2*Nc+1:3*Nc,1);

%% Check if we need to compute dH/dt
if length(xs_amp)>Nvar+Nvar_p
    Compute_dHdt = 1;
    Flag_xc      = 1;
else
    Compute_dHdt = 0;
    Flag_xc      = 0;
end

% Basic Kinetic Matrices
[R_KE Grad_R_KE ]               = Matrix_R_KE_KF(xs_str);
[SR CR OmR SK CK OmK SG CG OmG] = Matrix_V_KF(xs_str,xc_str,R_KE,PND);
  
% Compute matrices of the kinematics terms
[Ms Msc Mc]                             = Matrix_M_KF(SK,CK,OmK,SR,CR,OmR,SG,CG,xs_str,xc_str,PND);
[Ms_xs Ms_xc Msc_xs Msc_xc Mc_xs Mc_xc] = Grad_Matrix_M_KF(R_KE,Grad_R_KE,SR,CR,OmR,SK,CK,OmK,SG,CG,xs_str,xc_str,Flag_xc,PND);

% Compute relevant position vector (Earth frame components)   
[rQ rR rR_Edge rK rG ]     = Compute_Positions_KF(xs_str,xc_str,R_KE,PND);
% Compute velocities and angular velocities
vK           = SK*xs_p+CK*xc_p(1:Nc0,1); % Kite velocity components in SE
omegaK       = OmK*xs_p;                 % Kite angular velocity components in SK
for i=1:1:PND.Num.N 
   vR(:,i)      = squeeze(SR(:,:,i))*xs_p +  CR(:,:,i)*xc_p(1:Nc0,1); % Rod velocity components in the Earth Frame
   omegaR(:,i)  = squeeze(OmR(:,:,i))*xs_p;                           % Rod angular velocity in the Rod frame
end

if PND.Gen.Num==0
    vG = zeros(3,1);
    omegaG = zeros(3,1);
else
    for i=1:1:PND.Gen.Num
        vG(:,i)        = SG(:,:,i)*xs_p+CG(:,:,i)*xc_p(1:Nc0); % Generator velocities components in SE
        omegaG(:,i)    = squeeze(OmG(:,:,i))*xs_p;             % Generator angular velocity in SK                          
    end
end

% Compute the gradient of the potential
U_xs = Grad_U_KF(xs_str,xc_str,PND);

% Generalized Forces 
[QK FA_K MA_K MMC_K alfa_K beta_K] = Generalized_QK_KF(t,rK,vK,omegaK,R_KE,SK,OmK,xc,PND);
[QR FA_R]                          = Generalized_QR_KF(t,rR,vR,SR,xs_str,xc_str,PND);
if PND.Gen.Num==0
    QG    = zeros(Nvar_p,1);
    FA_G  = zeros(3,1);
    MA_G  = zeros(3,1);
    MMC_G = zeros(3,1);
else
    [QG FA_G MA_G MMC_G]           = Generalized_QG_KF(t,rG,vG,SG,OmG,R_KE,xc,PND);
end
Q                              = QK+QR+QG;
% Compute the RHS
AUX2 = zeros(Nvar_p,Nvar_p);
AUX3 = zeros(Nvar_p,Nc0);
AUX4 = zeros(Nvar_p,Nvar_p);
AUX5 = zeros(Nvar_p,Nc0);

RHS6 = zeros(Nvar_p,1);
RHS7 = zeros(Nvar_p,1);
RHS8 = zeros(Nvar_p,1);

for i=1:1:Nvar
    AUX2 = AUX2+xs_p(i)*Ms_xs(:,:,i);
    AUX3 = AUX3+xs_p(i)*Msc_xs(:,:,i);
    if i<=Nc0
        AUX4 = AUX4+xc_p(i)*Ms_xc(:,:,i);
        AUX5 = AUX5+xc_p(i)*Msc_xc(:,:,i);
    end
    RHS6(i,1) =  0.5*xs_p'*Ms_xs(:,:,i)*xs_p;
    RHS7(i,1) =  xs_p'*Msc_xs(:,:,i)*xc_p(1:Nc0,1);
    RHS8(i,1) =  0.5*xc_p(1:Nc0,1)'*Mc_xs(:,:,i)*xc_p(1:Nc0,1);
end


RHS1  = -Msc*xc_pp(1:Nc0,1);

RHS2  = -AUX2*xs_p;
RHS3  = -AUX3*xc_p(1:Nc0,1);
RHS4  = -AUX4*xs_p;
RHS5  = -AUX5*xc_p(1:Nc0,1);

RHS = Q+RHS1+RHS2+RHS3+RHS4+RHS5+RHS6+RHS7+RHS8-[U_xs;zeros(PND.Gen.Num,1)];

% Ordinary diferential Equations
DF(1:Nvar,1)             = xs_p(1:Nvar,1);
DF(Nvar+1:Nvar+Nvar_p,1) = Ms\RHS;

if Compute_dHdt==1
     Lt                  = Compute_dL_dt_KF(Msc,Mc,Ms_xc,Msc_xc,Mc_xc,xs_amp,xc_amp,PND);
     DF(Nvar+Nvar_p+1,1) = xs_p'*Q-Lt;
end

end