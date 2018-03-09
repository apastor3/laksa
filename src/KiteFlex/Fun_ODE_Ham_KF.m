function DF = Fun_ODE_KF(t,xham_amp)

global PND
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Right-Hand-Side of Hamilton's equations                        %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
% Input                                                                      %
% t           - > Time                                                       %
% X_Ham_amp   - > Hamiltonian variables X_Ham = [xs p]                       %
% Output                                                                     %
% DF          - > d X_Ham/dt                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Important Dimensions 
Nvar                              = 2*PND.Num.N+3;              % Number of variables 
Nvar_p                            = 2*PND.Num.N+3+PND.Gen.Num;  % Number of variables (dot) 
Nc                                = 4+PND.Gen.Num+3;            % Number of controlparameters
Nc0                               = 4;           % Number of control variables affecting the kinematics
NR                                = PND.Num.N; % Number of Bars
NG                                = PND.Gen.Num;
% Check if we need to compute dH/dt
if length(xham_amp)>Nvar+Nvar_p
    Compute_dHdt = 1;
    Flag_xc      = 1;
else
    Compute_dHdt = 0;
    Flag_xc      = 2;
end
%
% Compute the control
xc_amp                            = Fun_Control_KF(t,xham_amp,PND);
xc                                = xc_amp(0*Nc+1:1*Nc,1);
xc_p                              = xc_amp(1*Nc+1:2*Nc,1);
xc_pp                             = xc_amp(2*Nc+1:3*Nc,1);
% Structure form
[xc_str xc_p_str xc_pp_str]       = From_xc2Var_KF(xc_amp,PND);
% Recover Lagrangian variables
xs_amp(1:Nvar,1)                  = xham_amp(1:Nvar,1);
p                                 = xham_amp(Nvar+1:Nvar+Nvar_p,1);

% Recover state vector in structured form
xs_str.gamma = xs_amp(1:NR,:);
xs_str.varphi = xs_amp(NR+1:2*NR,:); 
xs_str.theta = xs_amp(2*NR+1,:);
xs_str.psi = xs_amp(2*NR+2,:);
xs_str.phi = xs_amp(2*NR+3,:);
% Compute basic matrices
[R_KE Grad_R_KE ]               = Matrix_R_KE_KF(xs_str);
[SR CR OmR SK CK OmK SG CG OmG] = Matrix_V_KF(xs_str,xc_str,R_KE,PND);

[Ms Msc Mc]                      = Matrix_M_KF(SK,CK,OmK,SR,CR,OmR,SG,CG,xs_str,xc_str,PND);

% Compute relevant position vector (Earth frame components)   
[rQ rR rR_Edge rK rG ]     = Compute_Positions_KF(xs_str,xc_str,R_KE,PND);

% Relevant matrices for the Hamiltonian
Hp                       = inv(Ms);
Hcp                      = Msc'*Hp;
%Hc                       = Hcp*Msc-Mc; % This is not need for the equations of motion

xs_p                     = Hp*(p-Msc*xc_p(1:Nc0,1));
% Recover the full state vector in Lagrangian coordinates 
xs_amp(Nvar+1:Nvar+Nvar_p,1)  = xs_p;

% Recover State and control vector in structured form
[xs_str xs_p_str]                 = From_xs2Var_KF(xs_amp,PND);
% Recover the State and control vector in column form
xs_p                              = xs_amp(Nvar+1:Nvar+Nvar_p,1);

% Compute matrices of the kinematics terms
[Ms_xs Ms_xc Msc_xs Msc_xc Mc_xs Mc_xc] = Grad_Matrix_M_KF(R_KE,Grad_R_KE,SR,CR,OmR,SK,CK,OmK,SG,CG,xs_str,xc_str,Flag_xc,PND);

% Compute tensors
Hp_xs  = zeros(Nvar_p,Nvar_p,Nvar);
Hcp_xs = zeros(Nc0,Nvar_p,Nvar);
Hc_xs  = zeros(Nc0,Nc0,Nvar);
for i=1:1:Nvar
    Hp_xs(:,:,i)  = -Hp*Ms_xs(:,:,i)*Hp;
    Hcp_xs(:,:,i) =  squeeze(Msc_xs(:,:,i))'*Hp+Msc'*Hp_xs(:,:,i);
    Hc_xs(:,:,i)  =  Hcp_xs(:,:,i)*Msc+Hcp*Msc_xs(:,:,i)-Mc_xs(:,:,i);
end


% Check Hamiltonian
%[Em Ec Ep H] = Compute_Energy(Sk,Ck,Omk,Si,Ci,Omi,PND.Tether.Sigma,PND.Inertia.ik,PND.Tether.Ups,xs_amp,xc_amp)
%H0 = 0.5*(p'*Hp*p-2*xc_p'*Hcp*p+xc_p'*Hc*xc_p)+Ep 

% Compute the gradient of the potential
U_xs = Grad_U_KF(xs_str,xc_str,PND);

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
Q          = QK+QR+QG;
% Compute the RHS
RHS = zeros(Nvar_p,1);
for i=1:1:Nvar
   RHS(i,1) = p'*Hp_xs(:,:,i)*p-2*xc_p(1:Nc0,1)'*Hcp_xs(:,:,i)*p+xc_p(1:Nc0,1)'*Hc_xs(:,:,i)*xc_p(1:Nc0,1); 
end
% Ordinary diferential Equations
DF(1:Nvar,1)             = Hp(1:Nvar,:)*p-Hcp(:,1:Nvar)'*xc_p(1:Nc0,1);
DF(Nvar+1:Nvar+Nvar_p,1) = -0.5*RHS-[U_xs;zeros(PND.Gen.Num,1)]+Q;

if Compute_dHdt==1
   Lt                  = Compute_dL_dt_KF(Msc,Mc,Ms_xc,Msc_xc,Mc_xc,xs_amp,xc_amp,PND);
   DF(Nvar+Nvar_p+1,1) = xs_p'*Q-Lt;
end


end