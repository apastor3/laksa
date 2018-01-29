function [T rR_Edge rR vR aR omegaR gR FA_R Tension rQ R_KE ...
          rK vK aK euler omegaK gK FA_K FR_K FG_K  MA_K MR_K MG_K MMC_K alfa_K beta_K... 
          rG vG aG omegaG gG FA_G  FK_G MA_G MK_G  MMC_G xc_out xs_target_out Error_C] = Fun_Post_KF(PD,PND,t,xs_amp,Flag_Dim)
            
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Postprocess the results                                        %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: PD        -> Physical parameters of the system                  %%
%         PND       -> Dimensionless parameters                           %%
%         t         -> Dimensionless time                                 %% 
%         xs_amp    -> Extended state vector                              %% 
%         Flag_Dim  -> set 0 for dimensionless outputs                    %%
%                      set 1 for outputs with dimensions                  %%
%                                                                         %%
% Outputs:  Dimensionless or with dimensions depending on Flag_Dim        %%
%          T        -> Time                                               %%
%          rR_Edge  -> Position of rods tips  (SE components)             %%
%          rR       -> Position of rods' centers (SE components)          %%
%          vR       -> Velocity of the rods (SE components)               %%
%          aR       -> Acceleration of the rods (SE components)           %%
%          omegaR   -> Angular velocity of the rods (SR components)       %%
%          gR       -> Angular acceleration of the rods (SR components)   %%
%          FA_R     -> Aerodynamic force upon the rods  (SE components)   %%
%          Tension  -> Tension at the rods  (SE components)               %%
%          rQ       -> Position of vector Q  (SE components)              %%
%          R_KE     -> SK-SE rotation matrix                              %%
%          rK       -> Kite position  (SE components)                     %%
%          vK       -> Kite velocity  (SE components)                     %%
%          aK       -> Kite acceleration  (SE components)                 %%
%          euler    -> Kite Euler angles                                  %%
%          omegaK   -> Kite angular velocity  (SK components)             %%
%          gK       -> Kite angular acceleration  (SK components)         %%
%          FA_K     -> Kite aerodynamic force  (SE components)            %%
%          FR_K     -> Rod force upon the kite  (SE components)           %%
%          FG_K     -> Rotor force upon the kite  (SE components)         %%
%          MA_K     -> Kite aerodynamic torque   (SK components)          %%
%          MR_K     -> Rotor torque upon the kite   (SK components)       %%
%          MG_K     -> Rods torque upon the kite  (SK components)         %%
%          MMC_K    -> Motor Cont. torque upon the kite (SK components)   %%
%          alfa_K   -> Kite angle of attack                               %%
%          beta_K   -> Kite sideslip angles                               %%
%          rG       -> Rotor position  (SE components)                    %%
%          vG       -> Rotor velocity (SE components)                     %%
%          aG       -> Rotor acceleration (SE components)                 %%
%          omegaG   -> Rotor angular velocity   (SK components)           %%
%          gG       -> Rotor angular acceleration  (SK components)        %%
%          FA_G     -> Rotor aerodynamic force (SE components)            %%
%          FK_G     -> Kite forceupon the rotors (SE components)          %%
%          MA_G     -> Rotor Aerodynamic torque   (SK components)         %%
%          MK_G     -> Kite torque upon the rotors  (SK components)       %%
%          MMC_G    -> Motor Cont. torque upon the kite   (SK components) %%
%          xc_out   -> Control vector                                     %%  
%          xs_target_out -> Target state vector for close loop            %%
%          Error_C  -> Error compared with classical formulation          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Important Dimensions
Nvar        = 2*PND.Num.N+3;              % Number of variables 
Nvar_p      = 2*PND.Num.N+3+PND.Gen.Num;  % Number of variables (dot) 
Nc          = 4+PND.Gen.Num+3;            % Number of control parameters
Nc0         = 4;                          % Number of control variables affecting the kinematics
NR          = PND.Num.N;                  % Number of Rods
NG          = PND.Gen.Num;                % Number of Generators
% Compute the control
xc_amp      = Fun_Control_KF(t,xs_amp,PND);
if PND.Control.Type == 5  %% Close loop
    % Target Trajectory
    [xs_target_out xs_p_target_out xs_pp_target_out ]   = Target_Trajectories(t,PND);
else  
   xs_target_out  = zeros(2*NR+3,1);
end

xc_out(1:4+NG+3,1) = xc_amp(1:4+NG+3,1);


         
% Recover State and control vector in structured form
[xs_str xs_p_str]             = From_xs2Var_KF(xs_amp,PND);
[xc_str xc_p_str xc_pp_str]   = From_xc2Var_KF(xc_amp,PND);
%% Time
T = t;
% Recover the State and control vector in column form
xs      = xs_amp(1:Nvar,1);
xs_p    = xs_amp(Nvar+1:Nvar+Nvar_p,1);
xc      = xc_amp(0*Nc+1:1*Nc,1);
xc_p    = xc_amp(1*Nc+1:2*Nc,1);
xc_pp   = xc_amp(2*Nc+1:3*Nc,1);

% Compute positions, velocities, angular velocities, and external forces and torques  
[rQ rR_Edge R_KE rR vR omegaR FA_R QR rK vK omegaK FA_K MA_K MMC_K alfa_K beta_K QK rG vG omegaG FA_G MA_G MMC_G QG DF RHS]=  Fun_ODE_Lag_Full_Output_KF(t,xs_amp,xc_amp,PND);
xs_pp   = DF(Nvar+1:Nvar+Nvar_p,1);
% Rotation matrices
[R_KE Grad_R_KE  ] = Matrix_R_KE_KF(xs_str);
[R_RE ]            = Matrix_R_RE_KF(xs_str,PND);
% Velocity matrices
[SR CR OmR SK CK OmK SG CG OmG] = Matrix_V_KF(xs_str,xc_str,R_KE,PND);
% Compute gradients of the velocity matrices
[SR_xs SR_xc CR_xs OmR_xs]       = Grad_Rod_KF(xs_str,xc_str,PND); 
[SK_xs SK_xc CK_xs CK_xc OmK_xs] = Grad_Kite_KF(xs_str,xc_str,R_KE,Grad_R_KE,PND); 
[SG_xs SG_xc CG_xs CG_xc ]       = Grad_Rotors_KF(xs_str,xc_str,R_KE,Grad_R_KE,SK_xs,SK_xc,CK_xs,CK_xc,PND);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute acceleration of the kite  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aux0_K        = zeros(3,Nvar_p); % (partial WK/partial xs)*dot{xs} 
Aux1_K        = zeros(3,Nvar_p); % (partial SK/partial xs)*dot{xs}
Aux2_K        = zeros(3,Nc0);  % (partial CK/partial xs)*dot{xs}
for k=1:1:Nvar
    Aux0_K = Aux0_K + squeeze(OmK_xs(:,:,k))*xs_p(k);  
    Aux1_K = Aux1_K + squeeze(SK_xs(:,:,k))*xs_p(k);
    Aux2_K = Aux2_K + squeeze(CK_xs(:,:,k))*xs_p(k);
end
Aux3_K = zeros(3,Nvar_p); % (partial SK/partial xc)*dot{xc}
Aux4_K = zeros(3,Nc0);    % (partial CK/partial xc)*dot{xc}
for k=1:1:Nc0
    Aux3_K(:,:) = Aux3_K(:,:) + squeeze(SK_xc(:,:,k))*xc_p(k);
    Aux4_K(:,:) = Aux4_K(:,:) + squeeze(CK_xc(:,:,k))*xc_p(k);
end
% Kite acceleration (SE components)
aK   = SK*xs_pp+(Aux1_K+Aux3_K)*xs_p+CK*xc_pp(1:Nc0,1)+(Aux2_K+Aux4_K)*xc_p(1:Nc0,1);
% Kite angular acceleration (SK components)
gK   = OmK*xs_pp+Aux0_K*xs_p;
% Angular momentum about the center of mass of the kite (SK components)
IK        = PND.Kite.ik;
LK        = IK*omegaK;
% Total force and torque upon the kite due to the Rotors and rod number  NR
MT_K      = IK*gK + cross(omegaK,LK)- MA_K - MMC_K; % Total torque due to the constraints forces (SK Components)
FT_K      = aK - [0 0 1]' - R_KE'*FA_K; % Total constraint force (SE Components)   
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute accelerations of the rotors     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aG   = zeros(3,1);
gG   = zeros(3,1);
MT_G = zeros(3,1); 
FT_G = zeros(3,1);
for i=1:1:NG 
    Aux1_G        = zeros(3,Nvar_p); % (partial SR/partial xs)*dot{xs}
    Aux2_G        = zeros(3,Nc0);  % (partial CR/partial xs)*dot{xs}
    for k=1:1:Nvar  
        Aux1_G = Aux1_G + squeeze(SG_xs(:,:,i,k))*xs_p(k);
        Aux2_G = Aux2_G + squeeze(CG_xs(:,:,i,k))*xs_p(k);
    end
    Aux3_G = zeros(3,Nvar_p); % (partial SK/partial xc)*dot{xc}
    Aux4_G = zeros(3,Nc0);    % (partial CK/partial xc)*dot{xc}
    for k=1:1:Nc0
        Aux3_G(:,:) = Aux3_G(:,:) + squeeze(SG_xc(:,:,i,k))*xc_p(k);
        Aux4_G(:,:) = Aux4_G(:,:) + squeeze(CG_xc(:,:,i,k))*xc_p(k);
    end
    % Rotor acceleration components in SE
    aG(:,i)   = SG(:,:,i)*xs_pp+(Aux1_G+Aux3_G)*xs_p+CG(:,:,i)*xc_pp(1:Nc0,1)+(Aux2_G+Aux4_G)*xc_p(1:Nc0,1);
    % Angular momentum about the center of mass of the rod (SR components)
    nu        = PND.Gen.nu(i); 
    IG        = PND.Gen.Sigma*PND.Gen.l^2*PND.Gen.iota; % Components in the Rotor Frame
    IG        = [IG(1,1)*(cos(nu))^2+IG(2,2)*(sin(nu))^2      0     (IG(2,2)-IG(1,1))*sin(nu)*cos(nu);...
                             0                               IG(2,2)         0;...
                (IG(2,2)-IG(1,1))*sin(nu)*cos(nu)             0      IG(1,1)*(sin(nu))^2+IG(2,2)*(cos(nu))^2]; % Components in the kite frame
    LG        = IG*omegaG(:,i);  % Rotor angular momentum (SK components)
   % Rotor angular acceleration for an observer linked to SG (components in SK)
    ux        = [cos(nu) 0 -sin(nu)]';  % Axis of rotation (components in SK)
    omega_GK  =  xs_p(2*NR+3+i)*ux;      % SR-SK angular velocity (components in SK)
    gG(:,i)   = xs_pp(2*NR+3+i)*ux + gK + cross(omegaK,omega_GK); %(components in SK)
   
    % Total force and torque about OG_j upon the rotors due to the kite
    MT_G(:,i) = IG*gG(:,i)+cross(omegaG(:,i),LG)-MA_G(:,i)-MMC_G(:,i); % Total torque due to the constraint forces (SK Components)
    FT_G(:,i) = PND.Gen.Sigma*(aG(:,i)-[0 0 1]')-R_KE'*FA_G(:,i); % Total constraint force (SE Components)       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute accelerations of the rods %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tensor og inertia of a rod about its center of mass  (SR components)
IR        =   PND.Tether.Sigma*(xc_str.lt)^3*PND.Tether.Ups;
IR_t      = 3*PND.Tether.Sigma*(xc_str.lt)^2*PND.Tether.Ups*xc_p_str.lt;
MassR     =   PND.Tether.Sigma*xc_str.lt;
MassR_t   =   PND.Tether.Sigma*xc_p_str.lt;
for i=1:1:NR
    Aux0_R        = zeros(3,Nvar_p); % (partial WR/partial xs)*dot{xs} 
    Aux1_R        = zeros(3,Nvar_p); % (partial SR/partial xs)*dot{xs}
    Aux2_R        = zeros(3,Nc0);  % (partial CR/partial xs)*dot{xs}
    for k=1:1:Nvar
        Aux0_R = Aux0_R + squeeze(OmR_xs(:,:,i,k))*xs_p(k);  
        Aux1_R = Aux1_R + squeeze(SR_xs(:,:,i,k))*xs_p(k);
        Aux2_R = Aux2_R + squeeze(CR_xs(:,:,i,k))*xs_p(k);
    end
    Aux3_R = zeros(3,Nvar_p); % (partial SR/partial xc)*dot{xc}
    for k=1:1:Nc0
        Aux3_R(:,:) = Aux3_R(:,:) + squeeze(SR_xc(:,:,i,k))*xc_p(k);
    end
    % Rods acceleration components in SE
    aR(:,i)   = SR(:,:,i)*xs_pp+(Aux1_R+Aux3_R)*xs_p+CR(:,:,i)*xc_pp(1:Nc0,1)+Aux2_R*xc_p(1:Nc0,1);
    % Angular acceleration components in the rod frame
    gR(:,i)   = OmR(:,:,i)*xs_pp+Aux0_R*xs_p;
    % Angular momentum about the center of mass of the rod (SR components)
    LR        = IR*omegaR(:,i);
    % Total force and torque upon the rods
    MT_R(:,i) = IR*gR(:,i)+IR_t*omegaR(:,i)+cross(omegaR(:,i),LR); % Total torque due to the constraint forces (SR Components)
    FT_R(:,i) = MassR*(aR(:,i)-[0 0 1]')-FA_R(:,i)+MassR_t*vR(:,i); % Total constraint force (SE Components)       
    
end


%% Compute all the constraints forces and torques
FG_K = zeros(3,1);
MG_K = zeros(3,1);
MK_G = zeros(3,NG);
FK_G = zeros(3,NG);
for i=1:1:NG
   % Moment about OG_j upon the rotor due to the kite (SK Components)
   MK_G(:,i) = MT_G(:,i);   
   % Force upon the rotor due to the kite (SE Components)
   FK_G(:,i) = FT_G(:,i);    
   % Force upong the kite due to the rotor (SE Components)
   FG_K      = FG_K - FT_G(:,i);
   % Moment about OK upon the  kite due to the rotors (SK components)
   MG_K      = MG_K - MK_G(:,i) - cross([PND.Gen.x(i) PND.Gen.y(i) PND.Gen.z(i)]',R_KE*FK_G(:,i)); 
end
% Force upon the kite due to the the last rod (SE components)
FR_K = FT_K - FG_K;
% Torque about OK upon the kite due to the last rod (SK components)
MR_K = cross(xc_str.lb*[cos(xc_str.delta)*cos(xc_str.eta)  cos(xc_str.delta)*sin(xc_str.eta)  sin(xc_str.delta)]',R_KE*FR_K);

% Force at the upper tip (point C+) of rod NR (SE componets)
FCp(:,NR)     = - FR_K;
for i=NR:-1:2
    % Force at the lower tip (point C-) of rod i (SE components)
    FCm(:,i)   =  FT_R(:,i) - FCp(:,i);
    % Force at the upper tip (point C+) of rod i-1 (SE components)
    FCp(:,i-1) =  -FCm(:,i);
end
FCm(:,1)          =  FT_R(:,1) - FCp(:,1);
Tension(:,1)      = -FCm(:,1);
Tension(:,2:NR+1) =  FCp(:,1:NR);

% Torques about OR_j due to the constraint forces  (SR components)
MR_R = zeros(3,NR);
for i=1:1:NR
    MR_R(:,i) = cross([-xc_str.lt/2 0 0 ]',R_RE(:,:,i)*FCp(:,i)) + cross([xc_str.lt/2 0 0 ]',R_RE(:,:,i)*FCm(:,i)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check for consistency with classical mechanics formulation  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Angular momentum equation of the kite
Error_K = max(abs(MT_K - (MR_K + MG_K)));

% Angular momentum equation of the rods
for i=1:1:NR
    Error_R = max(max(abs(MT_R(:,i)-MR_R(:,i))));
end
Error_C = max([Error_K Error_R]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Kinematics calculations                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Euler angles
euler(2,1)  = xs_amp(2*NR+1);%   -asin(R_KE(1,3));               % Theta (pithc)
euler(1,1)  = xs_amp(2*NR+2);% asin(R_KE(1,2)/cos(euler(2,1)));  % Psi   (yaw)
euler(3,1)  = xs_amp(2*NR+3);% asin(R_KE(2,3)/cos(euler(2,1)));  % Phi   (roll)
 
if Flag_Dim==1  % Put the outputs with dimensions
   
    % Characteristic Values
    L_ast     = PD.Tether.L;                          % Distance
    V_ast     = sqrt(PD.Tether.L*PD.Env.g);           % Velocity
    A_ast     = PD.Env.g;                             % Acceleration
    Omega_ast = sqrt(PD.Env.g/PD.Tether.L)*3600/2*pi; % angular velocity
    Gamma_ast = PD.Env.g/PD.Tether.L;                 % angular acceleration
     
    T_ast     = sqrt(PD.Tether.L/PD.Env.g);           % Time 
    F_ast     = PD.Inertia.m*PD.Env.g;                % Force
    M_ast     = PD.Inertia.m*PD.Env.g*PD.Tether.L;    % Torque
    
    T       = T*T_ast;
    rR_Edge = rR_Edge*L_ast;
    rR      = rR*L_ast;
    vR      = vR*V_ast;
    aR      = aR*A_ast;
    omegaR  = omegaR*sqrt(PD.Env.g/PD.Tether.L);
    gR      = gR*Gamma_ast;
    FA_R    = FA_R*F_ast;
    Tension = Tension*F_ast;
    rQ      = rQ*L_ast;
    R_KE    = R_KE;
    rK      = rK*L_ast; 
    vK      = vK*V_ast; 
    aK      = aK*A_ast; 
    euler   = euler*180/pi;
    omegaK  = omegaK*sqrt(PD.Env.g/PD.Tether.L); 
    gK      = gK*Gamma_ast; 
    FA_K    = FA_K*F_ast; 
    FR_K    = FR_K*F_ast; 
    FG_K    = FG_K*F_ast;
    MA_K    = MA_K*M_ast;
    MR_K    = MR_K*M_ast;  
    MG_K    = MG_K*M_ast;
    MMC_K   = MMC_K*M_ast;
    alfa_K  = alfa_K*180/pi;
    beta_K  = beta_K*180/pi;
    rG      = rG*L_ast;
    vG      = vG*V_ast;
    aG      = aG*A_ast;
    omegaG  = omegaG*sqrt(PD.Env.g/PD.Tether.L); 
    gG      = gG*Gamma_ast;
    FA_G    = FA_G*F_ast;
    FK_G    = FK_G*F_ast; 
    MA_G    = MA_G*M_ast; 
    MK_G    = MK_G*M_ast;
    MMC_G   = MMC_G*M_ast;
    
    xc_out(1:2,1)           = PD.Tether.L*xc_amp(1:2,1);
    xc_out(3:4,1)           = xc_amp(3:4,1)*180/pi;
    xc_out(4+1:4+NG,1)      = M_ast*xc_amp(4+1:4+NG,1);
    xc_out(4+NG+1:4+NG+3,1) = xc_amp(4+NG+1:4+NG+3,1)*180/pi;
    
    if PND.Control.Type == 5  %% Close loop
      xs_target_out(1:2*NR+3,1) = (180/pi)*xs_target_out(1:2*NR+3);
    end
    
end


end