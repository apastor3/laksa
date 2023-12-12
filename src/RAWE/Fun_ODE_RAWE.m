function [DF DG r_N alpha_N r_t alpha_t m_0 kappa tau Ten F_NL RHS_r RHS_alpha] = Fun_ODE_RAWE(tau0,xs,E_m_0_tt)

global PP NP s_fe s_nc
global B_alpha_inv Shape Shape_s Shape_ss Shape_sss
global phi phi_s phi_ss phi_sss varphi varphi_s varphi_ss
global P0_Inv Q0_Inv p0 P4 Q1 P_iN

% Load useful Parameters
N_fe = NP(1);
N_nc = NP(4);
Nu   = 3*(N_fe -1 ) + N_fe + 1;

% Load Physical Parameters
mu    = PP(1);
beta  = PP(2);
sigma = PP(3);
delta = PP(4);
nu    = PP(12);
vw    = PP(13);

% Compute the Boundary conditions
[r0 rN alpha_s_N r0_t rN_t rN_tt m_N_tt] = Compute_BC(tau0,PP);

% Recover the physical vectors at the nodes
[r_N alpha_N r_t alpha_t ] = From_xs_to_Physical(tau0,xs,NP,PP);


% Torque at s = 0 and Torque at s = 1
m_N     = alpha_s_N; 
m_0     = m_N    + xs(2*Nu+1);
m_0_tt  = m_N_tt + E_m_0_tt; 

% Compute vector V  for later computation of alpha
[V_alpha V_alpha_s V_alpha_ss V_alpha_tt] = Fun_Global_Valpha(NP,m_0,m_N,B_alpha_inv,Shape,Shape_s,Shape_ss,m_0_tt,m_N_tt);


% Compute m_s
%mm    = zeros(N_nc,1);
mm_s  = zeros(N_nc,1);
mm_tt = zeros(N_nc,1);
for i=1:1:NP(1)+1
 %   mm(:,1)    = mm(:,1)      +  V_alpha(i,:)'; 
    mm_s(:,1)  = mm_s(:,1)    +  V_alpha_s(i,:)';
    mm_tt(:,1) = mm_tt(:,1)   +  V_alpha_tt(i,:)'; 
end

% Evaluate the function at the NC points
[r_NC r_NC_s r_NC_ss r_NC_sss alpha_NC alpha_NC_s alpha_NC_ss ] = Compute_Functions_NC(NP,r_N,alpha_N,phi,phi_s,phi_ss,phi_sss,varphi,varphi_s,varphi_ss,V_alpha,V_alpha_s,V_alpha_ss); 


% Evaluate the velocity at the NC points
v_NC        = zeros(N_nc,3);  
for i=1:1:NP(1)+1
   v_NC(:,1) = v_NC(:,1) + r_t(i,1)*phi(i,:)';
   v_NC(:,2) = v_NC(:,2) + r_t(i,2)*phi(i,:)';
   v_NC(:,3) = v_NC(:,3) + r_t(i,3)*phi(i,:)';    
end

% Compute auxiliary functions
kappa   = zeros(N_nc,1);
gamma   = zeros(N_nc,1);
Mod_r_s = zeros(N_nc,1);
Ten     = zeros(N_nc,1);
tau     = zeros(N_nc,1);
F_NL    = zeros(N_nc,3);
F_A     = zeros(N_nc,3);

for i=1:1:N_nc
    % Central line Curvature
    kappa(i)   = sqrt(r_NC_ss(i,:)*r_NC_ss(i,:)');
    % Central line Torsion 
    if abs(kappa(i))>0
        gamma(i) = r_NC_s(i,:)*(cross(r_NC_ss(i,:),r_NC_sss(i,:)))'/kappa(i)^2;
    end
    % Modulus of the tension force
    Mod_r_s(i) = sqrt(r_NC_s(i,:)*r_NC_s(i,:)');
    Ten(i)       = PP(3)*(Mod_r_s(i) - 1);
    
    %  Total twist
    tau(i)     = gamma(i) +  alpha_NC_s(i);
    
    % Nonlinear Force
    F_NL(i,:)  = (Ten(i) - mu*kappa(i)^2)*r_NC_s(i,:) + beta*tau(i)*cross(r_NC_s(i,:),r_NC_ss(i,:));
    
    % Aerodynamic force
    va_vec     = v_NC(i,:) - vw*[-1 0 0];
    va         = sqrt(va_vec*va_vec');
    F_A(i,:)   = -nu*va*va_vec;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  RHS for  r-Equations               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RHS_r = zeros(N_fe-1,3);
for i=1:1:N_fe-1
    RHS_r(i,:) = p0(i+1)*[0 0 1] - P_iN(i)*rN_tt(:,1)';
    for j=1:1:3 
        AUX        = F_A(:,j).*phi(i+1,:)' - F_NL(:,j).*phi_s(i+1,:)';
        RHS_r(i,j) = RHS_r(i,j) - mu*P4(i+1,:)*r_N(:,j) +  trapz(s_nc,AUX);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  RHS for  alpha-Equations           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


RHS_alpha = zeros(N_fe+1,1);
for i=1:1:N_fe+1
    RHS_alpha(i,1) = delta*( tau(end)*varphi(i,end) - tau(1)*varphi(i,1) );
    RHS_alpha(i,1) = RHS_alpha(i,1) - delta*Q1(i,:)*alpha_N;
     
   % AUX            = delta*(gamma + mm_s).*varphi_s(i,:)' + mm_tt.*varphi(i,:)';
   % RHS_alpha(i,1) = RHS_alpha(i,1) - trapz(s_nc,AUX);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               RHS                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DF = zeros(2*Nu+1,1);

% dr/dtau = r_dot
DF(0*(N_fe-1)+1:1*(N_fe-1),1) = r_t(2:N_fe,1);
DF(1*(N_fe-1)+1:2*(N_fe-1),1) = r_t(2:N_fe,2);
DF(2*(N_fe-1)+1:3*(N_fe-1),1) = r_t(2:N_fe,3);

% Dalpha/dtau = alpha_t
DF(3*(N_fe-1)+1:Nu,1) = alpha_t;

% d^2 r/dtau^2
DF(Nu + 0*(N_fe-1)+1:Nu+1*(N_fe-1),1) = P0_Inv*RHS_r(:,1);
DF(Nu + 1*(N_fe-1)+1:Nu+2*(N_fe-1),1) = P0_Inv*RHS_r(:,2);
DF(Nu + 2*(N_fe-1)+1:Nu+3*(N_fe-1),1) = P0_Inv*RHS_r(:,3);

% d^2 alpha/dtau^2
DF(Nu + 3*(N_fe-1)+1:2*Nu,1) =  Q0_Inv*RHS_alpha(:,1);

% Control
DF(2*Nu +1)  = PP(10)*( alpha_t(1) - PP(9)) + PP(11)*DF(Nu + 3*(N_fe-1)+1);

DG            = [r_t(2:N_fe,1); r_t(2:N_fe,2);r_t(2:N_fe,3);alpha_t;RHS_r(:,1);RHS_r(:,2);RHS_r(:,3);RHS_alpha(:,1);DF(2*Nu +1)];

end