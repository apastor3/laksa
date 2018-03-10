function [rk vk ak omega_BE alfa_BE T_Bp T_Bm m_Bp m_Bm fa ma alfa beta ...
          Rp Rm Elong_p Elong_m DF RHS] = Fun_ODE_Full_Output_KS(t,X,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute all the relevant variables                             %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %%
% Inputs:  t               -> dimensionless time                       %%
%          X = [xs xs_dot] -> extended state vector                    %%
%          PND             -> Dimensionless parameters                 %%
%          FE              -> Finite element data                      %% 
% Outputs: rk -> Kite position vector (SE components)                  %%
%          vk -> Kite velocity vector (SE components)                  %%
%          ak -> Kite acceleration vector (SE components)              %%
%          omega_BE -> Kite angular velocity (SE components)           %%
%          alfa_BE  -> Kite angular acceleration (SE components)       %%
%          T_Bp     -> Tension upon the kite at B plus (SB components) %%
%          T_Bm     -> Tension upon the kite at B minus (SB components)%%
%          m_Bp     -> Torque about G due to the tension upon the      %%
%                      kite at B plus (SB components)                  %%
%          m_Bm     -> Torque about G due to the tension upon the      %% 
%                      kite at B minus (SB components)                 %%
%          fa       -> aerodynamic force (SB components)               %%
%          ma       -> aerodynamic torque about G (SB components)      %%
%          alfa     -> attack angle                                    %%
%          beta     -> sideslip angle                                  %%
%          Rp       -> Position vector of the elastic tether           %%
%                      attached to B plus                              %%
%          Rm        -> Position vector of the elastic tether          %%
%                       attached to B Minus                            %%
%          Elong_p   -> Elastic tether elongation                      %% 
%          Elong_m   -> Elastic tether elongation                      %%
%          RHS       -> Right-Hand side of the equations               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Data_Exp    % Emperimental Data
%%%%%%%%%%%%%%%%%%%%% Preliminary calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Recover the kite state vector
xk         = X(1:5,1);
xk_p       = X(6:10,1);

% Compute rotation matrices
R1E    = Fun_R1E_KS(xk);
R21    = Fun_R21_KS(xk);
RB2    = Fun_RB2_KS(xk);
R32    = Fun_R32_KS(xk);

RBE    = RB2*R21*R1E;
R2E    = R21*R1E;
R3E    = R32*R2E;

% Absolute Angular velocity in the Earth frame
[Phi Phi_xs ]   = Fun_Matrix_Omega_KS(xk);  
omega_body      = Phi*xk_p;
omega_BE        = RBE'*omega_body;                            % Omega of the Body frame with respect to the Earth frame

% Center of mass velocity in the Earth frame
[Ups_s Ups_s_xs] = Fun_Matrix_Upsilon_KS(xk,PND);
vk_body          = Ups_s*xk_p;
vk               = RBE'*vk_body;

% Matrix M
[Ms Ms_s  ]      = Fun_Matrix_M_KS(Phi,Phi_xs,Ups_s,Ups_s_xs,PND);

% Componentns of the kite position vector in Earth frame
OEF_OE          = -R3E'*[0 0 PND.Bar.Ls]';
OE_O2           = -R2E'*[0 0 PND.Tether.l]';  
O2_G            = -RBE'*[PND.Tether.XA 0 PND.Tether.ZA]';
rk              =  OEF_OE+OE_O2+O2_G;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Preliminary calculations about the kite to compute the       %%
%%%%      position and velocity of points B plus and B minus         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control vector 
[xc xc_p xc_pp] = Fun_Control_KS(t,X,PND,rk,vk,'');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PR              = xc(1,1);                     % Power ratio 
nu              = xc(2,1);                     % Bar deflection angle
PR_p            = xc_p(1,1);
nu_p            = xc_p(2,1);

d               = PND.Bar.Lds  + PR*(PND.Bar.Ls  - PND.Bar.Lds - PND.Bar.Lps );
d_p             =              PR_p*(PND.Bar.Ls  - PND.Bar.Lds - PND.Bar.Lps );

BarD_Kin_2      = [d d_p ]';

% Position vectors of the bar's points C0, C Plus and C Minus 
OE_C0    =  R3E'*[0,            0,                    d]'; 
C0_CP    =  R3E'*[0,  0.5*PND.Bar.Lc*cos(nu),  0.5*PND.Bar.Lc*sin(nu)]';
C0_CM    =  R3E'*[0, -0.5*PND.Bar.Lc*cos(nu), -0.5*PND.Bar.Lc*sin(nu)]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Compute the positions of the tips of the elastic tethers        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S_E components of the position vector of Points C_Plus and C_Minus 
Rp(1,:)       = OEF_OE + OE_C0 + C0_CP; % C_Plus
Rm(1,:)       = OEF_OE + OE_C0 + C0_CM; % C_Minus
% S_E Components  of the position vector of points B_Plus and B_Minus 
Rp(2,:)       = (rk + RBE'*[PND.Tether.XB,  PND.Tether.YB, PND.Tether.ZB]')'; % B Plus
Rm(2,:)       = (rk + RBE'*[PND.Tether.XB, -PND.Tether.YB, PND.Tether.ZB]')'; % B Minus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    tension and moment upon the kite due to the elastic tethers      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Tension_p TMod_p Elong_p Elong_t_p] = Tension_KS(xk,xk_p,xc,xc_p,+1,R2E,BarD_Kin_2,PND); % Tension acting on the kite at  B+
[Tension_m TMod_m Elong_m Elong_t_m] = Tension_KS(xk,xk_p,xc,xc_p,-1,R2E,BarD_Kin_2,PND); % Tension acting on the kite at  B-    
up0           = -(Rp(2,:)-Rp(1,:))';
up0           = up0/sqrt(up0'*up0);

um0           = -(Rm(2,:)-Rm(1,:))';
um0           = um0/sqrt(um0'*um0);

T_Bp          =  RBE*Tension_p;           % Tension at point B+ (component in the body frame) 
T_Bm          =  RBE*Tension_m;           % Tension at point B- (component in the body frame) 

% Compute the torque about G of the tensions of the flexible tethers
O_B_B_p =  [PND.Tether.XB  PND.Tether.YB  PND.Tether.ZB]';
O_B_B_m =  [PND.Tether.XB -PND.Tether.YB  PND.Tether.ZB]';
m_Bp    =  cross(O_B_B_p,T_Bp);     % Components in the Body frame
m_Bm    =  cross(O_B_B_m,T_Bm);     % Components in the Body frame 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     Compute Generalized Forces                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wind speed components in the Earth Frame
vw              = Fun_Wind(t,rk,PND);
% Vectors in the body frame
vw_body         = RBE*vw;

% Generalized forces
[Q f m alfa beta fa ma ] = Fun_Q_KS(vw_body,vk_body,omega_body,Ups_s,Phi,T_Bp,T_Bm,m_Bp,m_Bm,PND);     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     Compute Potential contribution                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U U_xs] = Fun_Potential_KS(xk,PND);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     Compute Right hand Side                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary operations
Aux1 = zeros(5,5); % (partial Ms/partial xs)*xs
Aux5 = zeros(5,1); %
for k=1:1:5
    Aux1(:,:) = Aux1 + Ms_s(:,:,k)*xk_p(k);
    Aux5(k,1) =   xk_p'*Ms_s(:,:,k)*xk_p;
end
RHS = Q-U_xs + 0.5*Aux5-Aux1*xk_p;

% Right-hand side of the Equations
DF(1:5,1)  = xk_p; 
DF(6:10,1) = Ms\RHS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Compute the RHS of the equations of motion               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xk_pp         = DF(6:10,1);
% Acceleration of the kite and relevant angular accelerations
[ak alfa_BE ] = Fun_Acce_KS(RBE,vk,xk_p,xk_pp,omega_BE,Phi,Phi_xs,Ups_s,Ups_s_xs);

if PND.Ctr.Type == 3
   [xc xc_p xc_pp] = Fun_Control_KS(t,X,PND,rk,vk,ak);
   DF(11,1)        = xc_p(2,1);
end

end