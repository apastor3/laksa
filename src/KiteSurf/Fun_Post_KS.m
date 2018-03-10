function [T RBE R2E R3E rk vk ak euler_k omega_k alfa_k Lambda FAP FAM MAP MAM FBP FBM MBP MBM ...
          FA MA W alfa beta Rp Rm Elong_p Elong_m xc Error0] = Fun_Post_KS(PD,t,X,Flag_Dim,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Postprocess the results                                        %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: t        -> Time                                             %%
%         X        -> Extended State Vector                            %%
%         Flag_Dim -> 0: outputs are dimensionless                     %%  
%                     1: outputs have dimensions)                      %%
%         PND      -> Dimensionless parameters                         %%
%                                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs: RBE      -> SB-SE Rotation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          T         -> Time                                           %%
%          RBE       -> SB-SE rotation matrix                          %%
%          R2E       -> S2-SE rotation matrix                          %%
%          R3E       -> S3-SE rotation matrix                          %%
%          rk        -> Earth components of  kite position vector      %%
%          vk        -> Earth components of  kite velocity vector      %%
%          ak        -> Earth components of  kite acceleration vector  %% 
%          Euler     -> Euler angles                                   %%
%          omega_k   -> Earth components of  kite angular velocity     %% 
%          alfa_k    -> Earth components of  kite angular acceleration %%
%                                                                      %%
%          Lambda    -> Modulus of the Inelastic tether tension        %%
%                       >0 implies traction                            %%                   
%          FAP       -> Force at A Plus (SE components)                %%
%          FAM       -> Force at A minus (SE components)               %%
%          MAP       -> Moment due to Force at A Plus (SE components)  %%
%          MAM       -> Moment due to Force at A Minus(SE components)  %%
%          FBP       -> Force at B Plus (SE components)                %%
%          FBM       -> Force at B minus (SE components)               %%
%          MBP       -> Moment due to Force at B Plus (SE components)  %%
%          MBM       -> Moment due to Force at B Minus(SE components)  %%
%          FA        -> Earth components of the Aerodynamic Force      %% 
%          MA        -> Earth components of the Aerodynamic Torque     %%
%                       about G                                        %% 
%          W         -> Earth components of the weight of the kite     %%
%          alfa      -> Attack angle                                   %%
%          beta      -> Sideslip angle                                 %% 
%                                                                      %%  
%          RP        -> Earth components of the position vector of     %%
%                         the nodes of the tether attached to B Plus   %% 
%          Tp        -> Earth components of the tension at B_Plus      %%
%          Rm        -> Earth components of the position vector of     %%
%                         the nodes of the tether attached to B Minus  %%
%          Tm        -> Earth components of the tension  at B_Minus    %%
%          Elong_p   -> Elastic tether elongation                      %% 
%          Elong_m   -> Elastic tether elongation                      %%
%          xc        -> Control vector                                 %%
%          Error0    -> Checking from Classical Mechanics              %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the computations by using dimensionless variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate the RHS while getting all the outputs
[rk vk ak omega_k alfa_k T_Bp T_Bm m_Bp m_Bm fa ma alfa beta ...
          Rp Rm Elong_p Elong_m DF RHS] = Fun_ODE_Full_Output_KS(t,X,PND);
                
% Recover the kite state vector and control
xk              = X(1:5,1);
[xc xc_p xc_pp] = Fun_Control_KS(t,X,PND,rk,vk,ak);
% Rotation matrices
R1E      = Fun_R1E_KS(xk);
R21      = Fun_R21_KS(xk);
RB2      = Fun_RB2_KS(xk);
R32      = Fun_R32_KS(xk);

RBE      = RB2*R21*R1E;
R2E      = R21*R1E;  
R3E      = R32*R2E;
%% Euler Angles
euler_k(2,1)  = -asin(RBE(1,3));                                                % Theta (cabeceo)
euler_k(1,1)  =  atan2(RBE(1,2)/cos(euler_k(2,1)),RBE(1,1)/cos(euler_k(2,1)));  % Psi (guiñada)
euler_k(3,1)  =  atan2(RBE(2,3)/cos(euler_k(2,1)),RBE(3,3)/cos(euler_k(2,1)));  % Phi (balance)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recover control and state variables%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
l0       = PND.Tether.l;
Ls       = PND.Bar.Ls;
%% Importante position vectors (components in the Earth system)
OE_O2    = -R2E'*[0 0  l0]';
O2_OB    = -RBE'*[PND.Tether.XA        0         PND.Tether.ZA]';

OB_Ap    =  RBE'*[PND.Tether.XA   PND.Tether.YA  PND.Tether.ZA]';
OB_Am    =  RBE'*[PND.Tether.XA  -PND.Tether.YA  PND.Tether.ZA]';

OE_Ap    =  OE_O2+O2_OB+OB_Ap;
OE_Am    =  OE_O2+O2_OB+OB_Am;

u_Ap     =  -OE_Ap/sqrt(OE_Ap'*OE_Ap);
u_Am     =  -OE_Am/sqrt(OE_Am'*OE_Am);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Forces and Moments acting on the kite %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weight of the Kite
W           =  [0 0 1]';    
% Tension at point B+
FBP         = RBE'*T_Bp;  % Tension at B+ (components in the Earth frame)
MBP         = RBE'*m_Bp;  % Torque (component in the Earth frame)
% Tension at point B-
FBM         = RBE'*T_Bm;   % Tension at B- (components in the Earth frame)
MBM         = RBE'*m_Bm;   % Torque (components in the Earth frame)
% Aerodynamic forces
FA          = RBE'*fa; % Earth Components of the aerodynamic force 
MA          = RBE'*ma; % Earth Components of the aerodynamic torque 
% Compute Tensions of the rigid tether 
M0          = [1  u_Ap'*u_Am; u_Ap'*u_Am 1 ];
b0(1,1)     = u_Ap'*(ak-FA-FBP-FBM-W); 
b0(2,1)     = u_Am'*(ak-FA-FBP-FBM-W); 
Aux         = inv(M0)*b0;
Lambda(1,1) = Aux(1);
Lambda(2,1) = Aux(2);

FAP         = Aux(1,1)*u_Ap;    % Force upon the kite due to the rigid tether (components in the Earth frame) 
FAM         = Aux(2,1)*u_Am;    % Force upon the kite due to the rigid tether (components in the Earth frame)
MAP         = cross(OB_Ap,FAP); % Torque (component in the Earth frame)
MAM         = cross(OB_Am,FAM); % Torque (component in the Earth frame)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check Calculations using classical mechanics  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L          = PND.Kite.iota*(RBE*omega_k);    % Angular momentum components in SB
Error(1,1) = max(abs(ak-(W+FA+FAP+FAM+FBP+FBM))); 
Error(2,1) = max(abs(PND.Kite.iota*RBE*alfa_k+cross(RBE*omega_k,L)-RBE*(MA+MAP+MAM+MBP+MBM)));

Error0     = max(abs(Error));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%           Outputs          %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T    =  t;             % Time

if Flag_Dim == 1
  
    L0       = PD.Tether.Ll;        % Reference length (m)
    g        = PD.Env.g;            % Graviational acceleration (m/s^2)
    M        = PD.Kite.m;           % Kite mass                  (kg)
    tast     = sqrt(PD.Tether.Ll/PD.Env.g);
    
    T        = t*tast;
    
    rk       = L0*rk;
    vk       = sqrt(g*L0)*vk;
    ak       = g*ak;
    euler_k  = 180/pi*euler_k;
    omega_k  = sqrt(g/L0)*omega_k; 
    alfa_k   = (g/L0)*alfa_k;
    
    Lambda   = Lambda*M*g;
    
    FAP      =  FAP*M*g; 
    FAM      =  FAM*M*g; 
    MAP      =  MAP*M*g*L0; 
    MAM      =  MAM*M*g*L0; 
    
    FBP      =  FBP*M*g; 
    FBM      =  FBM*M*g; 
    MBP      =  MBP*M*g*L0; 
    MBM      =  MBM*M*g*L0; 
    
    FA       =  FA*M*g; 
    MA       =  MA*M*g*L0; 
    W        =  W*M*g; 
    alfa     =  alfa*180/pi;
    beta     =  beta*180/pi;
    Rp       =  Rp*L0;
    Rm       =  Rm*L0;    
    
    xc(2,1)  = xc(2,1)*180/pi;
end


 
 
end
