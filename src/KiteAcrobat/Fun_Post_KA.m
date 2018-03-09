function [T RBE rk vk ak euler omega omega_p Lambda FAP FAM MAP MAM FA MA W alfa beta LP LM ] = Fun_Post_KA(PD,PND,t,u,Flag_Dim)
        

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Result Postprocessing                                          %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Inputs: PD        -> Physical parameters of the system            %%
%         PND       -> Dimensionless parameters                     %%
%         t         -> Dimensionless time                           %% 
%         u         -> Extended state vector                        %% 
%         Flag_Dim  -> set 0 for dimensionless outputs              %%
%                     set 1 for outputs with dimensions             %%
%                                                                   %%
% Outputs: t        -> Time                                         %%
%          Kinematics quantities of the kite                        %%
%                      RBE   -> Body-Earth rotation matrix          %%
%                      rk    -> Earth components of the position,   %%
%                      vk       velocity and acceleration vectors   %%
%                      ak       of the center of mass of the kite   %% 
%                      euler -> Euler angles                        %%
%                      omega -> Earth components of the kite        %%
%                      omega_p  angular velocity and acceleration   %%
%          Kite forces and moments                                  %%
%                      Lambda -> Tether Tensions along              %%
%                                tether directions                  %%
%                      FAP   -> Earth components of the Forces      %%
%                      FAM      exerted by the tethers linked at    %%
%                               points A plus and A Minus           %%
%                      MAP  ->  Earth components of the Torques     %%
%                      MAM      exerted by the tethers linked at    %%
%                               points A plus and A minus           %%
%                      FA   ->  Earth components of the Aerodynamic %% 
%                      MA       force and torque                    %% 
%                      W    -> Earth components of the kite weight  %%
%           Others                                                  %%
%                      alfa -> Angle of attack                      %%  
%                      beta -> Sideslip angle                       %%
%                      LP   -> Length of the A plus-tether          %%
%                      LM   -> Length of the A minus-tether          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recover the state vector and its derivative
xs     = u(1:4,1);
xs_p   = u(5:8,1);

% Recover the control vector
[xc xc_p xc_pp] = Fun_Control_KA(t,PND);

%% Time
T = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Kinematics calculations                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute rotation Matrices
R1E = Fun_R1E_KA(xs);
R21 = Fun_R21_KA(xs);
RB2 = Fun_RB2_KA(xs);
RBE = RB2*R21*R1E;
% Compute the kite position vector components in Earth frame
OEO2            =  R1E'*[0,  xc(1,1)*sin(xs(3)+xc(2,1)), -xc(1,1)*cos(xs(3)+xc(2,1))]'; 
O2G             = -RBE'*[PND.Tether.XA 0 PND.Tether.ZA]';
rk              =  OEO2+O2G;

% Euler angles
euler(2,1)  = -asin(RBE(1,3));                 % Theta (pithc)
euler(1,1)  = asin(RBE(1,2)/cos(euler(2,1)));  % Psi   (yaw)
euler(3,1)  = asin(RBE(2,3)/cos(euler(2,1)));  % Phi   (roll)

% Compute kite velocity vector componentns in the Earth frame
[Ups_s Ups_c Ups_s_xs Ups_s_xc Ups_c_xs Ups_c_xc] = Fun_Matrix_Upsilon_KA(xs,xc,PND);
vk_body = Ups_s*xs_p + Ups_c*xc_p;      % Components in the body frame 
vk      = RBE'*vk_body;                 % Components in the Earth frame

% Compute kite angular velocity vector components in the Earth frame
[Phi Phi_xs]    = Fun_Matrix_Omega_KA(xs);
omega_body      = Phi*xs_p;                  % Kite angular velocity components in Body axes
omega           = RBE'*omega_body;           % Kite angular velocity components in Earth axes

% Compute kite acceleration vector components in the Earth Frame
Aux1            = zeros(3,4); % (partial Upsilon_s/partial xs)*dot{xs}
Aux2            = zeros(3,2); % (partial Upsilon_c/partial xs)*dot{xs}
for k=1:1:length(xs(:,1))
    Aux1(:,:) = Aux1(:,:) + squeeze(Ups_s_xs(:,:,k))*xs_p(k);
    Aux2(:,:) = Aux2(:,:) + squeeze(Ups_c_xs(:,:,k))*xs_p(k);
end
Aux3 = zeros(3,4); % (partial Upsilon_s/partial xc)*dot{xc}
Aux4 = zeros(3,2); % (partial Upsilon_c/partial xc)*dot{xc}
for k=1:1:length(xc(:,1))
    Aux3(:,:) = Aux3(:,:) + squeeze(Ups_s_xc(:,:,k))*xc_p(k);
    Aux4(:,:) = Aux4(:,:) + squeeze(Ups_c_xc(:,:,k))*xc_p(k);
end
DF          = Fun_ODE_Lag_KA(t,u);
ak          = RBE'*(Ups_s*DF(5:8)+(Aux1+Aux3)*xs_p+Ups_c*xc_pp+(Aux2+Aux4)*xc_p+cross(omega_body,vk_body));

% Compute kite angular acceleration vector components in the Earth frame
Aux0 = zeros(3,4); % (Partial Phi/partial xs)*xsp
for k=1:1:length(xs)
    Aux0 = Aux0 + Phi_xs(:,:,k)*xs_p(k);
end
omega_p  = RBE'*(Phi*DF(5:8) + Aux0*xs_p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Forces and moments about G                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Aerodynamic forces components in the Earth frame
vw                  = Fun_Wind(t,rk,PND);
[Q FA MA alfa beta] = Fun_Q_KA(RBE*vw,vk_body,omega_body,Ups_s,Phi,PND);
FA                  = RBE'*FA;
MA                  = RBE'*MA;

% Weight components in the Earth frame
W                   = [0 0 1]'; 
% Compute unit vectors for the tensions of the tethers
GAp        = RBE'*[PND.Tether.XA  PND.Tether.YA PND.Tether.ZA]';              % Attachment point position vector  
GAm        = RBE'*[PND.Tether.XA -PND.Tether.YA PND.Tether.ZA]';              % Attachment point position vector
OEAp       = rk+GAp;                                                          % OE to A plus vector components in Earth frame
OEAm       = rk+GAm;                                                          % OE to A minus vector components in Earh frame
uAp        = -OEAp/sqrt(OEAp'*OEAp);                                          % Unit Vector Components in Earth frame  
uAm        = -OEAm/sqrt(OEAm'*OEAm);                                          % Unit Vector Components in Earth frame 
% Project Newton's Second Law
Mat    = [1 uAp'*uAm;uAp'*uAm 1 ];
b(1,1) = uAp'*(ak-FA-W );
b(2,1) = uAm'*(ak-FA-W );
% Tension along the tether directions
Sol         = inv(Mat)*b;
Lambda(1,1) = Sol(1);
Lambda(2,1) = Sol(2);
% Compute tension components in the Earth frame
FAP    = Sol(1)*uAp;
FAM    = Sol(2)*uAm;
% Compute torques about the center of mass
MAP    = cross(GAp,FAP);
MAM    = cross(GAm,FAM);

% Tether lengths
LP  = sqrt(xc(1,1)^2+PND.Tether.YA^2+2*xc(1,1)*PND.Tether.YA*sin(xc(2,1)));
LM  = sqrt(xc(1,1)^2+PND.Tether.YA^2-2*xc(1,1)*PND.Tether.YA*sin(xc(2,1)));


if Flag_Dim==1  % Put the outputs with dimensions
   
    % Characteristic Values
    L_ast    = PD.Tether.L0;                     % Distance
    V_ast    = sqrt(PD.Tether.L0*PD.Env.g);      % Velocity
    A_ast    = PD.Env.g;                         % Acceleration
    T_ast    = sqrt(PD.Tether.L0/PD.Env.g);      % Time 
    F_ast    = PD.Kite.m*PD.Env.g;               % Force
    M_ast    = PD.Kite.m*PD.Env.g*PD.Tether.L0;  % Torque
    
    T                = T_ast*t;
    rk       = L_ast*rk; 
    vk       = V_ast*vk;
    ak       = A_ast*ak;
    euler    = (180/pi)*euler;
    omega    = omega/T_ast;
    omega_p  = omega_p/T_ast^2;
    
    FAP      =  F_ast*FAP;
    FAM      =  F_ast*FAM;
    MAP      =  M_ast*MAP;
    MAM      =  M_ast*MAM;
    FA       =  F_ast*FA;
    MA       =  M_ast*MA;
    W        =  F_ast*W;
 
    alfa     =  alfa*180/pi;
    beta     =  beta*180/pi;
    LP       =  L_ast*LP;
    LM       =  L_ast*LM;
end





end