function PD = Fun_PD_Paper_Modes_KT

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %           
% Language  : Matlab                                                         %
% Synopsis  : Physical Parameters                                            %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: No inputs                                                %%
%                                                                  %%
% Outputs: structure PD with the physical parameters               %% 
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PD.Kite.N     =  10;              % Number of Kites        

% Kite Geometry and Inertia 
for i=1:1:PD.Kite.N
    PD.Kite.c(i)     =  1.5;               % Kite Chord                          (m)
    PD.Kite.b(i)     =  5.8;               % Kite span                           (m)
    PD.Kite.m(i)     =  4.0;               % Kite mass                           (kg)
    PD.Kite.h(i)     =  3.0;               % Kite height                         (m, only plotting purposes)
    PD.Kite.hg(i)    =  2.0;               % Z-Distance between G and kite tips  (m, only plotting purposes)
    PD.Kite.A(i)     =  14.4;             % Kite surface                        (m^2)  
    PD.Kite.Ix(i)    =  21.1;             % Ix                                  (kg m^2)
    PD.Kite.Iy(i)    =  4.7;              % Iy                                  (kg m^2)
    PD.Kite.Iz(i)    =  17.9;              % Iz                                  (kg m^2) 
    PD.Kite.Ixz(i)   =  0.;                % Ixz                                 (kg m^2)
    
    % Tether Length 
    PD.Tether.L(i)   = 100;                 % L0 main  line reference length       (m)
    
    % Attachment points
    PD.Tether.XA(i)  = 0.75;                % X-Body coordinate of Attachment point A (m) 
    PD.Tether.YA(i)  = 2.9;                 % Y-Body coordinate of Attachment point A (m) 
    PD.Tether.ZA(i)  = 2.;                  % Z-Body coordinate of Attachment point A (m)  
    
    PD.Tether.XC(i)  = 0.;                  % X-Body coordinate of Attachment point C (m) 
    PD.Tether.YC(i)  = 0.0;                 % Y-Body coordinate of Attachment point C (m) 
    PD.Tether.ZC(i)  = 0.0;                   % Z-Body coordinate of Attachment point C (m) 

end

% Force Aerodynamic coefficients
PD.Aero.Full      =  0;                % Set 1 to use full model     

PD.Aero.Cx0    = -0.065;         % Cx0                          (-)  
PD.Aero.Cxalfa =  0.18;          % Cx_alfa                      (-)  
PD.Aero.Cybeta = -1.6;           % Cy_beta                      (-)  
PD.Aero.Cz0    =  0.12;          % Cz0                          (-)  
PD.Aero.Czalfa =  -3.;          % Cz_alfa                      (-)  

% Torque Aerodynamic coefficients
PD.Aero.Clbeta =  0.1;          % Cl_beta                      (-)     
PD.Aero.Clp    = -0.15;          % Cl_p_tilde                   (-)    
PD.Aero.Cm0    =  0.13;          % Cm0                          (-)    
PD.Aero.Cmalfa = -0.76;          % Cm_alfa                      (-)    
PD.Aero.Cmq    = -0.17;          % Cm_q_tilde                   (-)      
PD.Aero.Cnbeta = -0.03;         % Cn_beta                      (-)    
PD.Aero.Cnr    = -0.002;         % Cn_r_tilde                   (-)       
PD.Aero.Vref   =  7;             % V_ref                        (m/s)     

% Control
PD.Aero.Cydelta_r = 0;
PD.Aero.Cmdelta_e = 0;
PD.Aero.Cldelta_a = 0;
PD.Aero.Cldelta_r = 0;
PD.Aero.Cndelta_r = 0;

% Aerodynamic Model Limits (only for postprocess checking purposes)
PD.Aero.alfa_s =  25;            % Stall angle                  (º)
PD.Aero.beta_m =  15;            % Maximum sideslip angle       (º)

% Environment parameters
PD.Env.g      =  9.81;           % Earth acceleration                (m/s^2)
PD.Env.ro     =  1.225;          % Air density                       (kg/m^3)
PD.Env.Type   =  0;              % 0 -> Wind Velocity is constant
                                 % 1 -> Wind Speed = Vw*(h/H0)^alfa*(1+eps*sin(Omega*t))
PD.Env.Vw     = 7;               % Wind Velocity                     (m/s)
PD.Env.alfa   = 0.14;            % Exponent of the wind speed law
PD.Env.H0     = 10;              % Height Scale                      (m)
PD.Env.eps    = 0.0;             % Wind Speed fluctuation level
PD.Env.Omega  = 0.0;             % Wind Speed fluctuation frequency  (rad/s)

% Control 
PD.Ctr.Type    = 0;               % 0 ->   Control Surfaces deflections are constant

PD.Ctr.delta_a = 0.0;           % Aileron deflection (deg)
PD.Ctr.delta_r = 0.0;           % Rudder deflection (deg)
PD.Ctr.delta_e = 0.0;           % Elevator deflection (deg)


% Numerical Parameters
PD.Num.RelTol  = 1e-3;           % Integrator Relative Tolerance
PD.Num.AbsTol  = 1e-6;           % Integrator Absolute Tolerance
PD.Num.NewTol  = 5e-12;          % Newton-Raphson Tolerance
PD.Num.MaxIt   = 100;            % Newton-Raphson Maximum number of iterations
PD.Num.DTmax   = 0.001;          % Maximum Time step
PD.Num.dh      = 1e-5;           % Numerical Jacobian step


end