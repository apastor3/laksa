function PD = Fun_PD_KA

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
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
% Kite Geometry and Inertia  
PD.Kite.c     =  1.5;            % Kite Chord                          (m)
PD.Kite.b     =  5.8;            % Kite span                           (m)
PD.Kite.m     =  4.0;            % Kite mass                           (kg)
PD.Kite.h     =  3.2;            % Kite height                         (m, only plotting purposes)
PD.Kite.hg    =  2.0;            % Z-Distance between G and kite tips  (m, only plotting purposes)
PD.Kite.A     =  14.4;           % Kite surface                        (m^2)  
PD.Kite.Ix    =  21.1;           % Ix                                  (kg m^2)
PD.Kite.Iy    =  4.66;           % Iy                                  (kg m^2)
PD.Kite.Iz    =  18.0;           % Iz                                  (kg m^2) 
PD.Kite.Ixz   =  0.0;            % Ixz                                 (kg m^2)

% Force Aerodynamic coefficients
PD.Aero.Full   =  0;             % Set 1 to use full model     

PD.Aero.Cx0    = -0.065;         % Cx0                          (-)  
PD.Aero.Cxalfa =  0.18;          % Cx_alfa                      (-)  
PD.Aero.Cybeta = -1.57;          % Cy_beta                      (-)  
PD.Aero.Cz0    =  0.116;         % Cz0                          (-)  
PD.Aero.Czalfa = -2.97;          % Cz_alfa                      (-)  

% Torque Aerodynamic coefficients
PD.Aero.Clbeta = -0.1;          % Cl_beta                      (-)     
PD.Aero.Clp    = -0.15;          % Cl_p_tilde                   (-)    
PD.Aero.Cm0    =  0.13;          % Cm0                          (-)    
PD.Aero.Cmalfa = -0.76;          % Cm_alfa                      (-)    
PD.Aero.Cmq    = -0.17;          % Cm_q_tilde                   (-)      
PD.Aero.Cnbeta = -0.027;         % Cn_beta                      (-)    
PD.Aero.Cnr    = -0.002;         % Cn_r_tilde                   (-)       
PD.Aero.Vref   =  7;             % V_ref                        (m/s)     

% Control ->  This control is not implemented in the current version of the code
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
PD.Env.eps    = 0.1;             % Wind Speed fluctuation level
PD.Env.Omega  = 0.6;             % Wind Speed fluctuation frequency  (rad/s)

% Tether characteristics
PD.Tether.L0  = 200;              % L0 main  line reference length       (m)
PD.Tether.XA  = PD.Kite.c/2;     % X-Body coordinate of Attachment point A (m) 
PD.Tether.YA  = PD.Kite.b/2;     % Y-Body coordinate of Attachment point A (m) 
PD.Tether.ZA  = 2;               % Z-Body coordinate of Attachment point A (m) 

% Control parameters 
PD.Ctr.Type   = 1;               % 0 ->  l = sqrt[1-(YA/L0)^2] and delta = 0
                                 % 1 ->  l = sqrt[1-(YA/L0)^2] and delta = delta1*sin(omega*t)                            
                                 % 2 ->  l = sqrt[1-(YA/L0)^2] + l1*sin(omega_l*t) and delta = delta1*sin(omega_delta*t) 
PD.Ctr.l1       = 1.0;           % l1                                (m)
PD.Ctr.Om_l     = 0.31;          % omega_l                           (rad/s)
PD.Ctr.delta1   = 8.0;           % delta1                            (deg)
PD.Ctr.Om_delta = 0.1;           % omega_delta                       (rad/s) 


% Numerical Parameters
PD.Num.RelTol  = 1e-3;           % Integrator Relative Tolerance
PD.Num.AbsTol  = 1e-6;           % Integrator Absolute Tolerance
PD.Num.NewTol  = 5e-12;           % Newton-Raphson Tolerance
PD.Num.MaxIt   = 100;            % Newton-Raphson Maximum number of iterations
PD.Num.DTmax   = 0.001;          % Maximum Time step
PD.Num.dh      = 1e-6;           % Numerical Jacobian step


end