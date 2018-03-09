function PND = Fun_PND_KA(PD)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute Dimensionless parameters                               %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: structure PD with the physical parameters                %%
%                                                                  %%
% Outputs: structure PND with the dimensionless parameters         %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kite Geometry and Inertia
PND.Kite.mu      =  PD.Env.ro*PD.Kite.A*PD.Tether.L0/(2*PD.Kite.m );    % mu
PND.Kite.c       =  PD.Kite.c/PD.Tether.L0;                             % eps_c
PND.Kite.b       =  PD.Kite.b/PD.Tether.L0;                             % eps_b
PND.Kite.h       =  PD.Kite.h/PD.Tether.L0;                             % eps_h (only plotting purposes)
PND.Kite.hg      =  PD.Kite.hg/PD.Tether.L0;                            % eps_hg (only plotting purposes)


iota_Ix          = PD.Kite.Ix/(PD.Kite.m*PD.Tether.L0^2);               % Ix/m*L0^2   
iota_Iy          = PD.Kite.Iy/(PD.Kite.m*PD.Tether.L0^2);               % Iy/m*L0^2   
iota_Iz          = PD.Kite.Iz/(PD.Kite.m*PD.Tether.L0^2);               % Iz/m*L0^2    
iota_Ixz         = PD.Kite.Ixz/(PD.Kite.m*PD.Tether.L0^2);              % Ixz/m*L0^2 

PND.Kite.iota    = [iota_Ix 0 iota_Ixz; 0 iota_Iy 0;iota_Ixz 0  iota_Iz];

% Force Aerodynamic coefficients
PND.Aero.Cx0     = PD.Aero.Cx0;                                         % Cx0  
PND.Aero.Cxalfa  = PD.Aero.Cxalfa;                                      % Cx_alfa                 
PND.Aero.Cybeta  = PD.Aero.Cybeta;                                      % Cy_beta                    
PND.Aero.Cz0     = PD.Aero.Cz0;                                         % Cz0                        
PND.Aero.Czalfa  = PD.Aero.Czalfa;                                      % Cz_alfa                    

% Torque Aerodynamic coefficients
PND.Aero.Clbeta  = PD.Aero.Clbeta;                                      % Cl_beta
PND.Aero.Clp     = PD.Aero.Clp;                                         % Cl_p_tilde        
PND.Aero.Cm0     = PD.Aero.Cm0;                                         % Cm0                         
PND.Aero.Cmalfa  = PD.Aero.Cmalfa;                                      % Cm_alfa                
PND.Aero.Cmq     = PD.Aero.Cmq;                                         % Cm_q_tilde                  
PND.Aero.Cnbeta  = PD.Aero.Cnbeta;                                      % Cn_beta                 
PND.Aero.Cnr     = PD.Aero.Cnr;                                         % Cn_r_tilde                 
PND.Aero.vt      = PD.Aero.Vref/sqrt(PD.Env.g*PD.Tether.L0);            % V_ref/sqrt(g*L0) 

% Aerodynamic Model Limits (only for postprocess checking purposes)
PND.Aero.alfa_s =  PD.Aero.alfa_s*pi/180;                               % Stall angle(rad)
PND.Aero.beta_m =  PD.Aero.beta_m*pi/180;                               % Maximum sideslip angle (rad)

% Enviromental Properties
PND.Env.Type     = PD.Env.Type;                                         % Wind law
PND.Env.vw       = PD.Env.Vw/sqrt(PD.Env.g*PD.Tether.L0);               % V_0 tilde 
PND.Env.alfa     = PD.Env.alfa;                                         % Exponent of the wind speed law
PND.Env.H0       = PD.Env.H0/PD.Tether.L0;                              % Height Scale
PND.Env.eps      = PD.Env.eps;                                          % Wind Speed fluctuation level
PND.Env.Omega    = PD.Env.Omega*sqrt(PD.Tether.L0/PD.Env.g);            % Normalized fluctuation frequency 

% Tether Properties
PND.Tether.XA    = PD.Tether.XA/PD.Tether.L0;                           % XA/L0 
PND.Tether.YA    = PD.Tether.YA/PD.Tether.L0;                           % YA/L0 
PND.Tether.ZA    = PD.Tether.ZA/PD.Tether.L0;                           % ZA/L0 

% Control parameters 
PND.Ctr.Type     = PD.Ctr.Type;                                         % 0 ->  l = sqrt(1-yA^2) and delta = 0
                                                                        % 1 ->  l = sqrt(1-yA^2) and delta = delta1*sin(omega*t)                                 
                                                                        % 2 ->  l = sqrt(1-yA^2) + l1*sin(omega_l*t) and delta = delta1*sin(omega_delta*t) 
PND.Ctr.l1       = PD.Ctr.l1/PD.Tether.L0;                              % l1/L0
PND.Ctr.Om_l     = PD.Ctr.Om_l*sqrt(PD.Tether.L0/PD.Env.g);             % Normalized omega_l1
PND.Ctr.delta1   = PD.Ctr.delta1*pi/180;                                % delta1 (rad)         
PND.Ctr.Om_delta = PD.Ctr.Om_delta*sqrt(PD.Tether.L0/PD.Env.g);         % Normalized omega_delta 
                
% Numerical Parameters
PND.Num.RelTol   = PD.Num.RelTol;        % Integrator Relative Tolerance
PND.Num.AbsTol   = PD.Num.AbsTol;        % Integrator Absolute Tolerance
PND.Num.NewTol   = PD.Num.NewTol;        % Newton-Raphson Tolerance
PND.Num.MaxIt    = PD.Num.MaxIt;         % Newton-Raphson Maximum number of iterations
PND.Num.DTmax    = PD.Num.DTmax;         % Maximum Time step
PND.Num.dh       = PD.Num.dh;            % Numerical Jacobian step


end