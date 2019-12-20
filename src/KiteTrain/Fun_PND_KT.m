function PND = Fun_PND_KT(PD)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and  Jose A. Serrano-Iglesia           %
% Language  : Matlab                                                         %
% Synopsis  : Compute Dimensionless parameters                               %
% Copyright : Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: structure PD with the physical parameters                %%
%                                                                  %%
% Outputs: structure PND with the dimensionless parameters         %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PND.Kite.N       = PD.Kite.N;                % Number of Kites        

% Characteristic variables
L0               =  PD.Tether.L(1);
M0               =  PD.Kite.m(1);
for i=1:1:PND.Kite.N
   % Kite Geometry and Inertia
    PND.Kite.mu(i)    =  PD.Env.ro*PD.Kite.A(i)*L0/(2*M0);                   % mu
    PND.Kite.c(i)     =  PD.Kite.c(i)/L0;                                    % eps_c
    PND.Kite.b(i)     =  PD.Kite.b(i)/L0;                                    % eps_b
    PND.Kite.h(i)     =  PD.Kite.h(i)/L0;                                    % eps_h (only plotting purposes)
    PND.Kite.hg(i)    =  PD.Kite.hg(i)/L0;                                   % eps_hg (only plotting purposes)
    PND.Kite.sigma(i) =  PD.Kite.m(i)/M0;                                    % Mass ratio
   
    iota_Ix        = PD.Kite.Ix(i)/(M0*L0^2);               % Ix/m*L0^2   
    iota_Iy        = PD.Kite.Iy(i)/(M0*L0^2);               % Iy/m*L0^2   
    iota_Iz        = PD.Kite.Iz(i)/(M0*L0^2);               % Iz/m*L0^2    
    iota_Ixz       = PD.Kite.Ixz(i)/(M0*L0^2);              % Ixz/m*L0^2 

    PND.Kite.iota(:,:,i) = [iota_Ix 0 iota_Ixz; 0 iota_Iy 0;iota_Ixz 0  iota_Iz];
    
    % Tether length
    PND.Tether.L(i)      = PD.Tether.L(i)/L0;
    
    % Attachment points 
    PND.Tether.XA(i)     = PD.Tether.XA(i)/L0;                           % XA/L0 
    PND.Tether.YA(i)     = PD.Tether.YA(i)/L0;                           % YA/L0 
    PND.Tether.ZA(i)     = PD.Tether.ZA(i)/L0;                           % ZA/L0 

    PND.Tether.XC(i)     = PD.Tether.XC(i)/L0;                           % XC/L0 
    PND.Tether.YC(i)     = PD.Tether.YC(i)/L0;                           % YC/L0 
    PND.Tether.ZC(i)     = PD.Tether.ZC(i)/L0;                           % ZC/L0 
    
end
% Aerodynamic Force and Torque Coefficients
PND.Aero.Full =  PD.Aero.Full;

if PD.Aero.Full ==  1
    PND.Aero.CX = PD.Aero.CX;
    PND.Aero.CY = PD.Aero.CY;
    PND.Aero.CZ = PD.Aero.CZ;
    
    PND.Aero.Cm = PD.Aero.Cm;
    PND.Aero.Cl = PD.Aero.Cl;
    PND.Aero.Cn = PD.Aero.Cn;
else    
    
    PND.Aero.CX = [0  PD.Aero.Cxalfa PD.Aero.Cx0...    % CX0 
                   0         0              0 ...         % CX_beta
                   0         0              0 ...         % CX_p
                   0         0              0 ...         % CX_q
                   0         0              0 ...         % CX_r 
                   0         0              0 ...         % CX_delta_aileron
                   0         0              0 ...         % CX_delta_elevator
                   0         0              0];           % CX_delta_rudder     

    PND.Aero.CY =  [0         0              0 ...         % CY0 
                   0         0        PD.Aero.Cybeta ...  % CY_beta
                   0         0              0  ...        % CY_p
                   0         0              0  ...        % CY_q
                   0         0              0  ...        % CY_r 
                   0         0              0  ...        % CY_delta_aileron
                   0         0       PD.Aero.Cydelta_r... % CY_delta_elevator
                   0         0              0  ];         % CY_delta_rudder     
             
    PND.Aero.CZ = [ 0    PD.Aero.Czalfa  PD.Aero.Cz0 ...   % CZ0 
                   0         0              0 ...         % CZ_beta
                   0         0              0 ...         % CZ_p
                   0         0              0 ...         % CZ_q
                   0         0              0 ...         % CZ_r 
                   0         0              0 ...         % CZ_delta_aileron
                   0         0              0 ...         % CZ_delta_elevator
                   0         0              0];           % CZ_delta_rudder     
              

    PND.Aero.Cm = [ 0      PD.Aero.Cmalfa  PD.Aero.Cm0 ... % Cm0 
                   0         0              0 ...         % Cm_beta
                   0         0              0 ...         % Cm_p
                   0         0            PD.Aero.Cmq...  % Cm_q
                   0         0              0 ...         % Cm_r 
                   0         0              0 ...         % Cm_delta_aileron
                   0         0       PD.Aero.Cmdelta_e ...% Cm_delta_elevator
                   0         0              0];           % Cm_delta_rudder     

     PND.Aero.Cl= [ 0         0              0 ...         % Cl0 
                   0         0         PD.Aero.Clbeta ... % Cl_beta
                   0         0            PD.Aero.Clp ... % Cl_p
                   0         0              0  ...        % Cl_q
                   0         0              0 ...         % Cl_r 
                   0         0      PD.Aero.Cldelta_a ... % Cl_delta_aileron
                   0         0              0 ...         % Cl_delta_elevator
                   0         0      PD.Aero.Cldelta_r];   % Cl_delta_rudder   
               
     PND.Aero.Cn= [ 0         0             0 ...         % Cn0 
                   0         0         PD.Aero.Cnbeta ... % Cn_beta
                   0         0              0         ... % Cn_p
                   0         0              0  ...        % Cn_q
                   0         0           PD.Aero.Cnr...   % Cn_r 
                   0         0              0 ... % Cn_delta_aileron
                   0         0              0 ...         % Cn_delta_elevator
                   0         0      PD.Aero.Cndelta_r ];  % Cn_delta_rudder         
end 

PND.Aero.vt      = PD.Aero.Vref/sqrt(PD.Env.g*L0);            % V_ref/sqrt(g*L0)     

% Aerodynamic Model Limits (only for postprocess checking purposes)
PND.Aero.alfa_s =  PD.Aero.alfa_s*pi/180;                               % Stall angle(rad)
PND.Aero.beta_m =  PD.Aero.beta_m*pi/180;                               % Maximum sideslip angle (rad)

% Enviromental Properties
PND.Env.Type     = PD.Env.Type;                                         % Wind law
PND.Env.vw       = PD.Env.Vw/sqrt(PD.Env.g*L0);                         % V_0 tilde 
PND.Env.alfa     = PD.Env.alfa;                                         % Exponent of the wind speed law
PND.Env.H0       = PD.Env.H0/L0;                                        % Height Scale
PND.Env.eps      = PD.Env.eps;                                          % Wind Speed fluctuation level
PND.Env.Omega    = PD.Env.Omega*sqrt(L0/PD.Env.g);                      % Normalized fluctuation frequency 

% Control 
PND.Ctr.Type     = PD.Ctr.Type;                 % 0 ->   Control Surfaces deflections are constant

PND.Ctr.delta_a  = PD.Ctr.delta_a*pi/180;       % Aileron deflection (rad)
PND.Ctr.delta_r  = PD.Ctr.delta_r*pi/180;       % Rudder deflection (rad)
PND.Ctr.delta_e  = PD.Ctr.delta_e*pi/180;       % Elevator deflection (rad)

% Numerical Parameters
PND.Num.RelTol   = PD.Num.RelTol;        % Integrator Relative Tolerance
PND.Num.AbsTol   = PD.Num.AbsTol;        % Integrator Absolute Tolerance
PND.Num.NewTol   = PD.Num.NewTol;        % Newton-Raphson Tolerance
PND.Num.MaxIt    = PD.Num.MaxIt;         % Newton-Raphson Maximum number of iterations
PND.Num.DTmax    = PD.Num.DTmax;         % Maximum Time step
PND.Num.dh       = PD.Num.dh;            % Numerical Jacobian step


end