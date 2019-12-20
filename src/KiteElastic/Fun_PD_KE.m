function PD = Fun_PD_Paper_10_Kite_KE

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Physical parameters                                            %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs: No inputs                                                %%
%                                                                    %%
%   Outputs: structure PD with the physical parameters               %% 
%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kite Geometry and Inertia  
PD.Kite.Num   =  10;                % Total Number of Kites
for i=1:1:PD.Kite.Num 
    PD.Kite.c(i)     =  1.5;              % Kite Chord                          (m)
    PD.Kite.b(i)     =  5.8;              % Kite span                           (m)
    PD.Kite.m(i)     =  4.0;              % Kite mass                           (kg)
    PD.Kite.h(i)     =  3.0;              % Kite height                         (m)
    PD.Kite.hg(i)    =  2.0;              % Z-Distance between G and kite tips  (m)
    PD.Kite.A(i)     =  14.4;             % Kite surface                        (m^2)  
    PD.Kite.Ix(i)    =  21.1;             % Ix                                  (kg m^2)
    PD.Kite.Iy(i)    =  4.7;              % Iy                                  (kg m^2)
    PD.Kite.Iz(i)    =  17.9;             % Iz                                  (kg m^2) 
    PD.Kite.Ixz(i)   =  0.;               % Ixz                                 (kg m^2)

    % Force Aerodynamic coefficients
    PD.Aero.Full(i)      =  0;                % Set 1 to use full model     

    PD.Aero.Cx0(i)    =  -0.065;               % Cx0                          (-)  
    PD.Aero.Cxalfa(i) =   0.18;              % Cx_alfa                      (-)  
    PD.Aero.Cybeta(i) =  -1.6;               % Cy_beta                      (-) 
    PD.Aero.Cz0(i)    =   0.12;             % Cz0                          (-)  
    PD.Aero.Czalfa(i) =  -3.0;              % Cz_alfa                      (-)  

    % Torque Aerodynamic coefficients
    PD.Aero.Clbeta(i) =  0.1;                             % Cl_beta                      (-)
    PD.Aero.Clp(i)    = -0.15;                             % Cl_p_tilde                   (-)      
    PD.Aero.Cm0(i)    =  0.13;                            % Cm0                          (-)      
    PD.Aero.Cmalfa(i) = -0.76;                            % Cm_alfa                      (-)      
    PD.Aero.Cmq(i)    = -0.17;                            % Cm_q_tilde                   (-)       
    PD.Aero.Cnbeta(i) = -0.03;                             % Cn_beta                      (-)      
    PD.Aero.Cnr(i)    = -0.002;                             % Cn_r_tilde                   (-)         

    % Control
    PD.Aero.Cydelta_r(i) = 0;
    PD.Aero.Cmdelta_e(i) = 0;
    PD.Aero.Cldelta_a(i) = 0;
    PD.Aero.Cldelta_r(i) = 0;
    PD.Aero.Cndelta_r(i) = 0;

    PD.Aero.Vref(i)   =  7;             % V_ref                        (m/s)    

    % Aerodynamic Model Limits (only for postprocess checking purposes)
    PD.Aero.alfa_s(i) =  25;            % Stall angle                  (º)
    PD.Aero.beta_m(i) =  15;            % Maximum sideslip angle       (º)
end

% Environment parameters
PD.Env.g      =  9.81;           % Earth acceleration                (m/s^2)
PD.Env.ro     =  1.225;           % Air density                       (kg/m^3)
PD.Env.Type   =  0;              % 0 -> Wind Velocity is constant
                                 % 1 -> Wind Speed = Vw*(h/H0)^alfa*(1+eps*sin(Omega*t))
PD.Env.Vw     = 7;               % Wind Velocity                     (m/s)
PD.Env.alfa   = 0.14;            % Exponent of the wind speed law
PD.Env.H0     = 10;              % Height Scale                      (m)
PD.Env.eps    = 0.0;             % Wind Speed fluctuation level
PD.Env.Omega  = 0.0;             % Wind Speed fluctuation frequency  (rad/s)

% Common Tether characteristics
PD.Tether.L0  = 100.0;                                          % Reference length   (m)    
PD.Tether.E   = 90e9; %10e9;                                    % Tether Young's Modulus                            (Pa)
PD.Tether.dt  = 2.0e-3;                                         % Tether diameter                                   (m)
PD.Tether.ft  = 0.01*sqrt(PD.Tether.L0/PD.Env.g);               % Dimensionless friction coefficient  
PD.Tether.ro  = 970;                                            % Tether density (kg/m^3)
PD.Tether.Cd  = 1.0;                                            % Tether Cd  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cretate list of tether   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PD.Tether.Num      =  2*PD.Kite.Num;             % TotalNumber of tethers
% Tethers linking the ground and kite number 1 
PD.Tether.L00(1)   =  PD.Tether.L0;       % Natural length of the tether 
PD.Tether.Mass(1)  =  3;                  % Number of masses used to discretize the tether
PD.Tether.Down(1)  =  0;                  % Number of the Kite attached to the Lowest tether point (0 means the Earth System)
PD.Tether.Dx(1)    =  0;                  % X coordinate of the down-point (SK coordinates)
PD.Tether.Dy(1)    =  0;                  % Y coordinate of the down-point (SK coordinates)
PD.Tether.Dz(1)    =  0;                  % Z coordinate of the down-point (SK coordinates)
PD.Tether.Up(1)    =  1;                  % Number of the Kite attached to the Upper tether point
PD.Tether.Ux(1)    =  0.75;             % X coordinate of the Up-point (SK coordinates)
PD.Tether.Uy(1)    =  2.9;               % Y coordinate of the Up-point (SK coordinates)
PD.Tether.Uz(1)    =  2;                  % Z coordinate of the Up-point (SK coordinates)

PD.Tether.L00(2)   =  PD.Tether.L0;       % Natural length of the tether 
PD.Tether.Mass(2)  =  PD.Tether.Mass(1);  % Number of masses used to discretize the tether
PD.Tether.Down(2)  =  0;                  % Number of the Kite attached to the Lowest tether point (0 means the Earth System)
PD.Tether.Dx(2)    =  0;                  % X coordinate of the down-point (SK coordinates)
PD.Tether.Dy(2)    =  0;                  % Y coordinate of the down-point (SK coordinates)
PD.Tether.Dz(2)    =  0;                  % Z coordinate of the down-point (SK coordinates)
PD.Tether.Up(2)    =  1;                  % Number of the Kite attached to the Upper tether point
PD.Tether.Ux(2)    =  PD.Tether.Ux(1);    % X coordinate of the Up-point (SK coordinates)
PD.Tether.Uy(2)    = -PD.Tether.Uy(1);    % Y coordinate of the Up-point (SK coordinates)
PD.Tether.Uz(2)    =  PD.Tether.Uz(1);    % Z coordinate of the Up-point (SK coordinates)

for i=2:1:PD.Kite.Num
    PD.Tether.L00(2*(i-1)+1)   =  PD.Tether.L0;       % Natural length of the tether 
    PD.Tether.Mass(2*(i-1)+1)  =  PD.Tether.Mass(1);  % Number of masses used to discretize the tether
    PD.Tether.Down(2*(i-1)+1)  =  i-1;                % Number of the Kite attached to the Lowest tether point (0 means the Earth System)
    PD.Tether.Dx(2*(i-1)+1)    =  0;                  % X coordinate of the down-point (SK coordinates)
    PD.Tether.Dy(2*(i-1)+1)    =  0;                  % Y coordinate of the down-point (SK coordinates)
    PD.Tether.Dz(2*(i-1)+1)    =  0;                  % Z coordinate of the down-point (SK coordinates)
    PD.Tether.Up(2*(i-1)+1)    =  i;                  % Number of the Kite attached to the Upper tether point
    PD.Tether.Ux(2*(i-1)+1)    =  PD.Tether.Ux(1);    % X coordinate of the Up-point (SK coordinates)
    PD.Tether.Uy(2*(i-1)+1)    =  PD.Tether.Uy(1);    % Y coordinate of the Up-point (SK coordinates)
    PD.Tether.Uz(2*(i-1)+1)    =  PD.Tether.Uz(1);    % Z coordinate of the Up-point (SK coordinates)

    PD.Tether.L00(2*i)   =  PD.Tether.L0;       % Natural length of the tether 
    PD.Tether.Mass(2*i)  =  PD.Tether.Mass(1);  % Number of masses used to discretize the tether
    PD.Tether.Down(2*i)  =  i-1;                  % Number of the Kite attached to the Lowest tether point (0 means the Earth System)
    PD.Tether.Dx(2*i)    =  0;                  % X coordinate of the down-point (SK coordinates)
    PD.Tether.Dy(2*i)    =  0;                  % Y coordinate of the down-point (SK coordinates)
    PD.Tether.Dz(2*i)    =  0;                  % Z coordinate of the down-point (SK coordinates)
    PD.Tether.Up(2*i)    =  i;                  % Number of the Kite attached to the Upper tether point
    PD.Tether.Ux(2*i)    =  PD.Tether.Ux(1);    % X coordinate of the Up-point (SK coordinates)
    PD.Tether.Uy(2*i)    = -PD.Tether.Uy(1);    % Y coordinate of the Up-point (SK coordinates)
    PD.Tether.Uz(2*i)    =  PD.Tether.Uz(1);    % Z coordinate of the Up-point (SK coordinates)   
end

% Control 
PD.Ctr.Type    = 0;               % 0 ->   Control Surfaces deflections are constant

for i=1:1:PD.Kite.Num    % Total Number of aircraft
                                        % delta = delta0 + delta1*cos(Om_delta*t)
    PD.Ctr.delta_a0(i) = 0.0;           % Aileron deflection (deg)
    PD.Ctr.delta_r0(i) = 0.0;           % Rudder deflection (deg)
    PD.Ctr.delta_e0(i) = 0.0;           % Elevator deflection (deg)
    
    PD.Ctr.delta_a1(i) = 0.0;           % Aileron deflection (deg)
    PD.Ctr.delta_r1(i) = 0.0;           % Rudder deflection (deg)
    PD.Ctr.delta_e1(i) = 0.0;           % Elevator deflection (deg)
    
    PD.Ctr.Om_delta_a(i) = 0.0;         % Angular frequency of Aileron deflection (rad/s)
    PD.Ctr.Om_delta_r(i) = 0.0;         % Angular frequency of Rudder deflection (rad/s)
    PD.Ctr.Om_delta_e(i) = 0.0;         % Angular frequency of Elevator deflection (rad/s) 
end

for i=1:1:PD.Tether.Num  % Total Number of tethers
                                        % L(t) = PD.Tether.L00*( 1 + PD.Ctr.eps_L * cos(Om_L*t))
    PD.Ctr.eps_L(i)      = 0.0;         % Length modulation (dimensionless)
    PD.Ctr.Om_L(i)       = 0.0;         % angular frequency of length modulation (rad/s) 
end


% Numerical Parameters
PD.Num.RelTol  = 1e-6;           % Integrator Relative Tolerance
PD.Num.AbsTol  = 1e-10;           % Integrator Absolute Tolerance
PD.Num.NewTol  = 5e-8;           % Newton-Raphson Tolerance
PD.Num.MaxIt   = 100;            % Newton-Raphson Maximum number of iterations
PD.Num.dh      = 1e-6;           % Numerical Jacobian step




end