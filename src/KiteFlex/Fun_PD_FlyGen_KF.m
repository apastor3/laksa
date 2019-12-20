function PD = Fun_PD_FlyGen_KF

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Physical Parameters of a Fly-Generation System                 %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %%
% Inputs: No inputs                                                       %%
%                                                                         %%
% Outputs: PD  -> Physical parameters of the system                       %% 
%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kite Geometry and Inertia  
PD.Inertia.c   =  1.0;                        % Drone wing Chord                      (m)
PD.Inertia.b   =  10.0;                       % Drone wing span                       (m)
PD.Inertia.h   =  2.54;                       % Drone height (only plotting purposes) (m)
PD.Inertia.m   =  80.0;                       % Drone mass                            (kg)
PD.Inertia.A   =  PD.Inertia.c*PD.Inertia.b;  % Drone Wing surface   (m^2)  

% Elevator 
PD.AerC.St     = 0.1*PD.Inertia.A;   % Tail Surface                       (m^2)
PD.AerC.lt     = 2.0;                  % Distance from G to the AC of the tail (m)
% Rudder
PD.AerC.Sv     = 0.25*PD.Inertia.A;   % Vertical stabilizer Surface                       (m^2)
PD.AerC.hv     = 1.0;                % Distance from G to the AC of the vertical stabilizer (m)
% Ailerons
PD.AerC.Sa     = 0.1*PD.Inertia.A;   % Aileron  Surface                   (m^2)
PD.AerC.ya     = 3.0;                % Distance from G to the AC of the ailerons (m)

% Equivalent aerea density to estimations the tensor of inertia
ro_S           = PD.Inertia.m/(PD.Inertia.A + PD.AerC.St + PD.AerC.Sv);
Ix             = PD.Inertia.m/12*PD.Inertia.b^2 +...
                 ro_S*PD.AerC.Sv*PD.AerC.hv^2;                                      % Ix (kg*m^2)
Iy             = ro_S*PD.Inertia.A/12*PD.Inertia.c^2 +...
                 ro_S*(PD.AerC.St+PD.AerC.Sv)*PD.AerC.lt^2;                         % Iy (kg*m^2)
Iz             = ro_S*PD.Inertia.A/12*(PD.Inertia.c^2 + PD.Inertia.b^2)+...
                 ro_S*(PD.AerC.St+PD.AerC.Sv)*PD.AerC.lt^2;                         % Iz (kg*m^2) 
Ixz            = ro_S*PD.AerC.Sv*PD.AerC.lt*PD.AerC.hv;                             % Ixz      (kg*m^2) 

PD.Inertia.IK  = [Ix       0           Ixz;...
                  0        Iy          0;...
                  Ixz      0           Iz];
% Environment parameters
PD.Env.g         =  9.81;         % Earth acceleration                 (m/s^2)
PD.Env.ro        =  1.225;        % Air density                        (kg/m^3)

PD.Env.Type   =  0;              % 0 -> Wind Velocity is constant
                                 % 1 -> Wind Speed = Vw*(h/H0)^alfa*(1+eps*sin(Omega*t))
PD.Env.Vw     = 15;              % Wind Velocity                     (m/s)
PD.Env.alfa   = 0.14;            % Exponent of the wind speed law
PD.Env.H0     = 10;              % Height Scale                      (m)
PD.Env.eps    = 0.1;             % Wind Speed fluctuation level
PD.Env.Omega  = 0.6;             % Wind Speed fluctuation frequency  (rad/s)

% Tether characteristics
PD.Tether.L      = 300;           % L0 main  line reference length     (m)
PD.Tether.ro     = 970;           % Tether density (kg/m^3)
PD.Tether.dt     = 6e-3;          % Tether diameter (m)
PD.Tether.Cd     = 1.0;           % Tether Cd  

% Bridle geometry
PD.Bridle.lb    = 3.0;            % l0 bridle line reference length    (m)
PD.Bridle.delta = 70*pi/180;      % delta                              (rad)
PD.Bridle.eta   = 0.;             % eta                                (rad)

% Force Aerodynamic coefficients
PD.Aero.Full   =  0;                % Set 1 to use full model     

PD.Aero.Cx0    = -0.065;         % Cx0                          (-)  % Groot
PD.Aero.Cxalfa =  0.18;          % Cx_alfa                      (-)  % Groot
PD.Aero.Cybeta = -1.57;          % Cy_beta                      (-)  % Williams
PD.Aero.Cz0    =  0.12;          % Cz0                          (-)  % Groot
PD.Aero.Czalfa = -2.97;          % Cz_alfa                      (-)  % Groot

% Torque Aerodynamic coefficients
PD.Aero.Clbeta = -0.49;          % Cl_beta                      (-)       % Williams
PD.Aero.Clp    = -0.15;          % Cl_p_tilde                   (-)       % Williams Clp = Clp * b/2Va
PD.Aero.Cm0    =  0.13;          % Cm0                          (-)       % Groot
PD.Aero.Cmalfa = -0.76;          % Cm_alfa                      (-)       % Groot
PD.Aero.Cmq    = -0.17;          % Cm_q_tilde                   (-)       % Williams     Cmq = Cmq*c/Va  
PD.Aero.Cnbeta = -0.027;         % Cn_beta                      (-)       % Williams
PD.Aero.Cnr    = -0.002;         % Cn_r_tilde                   (-)       % Williams Cn_r = Cnr*b/2Va
PD.Aero.Vref   =  7;             % V_ref                        (m/s)     % Groot

% Aerodynamic Model Limits (only for postprocess checking purposes)
PD.Aero.alfa_s =  25;            % Stall angle                  (º)
PD.Aero.beta_m =  15;            % Maximum sideslip angle       (º)


% Aerodynamic Control Surfaces
% Elevator
PD.AerC.at     = 5.5;                % Corrected Lift slope (1/rad)  
% Rudder
PD.AerC.av     = 5.5;                % Corrected Lift slope (1/rad)  
% Ailerons
PD.AerC.aa     = 5.5;                % Corrected Lift slope (1/rad)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Characteristics of the generators       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FlyGen/GroundGen
PD.Gen.Num      =  4; % Number of on-board Generators

% Inertia of the generators
PD.Gen.m        =  0.3; % kg, Mass of the three blade
PD.Gen.L        =  0.4; % m, length of a blade

% Position of the generators with respect to the center of mass
PD.Gen.x(1)     = PD.Inertia.c/2;
PD.Gen.x(2)     = PD.Inertia.c/2;
PD.Gen.x(3)     = PD.Inertia.c/2;
PD.Gen.x(4)     = PD.Inertia.c/2;


PD.Gen.y(1)     =  PD.Inertia.b/4;
PD.Gen.y(2)     = -PD.Inertia.b/4;
PD.Gen.y(3)     =  PD.Inertia.b/4;
PD.Gen.y(4)     = -PD.Inertia.b/4;

PD.Gen.z(1)     =  0.5;
PD.Gen.z(2)     =  0.5;
PD.Gen.z(3)     = -0.5;
PD.Gen.z(4)     = -0.5;

% Angle of the generators (rad) 
PD.Gen.nu(1)     =  0;
PD.Gen.nu(2)     =  0;
PD.Gen.nu(3)     =  0;
PD.Gen.nu(4)     =  0;

% Aerodynamic coefficients
PD.Gen.Cf        =  0.08;
PD.Gen.Cm        =  0.1;

% Control parameters
PD.Control.Type      = 3;         % 0 -> No control
                                  % 1 -> Periodic laws
                                  % 2 -> Target variables
                                  % 3 -> Equilibrium Calculations
                                 
PD.Control.T0        = 20;        % Control period 

PD.Control.eps_lt    = 2.5;       % lt = lt0*(1+eps_lt*sin(ome_lt*t))
PD.Control.ome_lt    = 0.5;       % 

PD.Control.eps_lb    = 1.25;      % lb = lb0*(1+eps_lb*sin(ome_lb*t))
PD.Control.ome_lb    = 0.2;       % 

PD.Control.eps_del   = 0.2;       % delta = delta0*(1+eps_del*sin(ome_del*t))
PD.Control.ome_del   = 0.15;      % 

PD.Control.eps_eta   = 0.12;      % eta   = eta0*(1+eps_eta*sin(ome_eta*t))
PD.Control.ome_eta   = 0.2;       %

PD.Control.TL        = 1.5;       % Characteristic Time in the steering manoeuvre 
PD.Control.TC        = 1.2;       % Characteristic Time in the steering manoeuvre 

% Close-Loop Control Parameters
PD.CL.Target.Type    =  10;                    % Type of Target trajectory   
PD.CL.Target.Control =  1;                     % Control mean 0-> Aerodynamic surfaces, 1-> bridle control

% Target state and control variables
PD.Target.Ome_Rotor = 3500;       % Target angular velocity of the rotors (rpm)
PD.Target.xi_Rotor  = 0;          % Target motor controller torque (N m)
PD.Target.delta_a   = 0.8;        % Target deflection of the aileron (deg)
PD.Target.delta_r   = 0.5;        % Target deflection of the rudder (deg)
PD.Target.delta_e   = 1.0;        % Target deflection of the elevator (deg)


%% Numerical Parameters
PD.Num.N       = 3;              % Number of bars
PD.Num.IntTol  = 1e-6;           % Integrator Tolerance
PD.Num.NewTol  = 5e-11;          % Newton-Raphson Tolerance
PD.Num.DT      = 0.001;          % Time step
PD.Num.NDT     = 100;            % Number of time steps
PD.Num.MaxIt   = 100;            % Newton-Raphson Maximum number of iterations
PD.Num.dh      = 1e-8;           % Numerical Jacobian step

end