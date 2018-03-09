function PD = Fun_PD_Close_Loop_KF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This subroutine belongs to KiteFlex   (KF)/LAKSA                       %%
%% Author: Gonzalo S?nchez-Arriaga and Alejandro Pastor-Rodr?guez         %%
%% August 20th, 2017                                                      %%
%% Universidad Carlos III de Madrid, Spain                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% See detailed description of the code in the  paper                     %%
%%                                                                        %%
%%                                                                        %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %%
% Inputs: No inputs                                                       %%
%                                                                         %%
% Outputs: PD  -> Physical parameters of the system                       %% 
%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kite Geometry and Inertia  
PD.Inertia.c    =  0.25;              % Drone wing Chord                      (m)
PD.Inertia.b    =  3.0;               % Drone wing span                       (m)
PD.Inertia.h    =  0.05;              % Drone height (only plotting purposes) (m)
PD.Inertia.m    =  2.0;               %0.9195;               % Drone mass                            (kg)
PD.Inertia.A    =  0.75;              % Drone Wing surface   (m^2)  

% Elevator 
PD.AerC.St      = 0.1*PD.Inertia.A;   % Tail Surface                       (m^2)
PD.AerC.lt      = 0.7;                % Distance from G to the AC of the tail (m)
% Rudder
PD.AerC.Sv      = 0.036*PD.Inertia.A; % Vertical stabilizer Surface                       (m^2)
PD.AerC.hv      = 0.05;               % Distance from G to the AC of the vertical stabilizer (m)
% Ailerons
PD.AerC.Sa      = 0.03*PD.Inertia.A; % Aileron  Surface                   (m^2)
PD.AerC.ya      = 1.0;                % Distance from G to the AC of the ailerons (m)

% Tensor of Inertia
Ix              =  0.20;             % Ix (kg*m^2)
Iy              =  0.78e-1;          % Iy (kg*m^2)
Iz              =  0.28;             % Iz (kg*m^2) 
Ixz             = -0.2e-2;           % Ixz      (kg*m^2) 

PD.Inertia.IK   = [Ix       0          Ixz;...
                  0        Iy          0;...
                  Ixz      0           Iz];
% Environment parameters
PD.Env.g        =  9.81;              % Earth acceleration                 (m/s^2)
PD.Env.ro       =  1.225;             % Air density                        (kg/m^3)

PD.Env.Type     =  0;                 % 0 -> Wind Velocity is constant
                                      % 1 -> Wind Speed = Vw*(h/H0)^alfa*(1+eps*sin(Omega*t))
PD.Env.Vw       =  7;                 % Wind Velocity                     (m/s)
PD.Env.alfa     =  0.14;              % Exponent of the wind speed law
PD.Env.H0       =  10;                % Height Scale                      (m)
PD.Env.eps      =  0.1;               % Wind Speed fluctuation level
PD.Env.Omega    =  0.6;               % Wind Speed fluctuation frequency  (rad/s)

% Tether characteristics
PD.Tether.L     = 30;                 % L0 main  line reference length     (m)
PD.Tether.ro    = 970;                % Tether density (kg/m^3)
PD.Tether.dt    = 2e-3;               % Tether diameter (m)
PD.Tether.Cd    = 1.0;                % Tether Cd  

% Bridle geometry
PD.Bridle.lb    = 3.0;                % l0 bridle line reference length    (m)
PD.Bridle.delta = 80*pi/180;          % delta                              (rad)
PD.Bridle.eta   = 0.;                 % eta                                (rad)

% Force Aerodynamic coefficients
PD.Aero.Vref   =  5;                  % V_ref                        (m/s)     % BD
PD.Aero.Cx0    = -0.025;              % Cx0                          (-)  % BD
PD.Aero.Cxalfa =  0.67;               % Cx_alfa                      (-) 
PD.Aero.Cybeta = -0.40;               % Cy_beta                      (-)  % BD
PD.Aero.Cz0    = -0.91;               % Cz0                          (-)  % BD
PD.Aero.Czalfa = -5.65;               % Cz_alfa                      (-)  % BD

% Torque Aerodynamic coefficients
PD.Aero.Clbeta = -0.26;               % Cl_beta                      (-)       % BD
PD.Aero.Clp    = -0.2;                % Cl_p_tilde                   (-)       % BD
PD.Aero.Cm0    =  0.074;              % Cm0                          (-)       
PD.Aero.Cmalfa = -0.86;               % Cm_alfa                      (-)       
PD.Aero.Cmq    = -0.3;                % Cm_q_tilde                   (-)           
PD.Aero.Cnbeta =  0.052;              % Cn_beta                      (-)       
PD.Aero.Cnr    = -0.02;              % Cn_r_tilde                   (-)       

% Aerodynamic Control Surfaces
% Elevator
PD.AerC.at     = 5.5;                % Corrected Lift slope (1/rad)  
% % Rudder
PD.AerC.av     = 5.5;                % Corrected Lift slope (1/rad)  
% % Ailerons
PD.AerC.aa     = 5.5;                % Corrected Lift slope (1/rad)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Characteristics of the generators       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FlyGen/GroundGen
PD.Gen.Num      =  2; % Number of on-board Generators

% Inertia of the generators
PD.Gen.m        =  0.3; % kg, Mass of the three blade
PD.Gen.L        =  0.2; % m, length of a blade

% Position of the generators with respect to the center of mass
PD.Gen.x(1)     = PD.Inertia.c/2;
PD.Gen.x(2)     = PD.Inertia.c/2;

PD.Gen.y(1)     =  PD.Inertia.b/4;
PD.Gen.y(2)     = -PD.Inertia.b/4;

PD.Gen.z(1)     =  0.0;
PD.Gen.z(2)     =  0.0;


% Angle of the generators (rad) 
PD.Gen.nu(1)     =  0;
PD.Gen.nu(2)     =  0;


% Aerodynamic coefficients
PD.Gen.Cf        =  0.08;
PD.Gen.Cm        =  0.1;

% Control parameters
PD.Control.Type      = 5;          % 0 -> No control
                                   % 1 -> Periodic laws
                                   % 2 -> Target variables
                                   % 3 -> Target deflection of the aileron
                                   % 4 -> Reel-In
                                   % 5 -> Close-Loop
                                  
PD.Control.T0        = 100;        % Control period 

PD.Control.eps_lt    = 1;          % lt = lt0*(1+eps_lt*sin(ome_lt*t))
PD.Control.ome_lt    = - 1e-6 ;        

PD.Control.eps_lb    = 0;          % lb = lb0*(1+eps_lb*sin(ome_lb*t))
PD.Control.ome_lb    = 0.2;        % 

PD.Control.eps_del   = 0;          % delta = delta0*(1+eps_del*sin(ome_del*t))
PD.Control.ome_del   = 0.15;       % 

PD.Control.eps_eta   = 0;          % eta   = eta0*(1+eps_eta*sin(ome_eta*t))
PD.Control.ome_eta   = 0.2;        %

PD.Control.eps_lin   = 0.01;       % lt = lt0*(1 - eps_lt*t)

PD.Control.TL        = 1.5;        % Characteristic Time in the steering manoeuvre 
PD.Control.TC        = 1.2;        % Characteristic Time in the steering manoeuvre 
         
% Close-Loop Control Parameters
PD.CL.Target.Type    =  0;                     % Type of Target trajectory   
PD.CL.Target.Control =  0;                     % Control mean 0-> Aerodynamic surfaces, 1-> bridle control

% Target state and control variables 
PD.Target.Ome_Rotor = 3500;        % Target angular velocity of the rotors (rpm)
PD.Target.xi_Rotor  = 0;           % Target motor controller torque (N m)
PD.Target.delta_a   = 0.0;         % Target deflection of the aileron (deg)
PD.Target.delta_r   = 0.0;         % Target deflection of the rudder (deg)
PD.Target.delta_e   = 0.0;         % Target deflection of the elevator (deg)


%% Numerical Parameters
PD.Num.N       = 3;              % Number of bars
PD.Num.IntTol  = 1e-7;           % Integrator Tolerance
PD.Num.NewTol  = 5e-11;          % Newton-Raphson Tolerance
PD.Num.DT      = 0.001;          % Time step
PD.Num.NDT     = 100;            % Number of time steps
PD.Num.MaxIt   = 100;            % Newton-Raphson Maximum number of iterations
PD.Num.dh      = 1e-8;           % Numerical Jacobian step

end