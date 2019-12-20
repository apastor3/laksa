function PND = Fun_PND_KF(PD)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Dimensionless Parameters                                       %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %%
% Inputs: PD: physical parameters                                         %%
%                                                                         %%
% Outputs: PND: Dimensionless parameters                                  %% 
%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kite/Drone Geometry and Inertia  

PND.Kite.mu      =  PD.Env.ro*PD.Inertia.A*PD.Tether.L/(2*PD.Inertia.m);  % mu
PND.Kite.c       =  PD.Inertia.c/PD.Tether.L;                             % eps_c
PND.Kite.b       =  PD.Inertia.b/PD.Tether.L;                             % eps_b
PND.Kite.h       =  PD.Inertia.h/PD.Tether.L;                             % h/L0 (plotting purposes only) 
PND.Kite.ik      =  PD.Inertia.IK/(PD.Inertia.m*PD.Tether.L^2);           % Kite tensor of inertia

% Environmental Variables
PND.Env.Type     = PD.Env.Type;                                       % Wind law
PND.Env.vw       = PD.Env.Vw/sqrt(PD.Env.g*PD.Tether.L);              % V_0 tilde 
PND.Env.alfa     = PD.Env.alfa;                                       % Exponent of the wind speed law
PND.Env.H0       = PD.Env.H0/PD.Tether.L;                             % Height Scale
PND.Env.eps      = PD.Env.eps;                                        % Wind Speed fluctuation level
PND.Env.Omega    = PD.Env.Omega*sqrt(PD.Tether.L/PD.Env.g);           % Normalized fluctuation frequency 

% Bridle Geomtry
PND.Bridle.lb    = PD.Bridle.lb/PD.Tether.L;                          % eps_l0
PND.Bridle.delta = PD.Bridle.delta;                                   % delta0
PND.Bridle.eta   = PD.Bridle.eta;                                     % nu0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rotor     Variables  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common variables
PND.Gen.Num   = PD.Gen.Num;                 % Number of on-board Generators
if PND.Gen.Num >0
    PND.Gen.Sigma = PD.Gen.m/PD.Inertia.m;      % Generator-to-Kite mass ratio
    PND.Gen.l     = PD.Gen.L/PD.Tether.L;       % Generator arm-to-initial tether length ratio  
    PND.Gen.iota  = [1/3 0 0;0 1/6 0; 0 0 1/6]; % Generator tensor of inertia (divided by Sigma*lG^2)
    PND.Gen.Cf    = PD.Gen.Cf;
    PND.Gen.Cm    = PD.Gen.Cm;

    PND.Gen_Chi   = PD.Env.ro*acos(-1)*PD.Gen.L^2*PD.Tether.L/(2*PD.Inertia.m); % Dimensionless Parameter for the aerodynamic force 

    for i=1:1:PND.Gen.Num
      PND.Gen.x(i) =  PD.Gen.x(i)/PD.Tether.L;
      PND.Gen.y(i) =  PD.Gen.y(i)/PD.Tether.L;
      PND.Gen.z(i) =  PD.Gen.z(i)/PD.Tether.L;

      PND.Gen.nu(i) =  PD.Gen.nu(i);
    end
end

% Aerodynamic Control Surfaces
%Elevator
Vt                 = PD.AerC.St*PD.AerC.lt/(PD.Inertia.A*PD.Inertia.c);
PD.Aero.Cmdelta_e = -Vt*PD.AerC.at;
% Rudder
Vvz                =  PD.AerC.Sv*PD.AerC.hv/(PD.Inertia.A*PD.Inertia.b);
Vvx                =  PD.AerC.Sv*PD.AerC.lt/(PD.Inertia.A*PD.Inertia.b);

PD.Aero.Cydelta_r =  PD.AerC.Sv/PD.Inertia.A*PD.AerC.av;
PD.Aero.Cldelta_r =  Vvz*PD.AerC.av;
PD.Aero.Cndelta_r = -Vvx*PD.AerC.av;
% Ailerons
Va                 =  PD.AerC.Sa*PD.AerC.ya/(PD.Inertia.A*PD.Inertia.b);
PD.Aero.Cldelta_a =  Va*PD.AerC.aa;

% Force Aerodynamic coefficients
PND.Aero.Full =  PD.Aero.Full;


PND.Aero.Cmdelta_e = PD.Aero.Cmdelta_e;
PND.Aero.Cydelta_r = PD.Aero.Cydelta_r;
PND.Aero.Cldelta_r = PD.Aero.Cldelta_r;
PND.Aero.Cndelta_r = PD.Aero.Cndelta_r;
PND.Aero.Cldelta_a = PD.Aero.Cldelta_a;

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

PND.Aero.vt      = PD.Aero.Vref/sqrt(PD.Env.g*PD.Tether.L);            % V_ref/sqrt(g*L0)     

% Aerodynamic Model Limits (only for postprocess checking purposes)
PND.Aero.alfa_s =  PD.Aero.alfa_s*pi/180;                               % Stall angle(rad)
PND.Aero.beta_m =  PD.Aero.beta_m*pi/180;   


% Tether parameters
Area_t    = pi*(PD.Tether.dt /2)^2;
masa_t    = PD.Tether.ro*Area_t*PD.Tether.L;

PND.Tether.Sigma = (masa_t/PD.Inertia.m);         % sigma = m_tether/m_kite
PND.Tether.Xi    = (PD.Tether.Cd /2)*(PD.Env.ro*PD.Tether.dt*PD.Tether.L^2/PD.Inertia.m); % chi = (Cd/2)*(ro_air*d_t*L0^2/m_kite) 
PND.Tether.Ups   = [0 0 0;0 1/12 0; 0 0 1/12];

% Control parameters
PND.Control.Type      = PD.Control.Type;       % 0 ->  No control
                                               % 1 ->  Periodic control
                                               % 2 ->  Target variables
                                               
PND.Control.T0        = PD.Control.T0;         % Control period         

PND.Control.eps_lt    = PD.Control.eps_lt;    % lt = lt0*(1+eps_lt*sin(ome_lt*t))
PND.Control.ome_lt    = PD.Control.ome_lt;    % 

PND.Control.eps_lb    = PD.Control.eps_lb;    % lb = lb0*(1+eps_lb*sin(ome_lb*t))
PND.Control.ome_lb    = PD.Control.ome_lb;    % 

PND.Control.eps_del   = PD.Control.eps_del;    % delta = delta0*(1+eps_del*sin(ome_del*t))
PND.Control.ome_del   = PD.Control.ome_del;    % 

PND.Control.eps_eta   = PD.Control.eps_eta;    % eta   = eps_eta*sin(ome_eta*t))
PND.Control.ome_eta   = PD.Control.ome_eta;    % 


PND.Control.TL   = PD.Control.TL;            %
PND.Control.TC   = PD.Control.TC;            %

% Close-Loop Control Parameters
PND.CL.Target.Type    =  PD.CL.Target.Type;                     % Type of Target trajectory 
PND.CL.Target.Control =  PD.CL.Target.Control;                  % Control mean 0-> Aerodynamic surfaces, 1-> bridle control

% Target state and control variables (used for closed loop or equilibrium calculations)
PND.Target.Ome_Rotor = PD.Target.Ome_Rotor*(2*pi/60)*sqrt(PD.Tether.L/PD.Env.g); % Target angular velocity of the rotors (-)
PND.Target.xi_Rotor  = PD.Target.xi_Rotor/(PD.Env.g*PD.Inertia.m*PD.Tether.L);   % Target motor controller torque (-)
PND.Target.delta_a   = PD.Target.delta_a*pi/180;                                 % Target deflection of the aileron (rad)
PND.Target.delta_r   = PD.Target.delta_r*pi/180;                                 % Target deflection of the rudder (rad)
PND.Target.delta_e   = PD.Target.delta_e*pi/180;                                 % Target deflection of the elevator (rad)

%% Numerical Parameters
PND.Num.N       = PD.Num.N;                % Number of bars
PND.Num.IntTol  = PD.Num.IntTol;           % Integrator Tolerance
PND.Num.NewTol  = PD.Num.NewTol;           % Newton-Raphson Tolerance
PND.Num.DT      = PD.Num.DT;               % Time step
PND.Num.NDT     = PD.Num.NDT;              % Number of time steps
PND.Num.MaxIt   = PD.Num.MaxIt;            % Newton-Raphson Maximum number of iterations
PND.Num.dh      = PD.Num.dh ;              % Numerical Jacobian step

end
