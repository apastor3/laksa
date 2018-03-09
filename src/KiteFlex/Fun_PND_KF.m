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
% Force Aerodynamic coefficients

PND.Aero.Cx0    = PD.Aero.Cx0;            % Cx0                          (-)  % Groot
PND.Aero.Cxalfa = PD.Aero.Cxalfa;         % Cx_alfa                      (-)  % Groot
PND.Aero.Cybeta = PD.Aero.Cybeta;         % Cy_beta                      (-)  % Williams
PND.Aero.Cz0    = PD.Aero.Cz0;            % Cz0                          (-)  % Groot
PND.Aero.Czalfa = PD.Aero.Czalfa;         % Cz_alfa                      (-)  % Groot

% Torque Aerodynamic coefficients
PND.Aero.Clbeta = PD.Aero.Clbeta;         % Cl_beta                      (-)       % Williams
PND.Aero.Clp    = PD.Aero.Clp;            % Cl_p_tilde                   (-)       % Williams Clp = Clp * b/2Va
PND.Aero.Cm0    = PD.Aero.Cm0;            % Cm0                          (-)       % Groot
PND.Aero.Cmalfa = PD.Aero.Cmalfa;         % Cm_alfa                      (-)       % Groot
PND.Aero.Cmq    = PD.Aero.Cmq;            % Cm_q_tilde                   (-)       % Williams     Cmq = Cmq*c/Va  
PND.Aero.Cnbeta = PD.Aero.Cnbeta;         % Cn_beta                      (-)       % Williams
PND.Aero.Cnr    = PD.Aero.Cnr;            % Cn_r_tilde                   (-)       % Williams Cn_r = Cnr*b/2Va
PND.Aero.vt     = PD.Aero.Vref/sqrt(PD.Env.g*PD.Tether.L);  % Vref/sqrt(g*L0)           (-)

% Aerodynamic Control Surfaces
%Elevator
Vt                 = PD.AerC.St*PD.AerC.lt/(PD.Inertia.A*PD.Inertia.c);
PND.Aero.Cmdelta_e = -Vt*PD.AerC.at;
% Rudder
Vvz                =  PD.AerC.Sv*PD.AerC.hv/(PD.Inertia.A*PD.Inertia.b);
Vvx                =  PD.AerC.Sv*PD.AerC.lt/(PD.Inertia.A*PD.Inertia.b);

PND.Aero.Cydelta_r =  PD.AerC.Sv/PD.Inertia.A*PD.AerC.av;
PND.Aero.Cldelta_r =  Vvz*PD.AerC.av;
PND.Aero.Cndelta_r = -Vvx*PD.AerC.av;
% Ailerons
Va                 =  PD.AerC.Sa*PD.AerC.ya/(PD.Inertia.A*PD.Inertia.b);
PND.Aero.Cldelta_a =  Va*PD.AerC.aa;

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
