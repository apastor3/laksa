function PND = Fun_PND_KS(PD)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Dimensionless Parameters                                       %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: structure PD with the physical parameters                  %%
%                                                                    %%
% Outputs: structure PND with the dimensionless parameters           %%
%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PND.Kite.mu      =  PD.Env.ro*PD.Kite.A*PD.Tether.Ll/(2*PD.Kite.m );     % mu
PND.Kite.c       =  PD.Kite.c/PD.Tether.Ll;                              % eps_c
PND.Kite.b       =  PD.Kite.b/PD.Tether.Ll;                              % eps_b
PND.Kite.h       =  PD.Kite.h/PD.Tether.Ll;                              % eps_h (only plotting purposes)
PND.Kite.hg      =  PD.Kite.hg/PD.Tether.Ll;                             % eps_g (only plotting purposes)

iota_Ix          = PD.Kite.Ix/(PD.Kite.m*PD.Tether.Ll^2);                % Ix/m *L0^2   
iota_Iy          = PD.Kite.Iy/(PD.Kite.m*PD.Tether.Ll^2);                % Iy/m *L0^2   
iota_Iz          = PD.Kite.Iz/(PD.Kite.m*PD.Tether.Ll^2);                % Iz/m *L0^2    
iota_Ixz         = PD.Kite.Ixz/(PD.Kite.m*PD.Tether.Ll^2);               % Ixz/m *L0^2 

PND.Kite.iota    = [iota_Ix 0 iota_Ixz; 0 iota_Iy 0;iota_Ixz 0  iota_Iz];

% Force Aerodynamic coefficients
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
PND.Aero.vt        = PD.Aero.Vref/sqrt(PD.Env.g*PD.Tether.Ll);            % V_ref/sqrt(g*L0)                 

% Aerodynamic Model Limits (only for postprocess checking purposes)
PND.Aero.alfa_s  =  PD.Aero.alfa_s*pi/180;                               % Stall angle(rad)
PND.Aero.beta_m  =  PD.Aero.beta_m*pi/180;                               % Maximum sideslip angle (rad)

% Enviromental Properties
PND.Env.Type     = PD.Env.Type;                                         % Wind law
PND.Env.vw       = PD.Env.Vw/sqrt(PD.Env.g*PD.Tether.Ll);               % V_0 tilde 
PND.Env.alfa     = PD.Env.alfa;                                         % Exponent of the wind speed law
PND.Env.H0       = PD.Env.H0/PD.Tether.Ll;                              % Height Scale
PND.Env.eps      = PD.Env.eps;                                          % Wind Speed fluctuation level
PND.Env.Omega    = PD.Env.Omega*sqrt(PD.Tether.Ll/PD.Env.g);            % Normalized fluctuation frequency 

% Attachment points and control bar length 
PND.Tether.XA    = PD.Tether.XA/PD.Tether.Ll;                            % XA/Ll 
PND.Tether.YA    = PD.Tether.YA/PD.Tether.Ll;                            % YA/Ll 
PND.Tether.ZA    = PD.Tether.ZA/PD.Tether.Ll;                            % ZA/Ll 
PND.Tether.XB    = PD.Tether.XB/PD.Tether.Ll;                            % XB/Ll 
PND.Tether.YB    = PD.Tether.YB/PD.Tether.Ll;                            % YB/Ll 
PND.Tether.ZB    = PD.Tether.ZB/PD.Tether.Ll;                            % ZB/Ll 

PND.Tether.l     = sqrt(1-PND.Tether.YA^2);                              % l -> OEO2/Ll

% Auxiliary calculation 
At               = pi*(PD.Tether.dt/2)^2;                                % Tether cross section
% Tether properties 
PND.Tether.seda  = PD.Tether.Ll/PD.Tether.Lt;                            % zeda = L0/L1
PND.Tether.kappa = PD.Tether.E *At/(PD.Kite.m*PD.Env.g);                 % nu = E_t A_t/mk*g 
PND.Tether.ft    = PD.Tether.ft*sqrt(PD.Env.g/PD.Tether.Lt);             % Tether damping coefficient ft
  
% Control Bar Geometry
PND.Bar.Lc       = PD.Bar.Lc/PD.Tether.Ll;                               % Lc/Ll
PND.Bar.Ls       = PD.Bar.Ls/PD.Tether.Ll;                               % Ls/Ll  
PND.Bar.Lds      = PD.Bar.Lds/PD.Tether.Ll;                              % Lds/Ll     
PND.Bar.Lps      = PD.Bar.Lps/PD.Tether.Ll;                              % Lps/Ll 

% Control Variables
PND.Ctr.Type   = PD.Ctr.Type;                                          % Control Type
PND.Ctr.Lam    = PD.Ctr.Lam*pi/180;                                    % lambda_0      (rad)
PND.Ctr.OmLam  = PD.Ctr.OmLam/sqrt(PD.Env.g/PD.Tether.Ll);             % omega_lambda
PND.Ctr.PR0    = PD.Ctr.PR0;                                           % PR0         
PND.Ctr.PRf    = PD.Ctr.PRf;                                           % PRf 
PND.Ctr.Tr     = PD.Ctr.Tr*sqrt(PD.Env.g/PD.Tether.Ll);                % tr  
PND.Ctr.Ts     = PD.Ctr.Ts*sqrt(PD.Env.g/PD.Tether.Ll);                % ts 
PND.Ctr.TL     = PD.Ctr.TL/sqrt(PD.Env.g/PD.Tether.Ll);                % Characteristic Time in the steering manoeuvre 
PND.Ctr.TC     = PD.Ctr.TC/sqrt(PD.Env.g/PD.Tether.Ll);                % Characteristic Time in the steering manoeuvre 

PND.Ctr.GD     = PD.Ctr.GD;
PND.Ctr.GP     = PD.Ctr.GP;
PND.Ctr.GI     = PD.Ctr.GI;

PND.Ctr.Y0     = PD.Ctr.Y0/PD.Tether.Ll;                % Y0       for Close Loop
PND.Ctr.Om     = PD.Ctr.Om/sqrt(PD.Env.g/PD.Tether.Ll); % omega_Y  for Close Loop

% Numerical Parameters
PND.Num.RelTol  = PD.Num.RelTol;           % Integrator Relative Tolerance
PND.Num.AbsTol  = PD.Num.AbsTol;           % Integrator Absolute Tolerance
PND.Num.NewTol  = PD.Num.NewTol;           % Newton-Raphson Tolerance
PND.Num.MaxIt   = PD.Num.MaxIt;            % Newton-Raphson Maximum number of iterations
PND.Num.dh      = PD.Num.dh;               % Numerical Jacobian step


end