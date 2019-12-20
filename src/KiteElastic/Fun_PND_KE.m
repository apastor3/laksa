function PND = Fun_PND_KE(PD)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Dimensionless Parameters                                       %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: structure PD with the physical parameters                  %%
%                                                                    %%
% Outputs: structure PND with the dimensionless parameters           %%
%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PND.Kite.Num     =  PD.Kite.Num;   


for i=1:1:PND.Kite.Num
    
    PND.Kite.mu(i)       =  PD.Env.ro*PD.Kite.A(i)*PD.Tether.L0/(2*PD.Kite.m(1));   % mu
   
    PND.Kite.Sigma(i)    =  PD.Kite.m(i)/PD.Kite.m(1);                              % Mass ratio
   
    PND.Kite.c(i)        =  PD.Kite.c(i)/PD.Tether.L0;                              % eps_c
    PND.Kite.b(i)        =  PD.Kite.b(i)/PD.Tether.L0;                              % eps_b
    PND.Kite.h(i)        =  PD.Kite.h(i)/PD.Tether.L0;                              % eps_h (only plotting purposes)
    PND.Kite.hg(i)       =  PD.Kite.hg(i)/PD.Tether.L0;                             % eps_g (only plotting purposes)


    iota_Ix              = PD.Kite.Ix(i)/(PD.Kite.m(1)*PD.Tether.L0^2);                % Ix/m *L0^2   
    iota_Iy              = PD.Kite.Iy(i)/(PD.Kite.m(1)*PD.Tether.L0^2);                % Iy/m *L0^2   
    iota_Iz              = PD.Kite.Iz(i)/(PD.Kite.m(1)*PD.Tether.L0^2);                % Iz/m *L0^2    
    iota_Ixz             = PD.Kite.Ixz(i)/(PD.Kite.m(1)*PD.Tether.L0^2);               % Ixz/m *L0^2 

    PND.Kite.iota(:,:,i) = [iota_Ix 0 iota_Ixz; 0 iota_Iy 0;iota_Ixz 0  iota_Iz];

    % Aerodynamic Force and Torque Coefficients
    PND.Aero.Full(i) =  PD.Aero.Full(i);

    if PD.Aero.Full(i) ==  1
        PND.Aero.CX(i) = PD.Aero.CX(i);
        PND.Aero.CY(i) = PD.Aero.CY(i);
        PND.Aero.CZ(i) = PD.Aero.CZ(i);

        PND.Aero.Cm(i) = PD.Aero.Cm(i);
        PND.Aero.Cl(i) = PD.Aero.Cl(i);
        PND.Aero.Cn(i) = PD.Aero.Cn(i);
    else    

        PND.Aero.CX(:,:,i) = [0  PD.Aero.Cxalfa(i) PD.Aero.Cx0(i)...    % CX0 
                       0         0              0 ...         % CX_beta
                       0         0              0 ...         % CX_p
                       0         0              0 ...         % CX_q
                       0         0              0 ...         % CX_r 
                       0         0              0 ...         % CX_delta_aileron
                       0         0              0 ...         % CX_delta_elevator
                       0         0              0];           % CX_delta_rudder     

        PND.Aero.CY(:,:,i) =  [0         0              0 ...         % CY0 
                       0         0        PD.Aero.Cybeta(i) ...  % CY_beta
                       0         0              0  ...        % CY_p
                       0         0              0  ...        % CY_q
                       0         0              0  ...        % CY_r 
                       0         0              0  ...        % CY_delta_aileron
                       0         0       PD.Aero.Cydelta_r(i)... % CY_delta_elevator
                       0         0              0  ];         % CY_delta_rudder     

        PND.Aero.CZ(:,:,i) = [ 0    PD.Aero.Czalfa(i)  PD.Aero.Cz0(i) ...   % CZ0 
                       0         0              0 ...         % CZ_beta
                       0         0              0 ...         % CZ_p
                       0         0              0 ...         % CZ_q
                       0         0              0 ...         % CZ_r 
                       0         0              0 ...         % CZ_delta_aileron
                       0         0              0 ...         % CZ_delta_elevator
                       0         0              0];           % CZ_delta_rudder     


        PND.Aero.Cm(:,:,i) = [ 0      PD.Aero.Cmalfa(i)  PD.Aero.Cm0(i) ... % Cm0 
                       0         0              0 ...         % Cm_beta
                       0         0              0 ...         % Cm_p
                       0         0            PD.Aero.Cmq(i)...  % Cm_q
                       0         0              0 ...         % Cm_r 
                       0         0              0 ...         % Cm_delta_aileron
                       0         0       PD.Aero.Cmdelta_e(i) ...% Cm_delta_elevator
                       0         0              0];           % Cm_delta_rudder     

         PND.Aero.Cl(:,:,i)= [ 0         0              0 ...         % Cl0 
                       0         0         PD.Aero.Clbeta(i) ... % Cl_beta
                       0         0            PD.Aero.Clp(i) ... % Cl_p
                       0         0              0  ...        % Cl_q
                       0         0              0 ...         % Cl_r 
                       0         0      PD.Aero.Cldelta_a(i) ... % Cl_delta_aileron
                       0         0              0 ...         % Cl_delta_elevator
                       0         0      PD.Aero.Cldelta_r(i)];   % Cl_delta_rudder   

         PND.Aero.Cn(:,:,i)= [ 0         0             0 ...         % Cn0 
                       0         0         PD.Aero.Cnbeta(i) ... % Cn_beta
                       0         0              0         ... % Cn_p
                       0         0              0  ...        % Cn_q
                       0         0           PD.Aero.Cnr(i)...   % Cn_r 
                       0         0              0 ... % Cn_delta_aileron
                       0         0              0 ...         % Cn_delta_elevator
                       0         0      PD.Aero.Cndelta_r(i) ];  % Cn_delta_rudder   

         PND.Aero.vt(i)      = PD.Aero.Vref(i)/sqrt(PD.Env.g*PD.Tether.L0);            % V_ref/sqrt(g*L0)                 
    end                           



    % Aerodynamic Model Limits (only for postprocess checking purposes)
    PND.Aero.alfa_s(i)  =  PD.Aero.alfa_s(i)*pi/180;                               % Stall angle(rad)
    PND.Aero.beta_m(i)  =  PD.Aero.beta_m(i)*pi/180;                               % Maximum sideslip angle (rad)
end
    
% Enviromental Properties
PND.Env.Type     = PD.Env.Type;                                         % Wind law
PND.Env.vw       = PD.Env.Vw/sqrt(PD.Env.g*PD.Tether.L0);               % V_0 tilde 
PND.Env.alfa     = PD.Env.alfa;                                         % Exponent of the wind speed law
PND.Env.H0       = PD.Env.H0/PD.Tether.L0;                              % Height Scale
PND.Env.eps      = PD.Env.eps;                                          % Wind Speed fluctuation level
PND.Env.Omega    = PD.Env.Omega*sqrt(PD.Tether.L0/PD.Env.g);            % Normalized fluctuation frequency 

% Auxiliary calculation 
At               = pi*(PD.Tether.dt/2)^2;                               % Tether cross section
% Common Tether Characteristics  
PND.Tether.Num   = PD.Tether.Num;                                       % Total number of tethers  
PND.Tether.nu    = PD.Tether.E*At/(PD.Kite.m(1)*PD.Env.g);              % nu = E_t A_t/mk*g 
PND.Tether.ft    = PD.Tether.ft*sqrt(PD.Env.g/PD.Tether.L0);            % Tether damping coefficient ft
  
for i=1:1:PND.Tether.Num  

  PND.Tether.Mass(i)  = PD.Tether.Mass(i);                 % Number of masses used to discretize the tether

  PND.Tether.seda0(i) = (PND.Tether.Mass(i)+1)*PD.Tether.L0/PD.Tether.L00(i);  % seda(t) = (NP+1)*L0/L*(1+PND.Ctr.eps_L*cos(PND.Ctr.Om_seda*t) ) 
  PND.Ctr.eps_L(i)    = PD.Ctr.eps_L(i);                                
  PND.Ctr.Om_L(i)     = PD.Ctr.Om_L(i)*sqrt(PD.Tether.L0/PD.Env.g);
  
  PND.Tether.Down(i)  = PD.Tether.Down(i);                 % Number of the Kite attached to the Lowest tether point (0 means the Earth System)
  PND.Tether.Dx(i)    = PD.Tether.Dx(i)/PD.Tether.L0;      % X coordinate of the down-point (SK coordinates)
  PND.Tether.Dy(i)    = PD.Tether.Dy(i)/PD.Tether.L0;      % Y coordinate of the down-point (SK coordinates)
  PND.Tether.Dz(i)    = PD.Tether.Dz(i)/PD.Tether.L0;      % Z coordinate of the down-point (SK coordinates)
  PND.Tether.Up(i)    = PD.Tether.Up(i);                   % Number of the Kite attached to the Upper tether point
  PND.Tether.Ux(i)    = PD.Tether.Ux(i)/PD.Tether.L0;      % X coordinate of the Up-point (SK coordinates)
  PND.Tether.Uy(i)    = PD.Tether.Uy(i)/PD.Tether.L0;      % Y coordinate of the Up-point (SK coordinates)
  PND.Tether.Uz(i)    = PD.Tether.Uz(i)/PD.Tether.L0;      % Z coordinate of the Up-point (SK coordinates)
  
  ATT                 =  PD.Tether.L00(i)*PD.Tether.dt;                                % Total tether transversal Area
  MTT                 =  acos(-1)*(PD.Tether.dt/2)^2*PD.Tether.L00(i)*PD.Tether.ro;    % Total tether mass
  Sigma0              =  MTT/PD.Kite.m(1);                                             % Mass ratio
  Xi0                 =  0.5*PD.Env.ro*ATT*PD.Tether.L0*PD.Tether.Cd/PD.Kite.m(1);     % Aerodynamic ratio
 
  for j=1:1:PND.Tether.Mass(i)
      if PND.Tether.Mass(i)==1
         PND.Mass.Sigma(i,j)  = Sigma0;                                     
         PND.Mass.Xi(i,j)     = Xi0;
      else
         if PND.Tether.Mass(i)==2
             PND.Mass.Sigma(i,j)  = Sigma0/2;                                     
             PND.Mass.Xi(i,j)     = Xi0/2;
         else
             if j==1 || j==PND.Tether.Mass(i)
                 PND.Mass.Sigma(i,j)  = 1.5*Sigma0/(PND.Tether.Mass(i)+1);                                     
                 PND.Mass.Xi(i,j)     =    1.5*Xi0/(PND.Tether.Mass(i)+1);
             else
                 PND.Mass.Sigma(i,j)  = Sigma0/(PND.Tether.Mass(i)+1);                                     
                 PND.Mass.Xi(i,j)     =    Xi0/(PND.Tether.Mass(i)+1);
             end
         end
      end
  end
end

for i=1:1:PD.Kite.Num    % Total Number of aircraft
                                        % delta = delta0 + delta1*cos(Om_delta*t)
    PND.Ctr.delta_a0(i)   = PD.Ctr.delta_a0(i)*pi/180;           % Aileron deflection (rad)
    PND.Ctr.delta_r0(i)   = PD.Ctr.delta_r0(i)*pi/180;           % Rudder deflection (rad)
    PND.Ctr.delta_e0(i)   = PD.Ctr.delta_e0(i)*pi/180;           % Elevator deflection (rad)
    
    PND.Ctr.delta_a1(i)   = PD.Ctr.delta_a1(i)*pi/180;           % Aileron deflection (rad)
    PND.Ctr.delta_r1(i)   = PD.Ctr.delta_r1(i)*pi/180;           % Rudder deflection (rad)
    PND.Ctr.delta_e1(i)   = PD.Ctr.delta_e1(i)*pi/180;           % Elevator deflection (rad)
    
    PND.Ctr.Om_delta_a(i) = PD.Ctr.Om_delta_a(i)*sqrt(PD.Tether.L0/PD.Env.g);         % Angular frequency of Aileron deflection (-)
    PND.Ctr.Om_delta_r(i) = PD.Ctr.Om_delta_r(i)*sqrt(PD.Tether.L0/PD.Env.g);         % Angular frequency of Rudder deflection (-)
    PND.Ctr.Om_delta_e(i) = PD.Ctr.Om_delta_e(i)*sqrt(PD.Tether.L0/PD.Env.g);         % Angular frequency of Elevator deflection (-) 
end

% Total number of masses of the simulators
PND.Tether.MassT      = sum(PND.Tether.Mass);

% Numerical Parameters
PND.Num.RelTol  = PD.Num.RelTol;           % Integrator Relative Tolerance
PND.Num.AbsTol  = PD.Num.AbsTol;           % Integrator Absolute Tolerance
PND.Num.NewTol  = PD.Num.NewTol;           % Newton-Raphson Tolerance
PND.Num.MaxIt   = PD.Num.MaxIt;            % Newton-Raphson Maximum number of iterations
PND.Num.dh      = PD.Num.dh;               % Numerical Jacobian step


end