function DF_Control = Fun_PID(t,xs_amp,Euler,Euler_p,Euler_pp,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Time derivative of the Lagrangian function                     %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
% Input                                                                      %
% t                         - > Kinematics matrices and gradients            %
% xs_amp                    - > State Vector                                 %
% Euler                     - > Euler angles (Theta, Psi, Phi)               %
% Euler_p                   - > First Derivative of Euler angles             %
% Euler_pp                  - > Second Derivative of Euler angles            %
% PND                       - > Dimensionless Parameters                     %
%                                                                            %
% Output                                                                     %
% DF_Control               - > Time derivative of control variables          %
%                              Bridle or Control Surfaces                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NR           = PND.Num.N;
Nvar         = 2*PND.Num.N+3;              % Number of variables 
Nvar_p       = 2*PND.Num.N+3+PND.Gen.Num;  % Number of variables (dot) 
  
% Variables
Theta        = Euler(1);
Psi          = Euler(2);
Phi          = Euler(3);
 
Theta_p      = Euler_p(1);
Psi_p        = Euler_p(2);
Phi_p        = Euler_p(3);

Theta_pp     = Euler_pp(1);
Psi_pp       = Euler_pp(2);
Phi_pp       = Euler_pp(3);
 
% Target Variables
[xs_target xs_p_target xs_pp_target ]   = Target_Trajectories(t,PND);
         
Theta_Target    = xs_target(2*NR+1);
Psi_Target      = xs_target(2*NR+2);
Phi_Target      = xs_target(2*NR+3);

Theta_p_Target  = xs_p_target(2*NR+1);
Psi_p_Target    = xs_p_target(2*NR+2);
Phi_p_Target    = xs_p_target(2*NR+3);

Theta_pp_Target = xs_pp_target(2*NR+1);
Psi_pp_Target   = xs_pp_target(2*NR+2);
Phi_pp_Target   = xs_pp_target(2*NR+3);

% Error Functions
E_Theta      = Theta_Target-Theta;
E_Psi        = Psi_Target-Psi;
E_Phi        = Phi_Target-Phi;

E_p_Theta    = Theta_p_Target-Theta_p;
E_p_Psi      =   Psi_p_Target-Psi_p;
E_p_Phi      =   Phi_p_Target-Phi_p;
 

E_pp_Theta    = Theta_pp_Target-Theta_pp;
E_pp_Psi      =   Psi_pp_Target-Psi_pp;
E_pp_Phi      =   Phi_pp_Target-Phi_pp;


%%%%%%%%%%%%%%%%%%%%
%%%% PID Controls %%
%%%%%%%%%%%%%%%%%%%%
if PND.CL.Target.Control == 0 % Use elevator, rudder and ailerons
    % Ailerons
    Ia = 20;%20;   
    Pa = 10;
    Da = 10; 

    DF_Control(1,1) = Ia*E_Phi + Pa*E_p_Phi+Da*E_pp_Phi; 

    % Rudder    % Poner números negativos
    Ir = -20; % -20 
    Pr = -10; 
    Dr = -10;
    DF_Control(2,1) = Ir*E_Psi + Pr*E_p_Psi+Dr*E_pp_Psi; 

    % Elevator
    Ie = -10; 
    Pe =   0;
    De =   0; 
    DF_Control(3,1) = Ie*E_Theta + Pe*E_p_Theta+De*E_pp_Theta;      
    
   
    % Estimate the next value of the deflections to avoid overpass
    % saturation
    DT = 0.01;
    delta_a_DT = xs_amp(Nvar+Nvar_p+1)+DF_Control(1,1)*DT;
    delta_r_DT = xs_amp(Nvar+Nvar_p+2)+DF_Control(2,1)*DT;
    delta_e_DT = xs_amp(Nvar+Nvar_p+3)+DF_Control(3,1)*DT;

    if abs(delta_a_DT)>30*pi/180
          DF_Control(1,1) = 0;
    end
    if abs(delta_r_DT)>30*pi/180
         DF_Control(2,1) = 0;
    end
    if abs(delta_e_DT)>30*pi/180
         DF_Control(3,1) = 0;
    end
    
    
else % use the bridle (lB, delta, eta )
    
    DF_Control(1,1) = 0; % d lB/dt
    DF_Control(2,1) = 0; % d delta/dt
    
    % eta
    Ie              =  -100;
    Pe              = -1;  % -20
    De              = -0;   % -2
    DF_Control(3,1) = Ie*E_Phi + Pe*E_p_Phi+De*E_pp_Phi;   
    
end

    



end