function xc_amp  = Fun_Control_KF(t,xs_amp,PND)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Control Laws                                                   %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%           t        -> Normalized Time                                   %%
%           PND      -> Dimensonless Parameters                           %%
%           xs_amp   -> Lagrangian or Hamiltonian state vector            %%
% Outputs                                                                 %%
%           xc_amp   -> Extended control vector xc_amp = [xc xc_p xc_pp]  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

NR       = PND.Num.N ;                 % Number of Rods
NG       = PND.Gen.Num;                % Number of Generators
NT       = 4+NG+3;                     % Control variables
Nvar     = 2*PND.Num.N+3;              % Number of variables 
Nvar_p   = 2*PND.Num.N+3+PND.Gen.Num;  % Number of variables (dot) 
 
% Initialize
xc_amp = zeros(3*NT,1);

% Main Tetherand bridle 
xc_amp(1,1) = 1/PND.Num.N;      % Bar length
xc_amp(2,1) = PND.Bridle.lb;    % Bridle length
xc_amp(3,1) = PND.Bridle.delta; % Bridle angle delta
xc_amp(4,1) = PND.Bridle.eta;   % Bridle angle eta
% Generators with target angular velocity
for i=1:1:PND.Gen.Num
    xc_amp(4+i,1) = PND.Target.xi_Rotor;   % Motor controller Torques
end

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PND.Control.Type==1
   
    % Variables
    xc_amp(1,1)  =      1/PND.Num.N*( 1 + PND.Control.eps_lt*sin(PND.Control.ome_lt*t));    % Bar length    
    xc_amp(4,1)  =                        PND.Control.eps_eta*sin(PND.Control.ome_eta*t) ;  % Bridle angle eta
    
    % First derivatives
    xc_amp(NT+1,1)  =      1/PND.Num.N*PND.Control.ome_lt*PND.Control.eps_lt*cos(PND.Control.ome_lt*t);      % Bar length 
    xc_amp(NT+4,1)  =                        PND.Control.ome_eta*PND.Control.eps_eta*cos(PND.Control.ome_eta*t);   % Bridle angle eta

    % Second derivatives
    xc_amp(2*NT+1,1)  =      -1/PND.Num.N*PND.Control.ome_lt^2*PND.Control.eps_lt*sin(PND.Control.ome_lt*t);      % Bar length  
    xc_amp(2*NT+4,1) =                      -PND.Control.ome_eta^2*PND.Control.eps_eta*sin(PND.Control.ome_eta*t);   % Bridle angle eta 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PND.Control.Type==2
     
    xc_amp(4+NG+1) = PND.Target.delta_a;    % Target deflection of the aileron (rad)
    xc_amp(4+NG+2) = PND.Target.delta_r;    % Target deflection of the rudder (rad)
    xc_amp(4+NG+3) = PND.Target.delta_e;    % Target deflection of the elevator (rad)
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PND.Control.Type==3
    xc_amp(4+NG+1) = PND.Target.delta_a;    % Target deflection of the aileron (rad)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PND.Control.Type==4% Reel-in maneuver
   
   % l(t) 
   xc_amp(1,1)       =   1/PND.Num.N*( 1 + PND.Control.eps_lt*t);    % Bar length     
   % First derivatives
   xc_amp(NT+1,1)    =   1/PND.Num.N*PND.Control.eps_lt;     % Derivative of Bar length 
   % Second derivatives
   xc_amp(2*NT+1,1)  =   0;  

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PND.Control.Type==5% Close-Loop -> Take the control variable from the state vector
    
    if PND.CL.Target.Control == 0 % Use elevator, rudder and ailerons
       xc_amp(4+PND.Gen.Num+1,1) = xs_amp(Nvar+Nvar_p+1);  % Aileron Deflection
       xc_amp(4+PND.Gen.Num+2,1) = xs_amp(Nvar+Nvar_p+2);  % Rudder Deflection
       xc_amp(4+PND.Gen.Num+3,1) = xs_amp(Nvar+Nvar_p+3);  % Elevator Deflection

    else % Use the bridle
        xc_amp(2,1) = xs_amp(Nvar+Nvar_p+1);    % Bridle length
        xc_amp(3,1) = xs_amp(Nvar+Nvar_p+2);    % Bridle angle delta
        xc_amp(4,1) = xs_amp(Nvar+Nvar_p+3);    % Bridle angle eta 
    end
    
    
end

if PND.Control.Type==6 % Open-Loop Make Figure of Eight
    
   TC    = PND.Control.TC;    % Circular  Time
   TL    = PND.Control.TL;    % Linear Time
   
   etaL  = PND.Control.eps_eta; % Phi angle during at the straight segments
   
   eta_t = 2*etaL/TL;
   
   T  = 2*(TL+TC);
   
   if t>T
       t = t-floor(t/T)*T;
   end
   
   if t<TC
        eta    = -etaL;
        eta_p  =   0;
   else  
       if t<TL+TC
            eta    = -etaL + eta_t*(t-TC);
            eta_p  =  eta_t;
       else
           if t<2*TC+TL
                eta    =  etaL;
                eta_p  =   0;
           else
                eta    =  etaL - eta_t*(t-(2*TC+TL));
                eta_p  =  -eta_t;
           end
       end
   end
   
   xc_amp(4,1)     = eta;  % Bridle angle eta
   xc_amp(NT+4,1)  = eta_p;   % Bridle angle eta

end
   
   
end