function [f m alfa beta]       = Aerokite(VA,omega,PND,delta_a,delta_r,delta_e)
         

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute aerodynamic force and torque upon the kite             %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  VA    -> components of the dimensionless airspeed       %%
%                   vector in the body frame                       %%
%          omega -> kite angular velocity components in body frame %% 
%          PND   -> dimensionless parameters of the system         %%
%          delta_a  -> aileron deflection (rad)                    %%
%          delta_r  -> rudder deflection (rad)                     %%
%          delta_e  -> elevator deflection (rad)                   %%
%                                                                  %%
% Outputs: f     -> components of the aerodynamic force of the     %%
%                   kite in the body frame                         %%
%          m     -> components of the aerodynamic torque about the %%
%                   center of mass of the kite  in the body frame  %%
%          alfa  -> angle of attack   (rad)                        %% 
%          beta  -> sideslip angle    (rad)                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
% Recover the parameters
 
mu      = PND.Kite.mu; 
eps_c   = PND.Kite.c;      % eps_c
eps_b   = PND.Kite.b;      % eps_b

Cx0     = PND.Aero.Cx0;    % Cx0                          (-)
Cx_alfa = PND.Aero.Cxalfa; % Cx_alfa                      (-)
Cy_beta = PND.Aero.Cybeta; % Cy_beta                      (-)
Cz0     = PND.Aero.Cz0;    % Cz0                          (-)
Cz_alfa = PND.Aero.Czalfa; % Cz_alfa                      (-)

Cl_beta = PND.Aero.Clbeta; % Cl_beta                      (-)
Cl_p    = PND.Aero.Clp;    % Cl_p                         (-)
Cm0     = PND.Aero.Cm0;    % Cm0                          (-)
Cm_alfa = PND.Aero.Cmalfa; % Cm_alfa                      (-)
Cm_q    = PND.Aero.Cmq;    % Cm_q                         (-)
Cn_beta = PND.Aero.Cnbeta; % Cn_beta                      (-)
Cn_r    = PND.Aero.Cnr;    % Cn_r                         (-)

Vref    = PND.Aero.vt ;    % Vref/sqrt(g*L0)   
 

% Control surfaces
if isfield(PND.Aero,'Cmdelta_e') ==1
    Cm_delta_e = PND.Aero.Cmdelta_e;
else
    Cm_delta_e = 0;
end
if isfield(PND.Aero,'Cydelta_r') ==1
    Cy_delta_r = PND.Aero.Cydelta_r;
else
    Cy_delta_r = 0;
end
if isfield(PND.Aero,'Cldelta_r') ==1
    Cl_delta_r = PND.Aero.Cldelta_r;
else
    Cl_delta_r = 0;
end
if isfield(PND.Aero,'Cndelta_r') ==1
    Cn_delta_r = PND.Aero.Cndelta_r;
else
    Cn_delta_r = 0;
end
if isfield(PND.Aero,'Cldelta_a') ==1
    Cl_delta_a = PND.Aero.Cldelta_a;
else
    Cl_delta_a = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Square of the aerodynamic velocity
VA2    = VA'*VA;
%% Attack and sideslip angles
alfa   = atan(VA(3)/VA(1));
beta   = asin(VA(2)/sqrt(VA2));

%% Angular velocity in standard flight mechanics normalization
p_tilde = eps_b*omega(1)/(2*Vref);
q_tilde = eps_c*omega(2)/(Vref);
r_tilde = eps_b*omega(3)/(2*Vref);
   
% Aerodynamic force and torque about the center of mass
f(1,1)   =  mu*VA2*(Cx0+Cx_alfa*alfa);
f(2,1)   =  mu*VA2*Cy_beta*beta;
f(3,1)   =  mu*VA2*(Cz0+Cz_alfa*alfa);

m(1,1)   =  mu*VA2*eps_b*(Cl_beta*beta + Cl_p*p_tilde + Cl_delta_a*delta_a +Cl_delta_r*delta_r);
m(2,1)   =  mu*VA2*eps_c*(Cm0+Cm_alfa*alfa  + Cm_q*q_tilde + Cm_delta_e*delta_e);
m(3,1)   =  mu*VA2*eps_b*(     Cn_beta*beta + Cn_r*r_tilde + Cn_delta_r*delta_r);


end