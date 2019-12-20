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
     
mu      = PND.Kite.mu; 
eps_c   = PND.Kite.c;      % eps_c
eps_b   = PND.Kite.b;      % eps_b
 

if PND.Aero.Full == 1
    Full = PND.Aero.Full;
else
    Full = 0;
    Vref = PND.Aero.vt;    % Vref/sqrt(g*L0)   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Square of the aerodynamic velocity
VA2    = VA'*VA;

if VA2 == 0
    alfa   = 0;
    beta   = 0;
    f      = zeros(3,1);
    m      = zeros(3,1); 
else
    %% Attack and sideslip angles
    alfa   = atan(VA(3)/VA(1));
    beta   = asin(VA(2)/sqrt(VA2));

    %% Angular velocity in standard flight mechanics normalization
    if Full ==0
        p_tilde = eps_b*omega(1)/(2*Vref);
        q_tilde = eps_c*omega(2)/(Vref);
        r_tilde = eps_b*omega(3)/(2*Vref);
    else
        p_tilde = eps_b*omega(1)/(2*sqrt(VA2));
        q_tilde = eps_c*omega(2)/(2*sqrt(VA2));
        r_tilde = eps_b*omega(3)/(2*sqrt(VA2));
    end
       
    [CX CY CZ Cl Cm Cn ] = Aero_Coefficients(alfa,beta,p_tilde,q_tilde,r_tilde,PND,delta_a,delta_r,delta_e);

    
    % Aerodynamic force and torque about the center of mass
    f(1,1)   =  mu*VA2*CX;
    f(2,1)   =  mu*VA2*CY;
    f(3,1)   =  mu*VA2*CZ;

    m(1,1)   =  mu*VA2*eps_b*Cl;
    m(2,1)   =  mu*VA2*eps_c*Cm;
    m(3,1)   =  mu*VA2*eps_b*Cn;
   
end

end