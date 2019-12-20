function [CX CY CZ Cl Cm Cn ] = Aero_Coefficients(alfa,beta,p,q,r,PND,delta_a,delta_r,delta_e)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute aerodynamic force and torque upon the kite             %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  alfa, beta -> angle of attck and sideslip angle         %%
%          p,q,r      -> normalized angular velocity components    %% 
%          PND   -> dimensionless parameters of the system         %%
%          delta_a  -> aileron deflection (rad)                    %%
%          delta_r  -> rudder deflection (rad)                     %%
%          delta_e  -> elevator deflection (rad)                   %%
%                                                                  %%
% Outputs: Cx,Cy,Cz -> aerodynamic force coeffcients               %%
%          Cl,Cm,Cn -> aerodynamic torque coefficient              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vec_alfa = [alfa^2 alfa 1]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CX  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CX0     = PND.Aero.CX(1:3)*Vec_alfa; 
CX_beta = PND.Aero.CX(4:6)*Vec_alfa;
CX_p    = PND.Aero.CX(7:9)*Vec_alfa;
CX_q    = PND.Aero.CX(10:12)*Vec_alfa;
CX_r    = PND.Aero.CX(13:15)*Vec_alfa;
CX_ail  = PND.Aero.CX(16:18)*Vec_alfa;
CX_ele  = PND.Aero.CX(19:21)*Vec_alfa;
CX_rud  = PND.Aero.CX(22:24)*Vec_alfa;

CX      = CX0 + CX_beta*beta + CX_p*p + CX_q*q  + CX_r*r + CX_ail*delta_a + CX_ele*delta_e  + CX_rud*delta_r;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CY0     = PND.Aero.CY(1:3)*Vec_alfa; 
CY_beta = PND.Aero.CY(4:6)*Vec_alfa;
CY_p    = PND.Aero.CY(7:9)*Vec_alfa;
CY_q    = PND.Aero.CY(10:12)*Vec_alfa;
CY_r    = PND.Aero.CY(13:15)*Vec_alfa;
CY_ail  = PND.Aero.CY(16:18)*Vec_alfa;
CY_ele  = PND.Aero.CY(19:21)*Vec_alfa;
CY_rud  = PND.Aero.CY(22:24)*Vec_alfa;

CY      = CY0 + CY_beta*beta + CY_p*p + CY_q*q  + CY_r*r + CY_ail*delta_a + CY_ele*delta_e  + CY_rud*delta_r;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CZ  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CZ0     = PND.Aero.CZ(1:3)*Vec_alfa; 
CZ_beta = PND.Aero.CZ(4:6)*Vec_alfa;
CZ_p    = PND.Aero.CZ(7:9)*Vec_alfa;
CZ_q    = PND.Aero.CZ(10:12)*Vec_alfa;
CZ_r    = PND.Aero.CZ(13:15)*Vec_alfa;
CZ_ail  = PND.Aero.CZ(16:18)*Vec_alfa;
CZ_ele  = PND.Aero.CZ(19:21)*Vec_alfa;
CZ_rud  = PND.Aero.CZ(22:24)*Vec_alfa;

CZ      = CZ0 + CZ_beta*beta + CZ_p*p + CZ_q*q  + CZ_r*r + CZ_ail*delta_a + CZ_ele*delta_e  + CZ_rud*delta_r;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Cl  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cl0     = PND.Aero.Cl(1:3)*Vec_alfa; 
Cl_beta = PND.Aero.Cl(4:6)*Vec_alfa;
Cl_p    = PND.Aero.Cl(7:9)*Vec_alfa;
Cl_q    = PND.Aero.Cl(10:12)*Vec_alfa;
Cl_r    = PND.Aero.Cl(13:15)*Vec_alfa;
Cl_ail  = PND.Aero.Cl(16:18)*Vec_alfa;
Cl_ele  = PND.Aero.Cl(19:21)*Vec_alfa;
Cl_rud  = PND.Aero.Cl(22:24)*Vec_alfa;

Cl      = Cl0 + Cl_beta*beta + Cl_p*p + Cl_q*q  + Cl_r*r + Cl_ail*delta_a + Cl_ele*delta_e  + Cl_rud*delta_r;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Cm  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cm0     = PND.Aero.Cm(1:3)*Vec_alfa; 
Cm_beta = PND.Aero.Cm(4:6)*Vec_alfa;
Cm_p    = PND.Aero.Cm(7:9)*Vec_alfa;
Cm_q    = PND.Aero.Cm(10:12)*Vec_alfa;
Cm_r    = PND.Aero.Cm(13:15)*Vec_alfa;
Cm_ail  = PND.Aero.Cm(16:18)*Vec_alfa;
Cm_ele  = PND.Aero.Cm(19:21)*Vec_alfa;
Cm_rud  = PND.Aero.Cm(22:24)*Vec_alfa;

Cm      = Cm0 + Cm_beta*beta + Cm_p*p + Cm_q*q  + Cm_r*r + Cm_ail*delta_a + Cm_ele*delta_e  + Cm_rud*delta_r;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Cn  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cn0     = PND.Aero.Cn(1:3)*Vec_alfa; 
Cn_beta = PND.Aero.Cn(4:6)*Vec_alfa;
Cn_p    = PND.Aero.Cn(7:9)*Vec_alfa;
Cn_q    = PND.Aero.Cn(10:12)*Vec_alfa;
Cn_r    = PND.Aero.Cn(13:15)*Vec_alfa;
Cn_ail  = PND.Aero.Cn(16:18)*Vec_alfa;
Cn_ele  = PND.Aero.Cn(19:21)*Vec_alfa;
Cn_rud  = PND.Aero.Cn(22:24)*Vec_alfa;

Cn      = Cn0 + Cn_beta*beta + Cn_p*p + Cn_q*q  + Cn_r*r + Cn_ail*delta_a + Cn_ele*delta_e  + Cn_rud*delta_r;   
