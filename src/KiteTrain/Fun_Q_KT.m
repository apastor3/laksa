function [Q f m alfa beta] = Fun_Q_KT(xc,vw,vg,omega,Ups,Phi,RBE,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and  Jose A. Serrano-Iglesia           %
% Language  : Matlab                                                         %
% Synopsis  : Generalized forces                                             %
% Copyright : Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %%
% Inputs:  xc    -> control vector                                    %% 
%          vw    -> components of the dimensionless wind velocity     %%
%                   vectors in the body frame                         %%
%          vg    -> components of the dimensionless velocity vector   %% 
%                   of the kite in the body frame                     %%
%          omega -> kite angular velocity components in body frame    %% 
%          Ups   -> see Fun_Upsilon_KT                                %%
%          Phi   -> see Fun_Omega_KT                                  %%
%          RBE   -> SB-SE rotation matrices                           %%
%          PND   -> dimensionless parameters of the system            %%
%                                                                     %%
% Outputs: Q     -> generalized forces                                %%
%          f     -> components of the aerodynamic forces of the       %%
%                   kites in the body frames                          %%
%          m     -> components of the aerodynamic torques about the   %%
%                   center of mass of the kites  in the body frames   %%
%          alfa  -> angles of attack  (rad)                           %% 
%          beta  -> sideslip angles   (rad)                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f    = zeros(3,PND.Kite.N);
m    = zeros(3,PND.Kite.N);
alfa = zeros(PND.Kite.N,1);
beta = zeros(PND.Kite.N,1);
Q    = zeros(4*PND.Kite.N,1);
for i=1:1:PND.Kite.N
   %% Components of the aerodinamic velocity vector in the body frame
   VA     = vg(:,i)-vw(:,i);
   
   %% Components in the body frame of the aerodynamic force and torque 
   PND0.Kite.mu     = PND.Kite.mu(i); 
   PND0.Kite.c      = PND.Kite.c(i);      % eps_c
   PND0.Kite.b      = PND.Kite.b(i);      % eps_b

   % Recover the parameters
   PND0.Aero.CX   = PND.Aero.CX;
   PND0.Aero.CY   = PND.Aero.CY;
   PND0.Aero.CZ   = PND.Aero.CZ;
    
   PND0.Aero.Cm   = PND.Aero.Cm;
   PND0.Aero.Cl   = PND.Aero.Cl;
   PND0.Aero.Cn   = PND.Aero.Cn;

   PND0.Aero.vt   = PND.Aero.vt;           % Vref/sqrt(g*L0)  
   PND0.Aero.Full = PND.Aero.Full;           % Vref/sqrt(g*L0)  
   
   
   [f(:,i) m(:,i) alfa(i,1) beta(i,1)] = Aerokite(VA,omega(:,i),PND0,xc(1,i),xc(2,i),xc(3,i));
   
   %% Components of the aerodynamic force in the Earth frame
   FAE      = squeeze(RBE(:,:,i))'*f(:,i);
   %% Generalized forces
   Q        = Q + (FAE'*Ups(:,:,i) + m(:,i)'*Phi(:,:,i))';
end

end