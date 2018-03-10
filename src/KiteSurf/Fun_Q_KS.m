function [Q f m alfa beta fa ma ] = Fun_Q_KS(vw,vg,omega,Ups_s,Phi,T_Bp,T_Bm,m_Bp,m_Bm,PND)
         

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Generalized forces                                             %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  vw    -> dimensionless wind velocity in body axes              %%
%         vg    -> dimensionless center of mass velocity in body axes    %%
%         omega -> dimensionless absolute angular velocity in body axes  %%
%         Ups_s -> Matrix (see function Fun_Matrix_Upsilon)              %% 
%         Phi   -> Matrix (see function Fun_Matrix_Omega)                %%
%         T_Bp  -> Tether force upon the kite at B Plus                  %%
%         T_Bm  -> Tether force upon the kite at B minus                 %%  
%         m_Bp  -> Torque about the center of mass due to T_Bp           %%
%         m_Bm  -> Torque about the center of mass due to T_Bm           %%
%         PND   -> Dimensionless parameter                               %%
%                                                                        %% 
% Outputs Q     -> Generalized forces                                    %%
%         f     -> total forces                                          %%
%         m     -> total torque about G                                  %%
%         alfa  -> angle of attack                                       %%
%         beta  -> sideslip angle                                        %%
%         fa    -> aerodynamic force                                     %%
%         ma    -> aerodynamic torque about G                            %%
%                                                                        %%
% % Note: T_Bp, T_Bm, m_Bp, m_Bm, f, m, fa, ma components are body frame %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This subroutine computes the generalized forces Q_i                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Aerodinamic velocity in the body frame
VA      = vg-vw;

%% Tension difference -> incorporated as a rudder deflection into the model 
delta_r    = sqrt(T_Bp'*T_Bp)-sqrt(T_Bm'*T_Bm);
%% Aerodynamic Forces and torques upon the kite
[fa ma alfa beta] = Aerokite(VA,omega,PND,0,delta_r,0);
% We sum the contribution of the tethers
f       = fa + T_Bp + T_Bm;
m       = ma + m_Bp + m_Bm;

%% Generalized forces
Q        = f'*Ups_s+m'*Phi;
Q        = Q';

end