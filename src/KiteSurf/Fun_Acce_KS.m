function [ak alfa_BE ] = Fun_Acce_KS(RBE,vg,xs_p,xs_pp,omega_BE,Phi,Phi_xs,Ups_s,Ups_s_xs)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Acceleration and angular accelerations                         %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                            %%
% Inputs:  RBE             -> Rotation matrix                                %%
%          vg              -> Velocity of the kite (Earth components)        %% 
%          xs_p, xs_pp     -> Time-derivatives of the kite state vector      %%
%          omega_BE        -> SB-SE angular velocity (Earth Components)      %%            
%          Phi,Phi_xs      -> Kinematic matrices (see Fun_Matrix_Omega)      %%
%          Ups_s,Ups_s_xs, -> Kinematic matrices (see Fun_Matrix_Upsilon)    %% 
%                                                                            %%
% Outputs: ak         -> Kite acceleration  (SE Components)                  %% 
%          alfa_BE    -> SB-SE angular acceleration (SE Components)          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%         Auxiliary Calculations                                                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute kite acceleration vector components in the Earth Frame -> ak
Aux1            = zeros(3,5); % (partial Upsilon_s/partial xs)*dot{xs}
for k=1:1:5
    Aux1(:,:) = Aux1(:,:) + squeeze(Ups_s_xs(:,:,k))*xs_p(k);
end
ak  = RBE'*(Ups_s*xs_pp+Aux1*xs_p)+cross(omega_BE,vg);

% Compute kite angular acceleration vector components in the Earth frame -> alfa_BE
Aux0 = zeros(3,5); % (Partial Phi/partial xs)*xsp
for k=1:1:5
    Aux0 = Aux0 + Phi_xs(:,:,k)*xs_p(k);
end
alfa_BE = RBE'*(Phi*xs_pp + Aux0*xs_p);  % Note that cross(omega,omega)=0


end
