function [QK FA_K MA_K MMC_K alfa beta] = Generalized_QK_KF(t,rK,vK,omegaK,R_KE,SK,OmK,xc,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Generalized forces upon the kite                               %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%----------------------------------------------------------------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%              t        -> dimensionlesstime                              %%
%              rK       -> Kite position vector (SE components)           %%
%              vK       -> Kite velocity vector (SE components)           %%
%              omegaK   -> Kite Angular velocity (SK components)          %%
%              R_KE     -> SE-SK Rotation matrix                          %%
%              SK,OmK   -> Kite kinematic matrices                        %%
%              xc       -> Control Vector                                 %%
%              PND      -> Dimensionless parameters                       %%
% Outputs      QK       -> Generalized forces                             %%
%              FA_K     -> Kite Aerodynamic force (SK components)         %%
%              MA_K     -> Kite Aerodynamic moment (SK components)        %%
%              MMC_K    -> Controller Torque upon the kite (SK components)%%
%              alfa     -> Angle of attack                                %%
%              beta     -> Sideslip angle                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------------------
% Project   : KiteFlex
% Author    : Alejandro Pastor-Rodriguez, Gonzalo Sanchez-Arriaga
% Language  : Matlab
% Synopsis  : Compute generalized forces due to the kite 
% Copyright : Universidad Carlos III de Madrid, 2016. All rights reserved
% Date      : June 2017
%------------------------------------------------------------------------------
% 
% Inputs
% Vw             -> Wind velocity
% SK,CK,OmK,R_KE -> Kinematic matrices
% xs_p,xc_p      -> Derivatives of the state and the control vectors
% Outputs
% Qk             -> Generalize forces
% Fk             -> Aerodynamic Forces upon the kite (Kite frame components)
% Mk             -> Aerodynamic Moment about the center of mass (Kite frame components)
% M_MC           -> Total Motor controller torque (Kite frame components)
% vK             -> Kite velocity
% omegaK         -> Kite absolute angular velocity
% alfa           -> Angle of attack
% beta           -> sideslip angle


Nc0 = 4;           % Number of control variables affecting the kinematics
NG  = PND.Gen.Num;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xc_aero = xc(4+NG+1:end,1);

% Recover Aerodynamic control surfaces deflection
delta_a = xc_aero(1,1);
delta_r = xc_aero(2,1);
delta_e = xc_aero(3,1);

% Kite Velocity in the body frame
vK     = R_KE*vK;
% Wind velocity 
Vw  = Fun_Wind(t,rK,PND); % Earth frame components 
Vw  = R_KE*Vw;            % Kite frame components
% Compute the aerodynamic velocity of the kite projected in the Body frame
VA_K    = vK-Vw;              
% Compute Generalized Forces in the body frame
[FA_K MA_K alfa beta] = Aerokite(VA_K,omegaK,PND,delta_a,delta_r,delta_e);

% Compute reaction torque from the motor controller
if NG>0
   nu    = PND.Gen.nu;
end
MMC_K = zeros(3,1);
for i=1:1:NG
  u_axis    = [cos(nu(i)) 0 -sin(nu(i))]';
  %Moment in the kite frame
  MMC_K   = MMC_K + xc(4+i)*u_axis;  
end
% Generalized Forces
QK       = ((R_KE'*FA_K)'*SK+(MA_K + MMC_K)'*OmK)';
end