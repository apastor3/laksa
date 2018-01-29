function [QG FA_G MA_G MMC_G] = Generalized_QG_KF(t,rG,vG,SG,OmG,R_KE,xc,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Generalized forces upon the rotors                             %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved     %
%----------------------------------------------------------------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%     t          -> Normalized time                                       %%  
%     rG         -> Position vector of the rotors (SE components)         %%
%     vG         -> Velocity of the rotors (SE components)                %%
%     SG,OmG     -> Kinematic Matrices                                    %% 
%     R_KE       -> Rotation matrix                                       %%
%     xc         -> control vector                                        %%
%     PND        -> Dimensionless parameters                              %%
% Outputs                                                                 %%
%     QG         -> Generalize forces                                     %%
%     FA_G       -> Rotor Aerodynamic Forces (SK components )         %%
%     MA_G       -> Rotor Aerodynamic Moment (SK components)          %%
%     MMC_G      -> Motor controller torques (SK components)              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NR  = PND.Num.N ;   % Number of Rods
NG  = PND.Gen.Num;  % Number of Generators
Nc0 = 4;           % Number of control variables affecting the kinematics


QG = zeros(2*NR+3+NG,1);

nu    = PND.Gen.nu;
Cf    = PND.Gen.Cf;
Cm    = PND.Gen.Cm;

for i=1:1:NG  
    % Wind velocity
    Vw        = Fun_Wind(t,rG(:,i),PND); % Earth frame components 
    % Compute the aerodynamic velocity projected in the Earth frame
    VAG       = vG(:,i)-Vw;        % Components in the Earth frame
    % Compute the aerodynamic velocity components in the kite frame
    VAG       = R_KE*VAG;
    % Compute the aerodynamic velocity normal to the plane of the blades
    u_axis    = [cos(nu(i)) 0 -sin(nu(i))]';
    VAG_perp  =  u_axis'*VAG; 
    % Force and torques in the kite frame
    FA_G(:,i)   =  -PND.Gen_Chi*Cf*VAG_perp^2*u_axis;
    MA_G(:,i)   =   PND.Gen.l*PND.Gen_Chi*Cm*VAG_perp^2*u_axis;
    MMC_G(:,i)  =  -xc(4+i)*u_axis;
    % Generalized Forces   
    QG        = QG + ((R_KE'*FA_G(:,i))'*squeeze(SG(:,:,i))+( MA_G(:,i) + MMC_G(:,i) )'*OmG(:,:,i))';
    
end


if NG ==0
    FG = [0 0 0]';
    MG = [0 0 0]';
end

end