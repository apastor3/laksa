function [QR FA_R]   = Generalized_QR_KF(t,rR,VR,SR,xs_str,xc_str,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Generalized forces upon the rods                               %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%     t              -> Dimensionless time                                %%
%     rR             -> Center of mass Position vector of the rods        %% 
%                       (SE components)                                   %%
%     VR             -> Velocity of the rods (SE components)              %%
%     SR             -> Kinematic matrix                                  %%
%     xs_str         -> State vector vector                               %%
%     xc_str         -> Control vector                                    %%
%     PND            -> Dimensionless parameters                          %%
% Outputs                                                                 %%
%     QR             -> Generalize forces                                 %%
%     FA_R           -> Aerodynamic Forces upon the rods (SE components)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators
Nc0 = 4;           % Number of control variables affecting the kinematics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the position, velocity and angular velocity of the bars, projected in the Earth Frame
Ri     = zeros(3,NR+1);           % Positions vectors of the nodes
Vi     = zeros(3,NR);             % Velocity vectors of the center of mass
ti     = zeros(3,NR);             % Tangent vectors 
QR     = zeros(2*NR+3+NG,1);      % Generalized force
for i=1:1:NR
    % Tangent Vector
    Ri(:,i+1) = Ri(:,i)-xc_str.lt*[cos(xs_str.gamma(i))*cos(xs_str.varphi(i));cos(xs_str.gamma(i))*sin(xs_str.varphi(i));sin(xs_str.gamma(i))];
    ti        = Ri(:,i+1)-Ri(:,i);
    ti        = ti/sqrt(ti'*ti); 
    % Wind velocity
    Vw        = Fun_Wind(t,rR(:,i),PND); % Earth frame components 
    % Aerodynamic Velocity    
    VA_R      = VR(:,i)-Vw;
    % Aerodynamic velocity normal to the bar
    Vni       =  VA_R-( VA_R'*ti)*ti;
    % Force
    FA_R(:,i)   = -PND.Tether.Xi*xc_str.lt*sqrt(Vni'*Vni)*Vni;
    % Generalized Force
    QR        = QR+ (FA_R(:,i)'*squeeze(SR(:,:,i)))';
    
end 


end
