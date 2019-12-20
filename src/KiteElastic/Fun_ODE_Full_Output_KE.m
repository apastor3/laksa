function [R_KE rK vK eulerK omegaK F_Tether M_Tether FA MA alfa beta  rM DF T_U T_D] = Fun_ODE_Full_Output_KE(t,X,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga                                        %
% Language  : Matlab                                                         %
% Synopsis  : Compute all relevant quantities                                %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:                                                                  %
%           t   -> Time                                                      %
%           X   -> Extended state vector                                     %
%           PND -> Dimensionless parameters                                  %
%   Outputs:                                                                 %
%           R_KE     -> Kite-Earth Rotation matrices                         % 
%           rK       -> Kite position vectors (SE components)                % 
%           vK       -> Kite velocity vectors (SB components)                % 
%           eulerK   -> Euler angles                                         % 
%           omegaK   -> Kite Angular velocities (SB components)              % 
%           F_Tether -> Tether Forces upon the kites                         % 
%           M_Tether -> Tether torque upon the kites                         % 
%           FA       -> Aerodynamic force upon the kites                     % 
%           MA       -> Aerodynamic force upon the kites                     % 
%           alfa     -> Angles of attack of the kites                        % 
%           beta     -> Sideslip angles of the kites                         % 
%           rM       -> Position vectors of the masses (SE components)       % 
%           DF       -> Time derivative of the extended state vector         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do calculations for each kite      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rK     = zeros(3,PND.Kite.Num);   % Position vector (SE components)
vK     = zeros(3,PND.Kite.Num);   % Velocity vector (SB components)
eulerK = zeros(3,PND.Kite.Num);   % Euler angles 
omegaK = zeros(3,PND.Kite.Num);   % Angular velocity (SB components)
R_KE   = zeros(3,3,PND.Kite.Num); % Rotation Matrices
vw     = zeros(3,PND.Kite.Num);   % Wind velocity
VA     = zeros(3,PND.Kite.Num);   % Aerodynamic velocity
FA     = zeros(3,PND.Kite.Num);   % Aerodynamic force
MA     = zeros(3,PND.Kite.Num);   % Aerodynamic torque
alfa   = zeros(PND.Kite.Num,1);   % Angle of attack
beta   = zeros(PND.Kite.Num,1);   % Sideslip angle

p      = zeros(PND.Kite.Num,1);
q      = zeros(PND.Kite.Num,1);
r      = zeros(PND.Kite.Num,1);
Phi    = zeros(PND.Kite.Num,1);
Theta  = zeros(PND.Kite.Num,1);
Psi    = zeros(PND.Kite.Num,1);

xc  = Fun_Control_KE(t,PND);

for i=1:1:PND.Kite.Num
    % Kinemtic quantities
    rK(:,i)     = X(12*(i-1)+1:12*(i-1)+3,1);   % Position vector (SE components)
    vK(:,i)     = X(12*(i-1)+4:12*(i-1)+6,1);   % Velocity vector (SB components)
    eulerK(:,i) = X(12*(i-1)+7:12*(i-1)+9,1);   % Euler angles 
    omegaK(:,i) = X(12*(i-1)+10:12*(i-1)+12,1); % Angular velocity (SB components)
    % Euler Angles
    Phi(i)      = eulerK(1,i);
    Theta(i)    = eulerK(2,i);
    Psi(i)      = eulerK(3,i);
    % Angular velocity (SB components)
    p(i)        = omegaK(1,i); 
    q(i)        = omegaK(2,i); 
    r(i)        = omegaK(3,i); 
    % Rotation matrix
    R_KE(:,:,i) = Matrix_R_KE_KE(eulerK(:,i));
    
    % Aerodynamic velocity (SK components)
    vw     = Fun_Wind(t,rK(:,i),PND); % Wind velocity (SE components)
    VA     = vK(:,i)-squeeze(R_KE(:,:,i))*vw;
    
    % Compute Aerodynamic Forces
    delta_a = xc(3*(i-1)+1,1);  
    delta_r = xc(3*(i-1)+2,1);
    delta_e = xc(3*(i-1)+3,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% Components in the body frame of the aerodynamic force and torque 
   
    PND0.Aero.vt   = PND.Aero.vt(i);           % Vref/sqrt(g*L0)  
    PND0.Aero.Full = PND.Aero.Full(i);         % Aerodynamic Model
    
    PND0.Kite.mu     = PND.Kite.mu(i); 
    PND0.Kite.c      = PND.Kite.c(i);      % eps_c
    PND0.Kite.b      = PND.Kite.b(i);      % eps_b

    % Recover the parameters
    if PND.Aero.Full(i) ==  1
        PND0.Aero.CX   = PND.Aero.CX(i);
        PND0.Aero.CY   = PND.Aero.CY(i);
        PND0.Aero.CZ   = PND.Aero.CZ(i);
    
        PND0.Aero.Cm   = PND.Aero.Cm(i);
        PND0.Aero.Cl   = PND.Aero.Cl(i);
        PND0.Aero.Cn   = PND.Aero.Cn(i);
    else
        PND0.Aero.CX   = PND.Aero.CX(:,:,i);
        PND0.Aero.CY   = PND.Aero.CY(:,:,i);
        PND0.Aero.CZ   = PND.Aero.CZ(:,:,i);
    
        PND0.Aero.Cm   = PND.Aero.Cm(:,:,i);
        PND0.Aero.Cl   = PND.Aero.Cl(:,:,i);
        PND0.Aero.Cn   = PND.Aero.Cn(:,:,i); 
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%
    [FA(:,i) MA(:,i) alfa(i) beta(i)] = Aerokite(VA,omegaK(:,i),PND0,delta_a,delta_r,delta_e);
    
    % Gravitational acceleration
    aW(:,i)    = squeeze(R_KE(:,:,i))*[0 0 1]';
    
    % Normalized angular momentum
    L_G(:,i) = PND.Kite.iota(:,:,i)*omegaK(:,i);
end
% Recover Mass position (rM) and velocity (vM) and compute the weight and aerodynamic force (F_Mass)
Ind  = 12*PND.Kite.Num;
Ind0 = 12*PND.Kite.Num + 3*PND.Tether.MassT; 

for i=1:1:PND.Tether.Num  
    for j=1:1:PND.Tether.Mass(i)
        rM(:,i,j)   = X(Ind+1:Ind+3,1); % Masses position vectors
        Ind         = Ind+3;
        vM(:,i,j)   = X(Ind0+1:Ind0+3,1); % Masses velocity vectors
        Ind0        = Ind0+3;
        
        vw          = Fun_Wind(t,rM(:,i,j),PND);         % Wind velocity (SE components)
        VA          = vM(:,i,j)-vw;                      % Aerodynamic velocity of the mass (SE components)
        
        F_Mass(:,i,j) = PND.Mass.Sigma(i,j)*[0 0 1]' - PND.Mass.Xi(i,j)*sqrt(VA'*VA)*VA;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do calculations for each tether %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_Tether = zeros(3,PND.Kite.Num);  % Total elastic Force upon kites (SB components)
M_Tether = zeros(3,PND.Kite.Num);  % Total elastic Torque upon kites (SB components)

T_U = zeros(3,PND.Kite.Num);  % Total elastic Force at "up"  kite's  points (SB components)
T_D = zeros(3,PND.Kite.Num);  % Total elastic Force at "down" kite's points (SB components)

for i=1:1:PND.Tether.Num  % Do all the tethers
    
    %seda = (PND.Tether.Mass(i)+1)*PND.Tether.seda(i); 
    
    seda   =  xc(3*PND.Kite.Num + 2*(i-1)+1);
    seda_p =  xc(3*PND.Kite.Num + 2*(i-1)+2);
   
    
    for j=1:1:PND.Tether.Mass(i)+1 % Number of springs for  tether "i"
        Index_Up   = [0 0 0];     % Kite/Mass, index, index
        Index_Down = [0 0 0];
        
        if j==1 %The lower point of this springs is a kite or the ground  
            %Position vector of the Lower point  (SE components)
            Kite_Down = PND.Tether.Down(i);
            if Kite_Down==0
                rD = [0 0 0]';
                vD = [0 0 0]';
            else
               GD         = [PND.Tether.Dx(i)  PND.Tether.Dy(i)   PND.Tether.Dz(i)]';    % SB components
               rD         = rK(:,Kite_Down) + squeeze(R_KE(:,:,Kite_Down))'*GD;          % SE components
               vD         = squeeze(R_KE(:,:,Kite_Down))'*(vK(:,Kite_Down) + cross(omegaK(:,Kite_Down),GD)); % SE components
               Index_Down = [1 Kite_Down 0];
            end
            if PND.Tether.Mass(i)==0
                %Position vector of the Upper point  (SE components)
                Kite_Up    = PND.Tether.Up(i); 
                GU         = [PND.Tether.Ux(i)  PND.Tether.Uy(i)   PND.Tether.Uz(i)]';    % SB components
                rU         = rK(:,Kite_Up) + squeeze(R_KE(:,:,Kite_Up))'*GU;              % SE components
                vU         = squeeze(R_KE(:,:,Kite_Up))'*(vK(:,Kite_Up) + cross(omegaK(:,Kite_Up),GU)); % SE components
                Index_Up   = [1 Kite_Up 0];
            else
                rU        = rM(:,i,j);    % Masses position vectors (SE components)
                vU        = vM(:,i,j);    % Mass velocity vector (SE components)
                Index_Up   = [2 i j];
            end
        else
            rD         = rM(:,i,j-1);          % Masses position vectors (SE components)
            vD         = vM(:,i,j-1);          % Mass velocity vector (SE components)
            Index_Down = [2 i j-1];
            if j == PND.Tether.Mass(i)+1
                Kite_Up    = PND.Tether.Up(i); 
                GU         = [PND.Tether.Ux(i)  PND.Tether.Uy(i)   PND.Tether.Uz(i)]';    % SB components
                rU         = rK(:,Kite_Up) + squeeze(R_KE(:,:,Kite_Up))'*GU;              % SE components
                vU         = squeeze(R_KE(:,:,Kite_Up))'*(vK(:,Kite_Up) + cross(omegaK(:,Kite_Up),GU)); % SE component
                Index_Up   = [1 Kite_Up 0];
            else
                rU         = rM(:,i,j);          % Masses position vectors (SE components)
                vU         = vM(:,i,j);          % Mass velocity vector (SE components)    
                Index_Up   = [2 i j];
            end
        end
        % Compute elongation
        DU        = rU-rD;                 % Vector with origin at the Down-Point and tip at the Up-Point
        Elong     = seda*sqrt(DU'*DU)-1;
        % Time Derivative of the elongation 
        Vrel      = vU-vD; % Earth components
        Elong_t   = seda*(DU'*Vrel)/sqrt(DU'*DU)+ seda_p*sqrt(DU'*DU);
        % Assign the forces to the kites and the masses 
        if Elong >0
           Mod                   = PND.Tether.nu*(Elong + PND.Tether.ft*Elong_t);
           if Index_Down(1)==1 %Assign force and moment to kite
                F0_Down               = squeeze(R_KE(:,:,Kite_Down))*Mod*DU/sqrt(DU'*DU);   % SB components
                F_Tether(:,Kite_Down) = F_Tether(:,Kite_Down) + F0_Down;                % SB components
                M_Tether(:,Kite_Down) = M_Tether(:,Kite_Down) + cross(GD,F0_Down);      % SB components        
                % Make copy for outputs
                T_D(:,Kite_Down)      = T_D(:,Kite_Down) + F0_Down;     % SB components of tension at an "up point"
           end
           if Index_Up(1)==1 %Assign force and moment to kite
                F0_Up                     =-squeeze(R_KE(:,:,Kite_Up))*Mod*DU/sqrt(DU'*DU); % SB components
                F_Tether(:,Kite_Up)       = F_Tether(:,Kite_Up) + F0_Up;                  % SB components
                M_Tether(:,Kite_Up)       = M_Tether(:,Kite_Up) + cross(GU,F0_Up);        % SB components
                % Make  copy for outputs
                T_U(:,Kite_Up)       = T_U(:,Kite_Up) + F0_Up;     % SB components of tension at an "up point"
           end
           if Index_Down(1)==2 %Assign force to mass
                F0_Down               = Mod*DU/sqrt(DU'*DU);   % SE components
                F_Mass(:,Index_Down(2),Index_Down(3)) = F_Mass(:,Index_Down(2),Index_Down(3))+F0_Down; % SE components
           end
           if Index_Up(1)==2 %Assign force to mass
                F0_Up                 = -Mod*DU/sqrt(DU'*DU);   % SE components
                F_Mass(:,Index_Up(2),Index_Up(3)) = F_Mass(:,Index_Up(2),Index_Up(3))+F0_Up; % SE components
           end
         
        end % Close Elong>0
    end  % Close spring loop 
end % Close tether Loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%      Compute Right hand side       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DF = zeros(12*PND.Kite.Num+3*PND.Tether.MassT,1);
% Aircraft
for i=1:1:PND.Kite.Num
    % Total acceleration and total torque  (SK components)
    aT = aW(:,i) + (FA(:,i) + F_Tether(:,i))/PND.Kite.Sigma(i); 
    MT =            MA(:,i) + M_Tether(:,i);

    % Kinematics relations 
    DF(12*(i-1)+1:12*(i-1)+3,1) = squeeze(R_KE(:,:,i))'*vK(:,i);   
    % F = ma
    DF(12*(i-1)+4:12*(i-1)+6,1) = aT - cross(omegaK(:,i),vK(:,i));
    % Kinematics relations
    DF(12*(i-1)+7,1)   = p(i) + (q(i)*sin(Phi(i))+r(i)*cos(Phi(i)))*tan(Theta(i));
    DF(12*(i-1)+8,1)   = q(i)*cos(Phi(i)) - r(i)*sin(Phi(i));
    DF(12*(i-1)+9,1)   = (q(i)*sin(Phi(i))+r(i)*cos(Phi(i)))/cos(Theta(i));
    % Angular momentum
    DF(12*(i-1)+10:12*(i-1)+12,1) = PND.Kite.iota(:,:,i)\(MT-cross(omegaK(:,i),L_G(:,i)));
end

% Tether masses
Ind  = 12*PND.Kite.Num;
Ind0 = 12*PND.Kite.Num + 3*PND.Tether.MassT;
for i=1:1:PND.Tether.Num  % Do all the tethers
    for j=1:1:PND.Tether.Mass(i) % Mass loop 
        DF(Ind+1:Ind+3,1)   = vM(:,i,j);
        DF(Ind0+1:Ind0+3,1) = F_Mass(:,i,j)/PND.Mass.Sigma(i,j);
        Ind                 = Ind+3;
        Ind0                = Ind0+3;
    end
end
