function [X_Amp Error Flag] = Equilibrium_KE(X0,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute equilibrium state                                      %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:                                                                  %
%           X0  -> Intial guess (set to zero if an initial guess is unknown) %
%           PND -> Dimensionless Parameters                                  %
%   Outputs:                                                                 %
%           X_Amp -> State vector at the equilibrium                         % 
%           Error -> Equilibrium Error                                       %
%           Flag  -> 1 (successful) or 0 (calculation failed)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate an initial guess for the Newton Method
if X0==0
    theta = 15*pi/180;
    X0    = zeros(3*PND.Kite.Num,1);
    for i=1:1:PND.Kite.Num % Kites
        seda = 1e9;
        for j=1:1:PND.Tether.Num % Select the tether with maximum length
            if PND.Tether.Down(j)==i-1 && PND.Tether.Up(j)==i
                if PND.Tether.seda0(j)<seda
                   Length = PND.Tether.seda0(j)/(PND.Tether.Mass(i)+1); 
                end
            end
        end
        xg    = -1.08/Length*cos(75*pi/180);
        zg    = -1.08/Length*sin(75*pi/180);
        if i==1
           X0(3*(i-1)+1) = xg;
           X0(3*(i-1)+2) = zg;
        else
           X0(3*(i-1)+1) = X0(3*(i-2)+1) + xg;
           X0(3*(i-1)+2) = X0(3*(i-2)+2) + zg; 
        end
        X0(3*(i-1)+3) = theta;
        rK(:,i)       = [X0(3*(i-1)+1) 0 X0(3*(i-1)+2)]';
    end
    Ind = 3*PND.Kite.Num;
    for i=1:1:PND.Tether.Num
        R_KE  = Matrix_R_KE_KE([0 theta 0]'); 
        % Low Point
        if PND.Tether.Down(i)==0
            rD = [0 0 0]';       
        else
            GD    = [PND.Tether.Dx(i)  PND.Tether.Dy(i)   PND.Tether.Dz(i)]';          % SB components
            rD    = rK(:,PND.Tether.Down(i)) + R_KE'*GD;                               % SE components
        end
        % Up Point       
        GU    = [PND.Tether.Ux(i)  PND.Tether.Uy(i)   PND.Tether.Uz(i)]';          % SB components
        rU    = rK(:,PND.Tether.Up(i)) +   R_KE'*GU;                               % SE components
         
        for j=1:1:PND.Tether.Mass(i) % Masses
            X0(Ind+1:Ind+3) = rD+j*(rU-rD)/(PND.Tether.Mass(i)+1);
            Ind             = Ind+3;
        end
    end
end
clear Ind rU GU rD GD R_KE rK xg zg seda
% Call Newton-Method
[X Error Flag] = my_fzero(@Equilibrium,X0,PND); 
% Recover the full State-Vector
for i=1:1:PND.Kite.Num
   x                        = X(3*(i-1)+1,1);
   z                        = X(3*(i-1)+2,1);
   Theta                    = X(3*(i-1)+3,1);
   rK(:,i)                  = [x 0 z]';     % Position vector (SE components)
   v                        = zeros(3,1);   % Velocity vector (SB components)
   eulerK                   = [0 Theta 0]'; % Euler angles 
   omega                    = zeros(3,1);   % Angular velocity (SB components)
   R_KE(:,:,i)              = Matrix_R_KE_KE(eulerK);
   X_Amp(12*(i-1)+1:12*i,1) = [rK(:,i);v;eulerK;omega];   
end
% Add masses positions
X_Amp = [X_Amp;X(3*PND.Kite.Num+1:3*PND.Kite.Num+3*PND.Tether.MassT)];
% Add masses velocities
X_Amp = [X_Amp;zeros(3*PND.Tether.MassT,1)];

     function DF = Equilibrium(t,X_Eq)

            % Kites variables
            
            for i=1:1:PND.Kite.Num
                  x                        = X_Eq(3*(i-1)+1,1);
                  z                        = X_Eq(3*(i-1)+2,1);
                  Theta                    = X_Eq(3*(i-1)+3,1); 
                  r                        = [x 0 z]';     % Position vector (SE components)
                  v                        = zeros(3,1);   % Velocity vector (SB components)
                  euler                    = [0 Theta 0]'; % Euler angles 
                  omega                    = zeros(3,1);   % Angular velocity (SB components)
                  X_Amp(12*(i-1)+1:12*i,1) = [r;v;euler;omega];
            end
            % Discrete masses positions 
            X_Amp = [X_Amp;X_Eq(3*PND.Kite.Num+1:3*PND.Kite.Num+3*PND.Tether.MassT)];
            % Discrete masses velocities
            X_Amp = [X_Amp;zeros(3*PND.Tether.MassT,1)];        
            % Call de full RHS
            DF0   = Fun_ODE_KE(0,X_Amp);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Compute the reduced RHS %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kites
            for i=1:1:PND.Kite.Num
                  DF(3*(i-1)+1,1) = DF0(12*(i-1)+4,1);  % dvx/dt = 0
                  DF(3*(i-1)+2,1) = DF0(12*(i-1)+6,1);  % dvz/dt = 0
                  DF(3*(i-1)+3,1) = DF0(12*(i-1)+11,1); % dq/dt  = 0
            end
            % Masses
            DF(3*PND.Kite.Num+1:3*(PND.Kite.Num+PND.Tether.MassT),1) = DF0(12*PND.Kite.Num+3*PND.Tether.MassT+1:12*PND.Kite.Num+6*PND.Tether.MassT);
            clear X_Amp
     end        
end