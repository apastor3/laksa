function [u0  Error Flag PND delta_a_Eq xi_Rotor_Eq]=Equilibrium_FlyGen_KF(u0,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute equilibrium in a Fly-Generation system                 %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %%
% Inputs:  u0       -> Initial guess                                      %%
%          PND      -> Dimensionless parameters                           %%
%                                                                         %%
% Outputs: u0         -> Vector that makes | function(0,u0) <Error        %%
%          Error      -> Error of the solution                            %%
%          Flag       -> 1-> Success, 0 -> The method did not converge    %% 
%          delta_a_Eq -> 1-> aileron deflection at the equilibrium        %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save original control and wind types
Control_Type = PND.Control.Type  ;
Wind_Type    = PND.Env.Type;
delta_a      = PND.Target.delta_a; 
xi_rotor     = PND.Target.xi_Rotor; 
 
% Set steady wind and target control
PND.Control.Type   = 3;                 % Set Equilibrium control configuration (aileron deflection is computed to compensate the reaction torque )
PND.Env.Type       = 0;

if u0 == 0 % The user did not insert any initial guess 
    % State variables
    u0 =  pi/4*ones(PND.Num.N,1);
    u0 = [u0;zeros(PND.Num.N,1)];
    u0 = [u0;15*pi/180;0;0];
    % Time derivative of xs + angular velocity of the rotors
    u0 = [u0;zeros(2*PND.Num.N+3,1)];    
    u0 = [u0;PND.Target.Ome_Rotor*ones(PND.Gen.Num,1)];
    
    
end

% Extract the unknowns and add the motor controller torque and aileron deflection
PND.Target.xi_Rotor =  PND.Gen.l*PND.Gen_Chi*PND.Gen.Cm*PND.Env.vw^2;
PND.Target.delta_a  =  -4*PND.Target.xi_Rotor/(PND.Kite.mu*PND.Env.vw^2*PND.Kite.b*PND.Aero.Cldelta_a);

x0_Red              = [u0(1:PND.Num.N); u0(2*PND.Num.N+1);PND.Target.xi_Rotor;PND.Target.delta_a]; 
% Look for equilibrium
[x0_Red Error Flag] = my_fzero(@Fun_ODE_KF_EQ_FlyGen,x0_Red,PND);

if Flag ==1  
    PND.Target.xi_Rotor =  x0_Red(PND.Num.N+2);
    PND.Target.delta_a  =  x0_Red(PND.Num.N+3);    
    u0  = [x0_Red(1:PND.Num.N); zeros(PND.Num.N,1); x0_Red(PND.Num.N+1); 0; 0];
    u0  = [u0;zeros(2*PND.Num.N+3,1);PND.Target.Ome_Rotor*ones(PND.Gen.Num,1)];
    xc_amp  = Fun_Control_KF(0,[],PND);
    [rQ rR_Edge R_KE rR vR omegaR FA_R QR rK vK omegaK FA_K MA_K MMC_K alfa_K beta_K QK rG vG omegaG FA_G MA_G MMC_G QG DF RHS]=  Fun_ODE_Lag_Full_Output_KF(0,u0,xc_amp,PND);
    Error = max(sqrt(RHS'*RHS));
    display(['Equilibrium Error = ' num2str(Error)])
else
    u0     = 0;
    Error  = 1e9;
    Flag   = 0
end

% Set original variables
PND.Control.Type    = Control_Type;
PND.Env.Type        = Wind_Type;
delta_a_Eq          = PND.Target.delta_a;
xi_Rotor_Eq         = PND.Target.xi_Rotor;
PND.Target.delta_a  = delta_a;
PND.Target.xi_Rotor = xi_rotor;

    function DF0 = Fun_ODE_KF_EQ_FlyGen(t,X)
        
        PND.Target.xi_Rotor =  X(PND.Num.N+2);
        PND.Target.delta_a  =  X(PND.Num.N+3);

        w00     = [X(1:PND.Num.N); zeros(PND.Num.N,1); X(PND.Num.N+1); 0; 0];
        w00     = [w00;zeros(2*PND.Num.N+3,1);PND.Target.Ome_Rotor*ones(PND.Gen.Num,1)];   
        xc_amp  = Fun_Control_KF(t,[],PND);
        [rQ rR_Edge R_KE rR vR omegaR FA_R QR rK vK omegaK FA_K MA_K MMC_K alfa_K beta_K QK rG vG omegaG FA_G MA_G MMC_G QG DF RHS]=  Fun_ODE_Lag_Full_Output_KF(0,w00,xc_amp,PND);
   
        
        
        DF0(1:PND.Num.N,1) = RHS(1:PND.Num.N,1);     % Rods equilibrium
        DF0(PND.Num.N+1,1) = RHS(2*PND.Num.N+1,1);   % Pitch angle equilibrium
        DF0(PND.Num.N+2,1) = RHS(2*PND.Num.N+3,1);   % Phi angle equilibrium
        if PND.Gen.Num>0
           DF0(PND.Num.N+3,1) = RHS(2*PND.Num.N+3+1,1); % 1 rotor equilibrium
        end
    end



end