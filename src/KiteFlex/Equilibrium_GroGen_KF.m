function [u0  Error Flag PND ]=Equilibrium_GroGen_KF(u0,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute Equilibrium in a Ground-Generation System              %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%----------------------------------------------------------------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %%
% Inputs:  u0       -> Initial guess                                      %%
%          PND      -> Dimensionless parameters                           %%
%                                                                         %%
% Outputs: u0       -> Vector that makes | function(0,u0) <Error          %%
%          Error    -> Error of the solution                              %%
%          Flag     -> 1-> Success, 0 -> The method did not converge      %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save original  wind type
Wind_Type    = PND.Env.Type;

% Set steady wind and target control%
PND.Env.Type       = 0;

if u0 == 0 % The user did not insert any initial guess 
    % State variables
    u0 =  pi/4*ones(PND.Num.N,1);
    u0 = [u0;zeros(PND.Num.N,1)];
    u0 = [u0;15*pi/180;0;0];
    % Time derivative of xs 
    u0 = [u0;zeros(2*PND.Num.N+3,1)];        
end

% Vector with longitudinal variables
x0_Red              = [u0(1:PND.Num.N); u0(2*PND.Num.N+1)]; 
% Look for equilibrium
[x0_Red Error Flag] = my_fzero(@Fun_ODE_KF_EQ_GroGen,x0_Red,PND);

if Flag ==1    
    u0  = [x0_Red(1:PND.Num.N); zeros(PND.Num.N,1); x0_Red(PND.Num.N+1); 0; 0];
    u0  = [u0;zeros(2*PND.Num.N+3,1)];
    xc_amp   = Fun_Control_KF(0,u0,PND);
       
    [rQ rR_Edge R_KE rR vR omegaR FA_R QR rK vK omegaK FA_K MA_K MMC_K alfa_K beta_K QK rG vG omegaG FA_G MA_G MMC_G QG DF RHS]=  Fun_ODE_Lag_Full_Output_KF(0,u0,xc_amp,PND); 
    Error = max(sqrt(RHS'*RHS))
    display(['Equilibrium Error = ' num2str(Error)])
else
    u0     = 0;
    Error  = 1e9;
    Flag   = 0
end

% Set original variables
PND.Env.Type       = Wind_Type;

    function DF0 = Fun_ODE_KF_EQ_GroGen(t,X)
        
        w00  = [X(1:PND.Num.N); zeros(PND.Num.N,1); X(PND.Num.N+1); 0; 0];
        w00  = [w00;zeros(2*PND.Num.N+3,1)];   
        xc_amp   = Fun_Control_KF(t,w00,PND);
        [rQ rR_Edge R_KE rR vR omegaR FA_R QR rK vK omegaK FA_K MA_K MMC_K alfa_K beta_K QK rG vG omegaG FA_G MA_G MMC_G QG DF RHS]= Fun_ODE_Lag_Full_Output_KF(0,w00,xc_amp,PND);
    
        DF0(1:PND.Num.N,1) = RHS(1:PND.Num.N,1);     % Rods equilibrium
        DF0(PND.Num.N+1,1) = RHS(2*PND.Num.N+1,1);   % Pitch angle equilibrium
    end



end