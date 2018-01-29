function DF = Fun_ODE_Lag_KF(t,xs_amp)

global PND
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Right-Han-Side of Lagrangian equations                         %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%              t                         -> normalized time               %%
%              Extended state vector     -> xs_amp = [xs xs_p]            %%
% Outputs      Time derivative of xs_amp -> DF = d xs_amp/ dt             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xc_amp  = Fun_Control_KF(t,xs_amp,PND);


% Equations of motion of the AWE System
[rQ rR_Edge R_KE rR vR omegaR FA_R QR rK vK omegaK FA_K MA_K MMC_K alfa_K beta_K QK rG vG omegaG FA_G MA_G MMC_G QG DF RHS]=  Fun_ODE_Lag_Full_Output_KF(t,xs_amp,xc_amp,PND);

% Close Loop-Control
if PND.Control.Type == 5 % Use close loop
    
    Nvar     = 2*PND.Num.N+3;              % Number of variables 
    Nvar_p   = 2*PND.Num.N+3+PND.Gen.Num;  % Number of variables (dot) 
  
    Euler    = xs_amp(2*PND.Num.N+1:2*PND.Num.N+3);    % Euler Angles
    Euler_p  = xs_amp(4*PND.Num.N+4:4*PND.Num.N+6);    % First Derivative of Euler Angles
    Euler_pp =     DF(4*PND.Num.N+4:4*PND.Num.N+6);    % Second Derivative of Euler Angles
     
    DF(Nvar+Nvar_p+1:Nvar+Nvar_p+3) = Fun_PID(t,xs_amp,Euler,Euler_p,Euler_pp,PND);
   
end
    
end