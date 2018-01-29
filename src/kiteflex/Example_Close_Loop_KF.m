%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Example for a ground-generation system                         %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %%
% Inputs: No inputs                                                       %%
%                                                                         %%
% Outputs: Integration Results are placed in the workspace                %% 
%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

global PND

% Add the path of the common folder
addpath('../common/')

Flag_Dim = 1;
% Construct the dimensionless parameters
PD = Fun_PD_Close_Loop_KF;
% Construct the dimensionless parameters
PND             = Fun_PND_KF(PD);


% Compute Equilibrium
[u0  Error Flag PND delta_a_Eq xi_Rotor_Eq]=Equilibrium_FlyGen_KF(0,PND);
PND.Target.delta_a   = delta_a_Eq;
PND.Target.xi_Rotor  = xi_Rotor_Eq;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Open Loop Configuration                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Open Loop Results')
% Set control to the target variables (elevator deflection and motor controller torque)
PND.Control.Type     = 2;  
% Check that the right-hand-side vanishes (u0 is equilibrium state)   
Error                = max(abs(Fun_ODE_Lag_KF(0,u0)))
% Compute Jacobian
Jnum                 = Jacobian('Fun_ODE_Lag_KF',0,u0,PND);
% Compute eigenvalues
Val                  = eig(Jnum)


Pert = [0.231594386708524    0.488897743920167    0.624060088173690    0.679135540865748    0.395515215668593...
        0.367436648544477    0.987982003161633    0.037738866239552    0.885168008202475    0.913286827639239...
        0.796183873585212    0.098712278655574    0.261871183870716    0.335356839962797    0.679727951377338...
        0.136553137355370    0.721227498581740    0.106761861607241    0.653757348668560    0.494173936639270...
         ]';
if 1==2
% Integrate the equations of motion
[T X] = ode45('Fun_ODE_Lag_KF',[0:0.01:25.0],u0+0.001*Pert);
   

for i=1:1:length(T)
       [Tout(i) rR_Edge(:,:,i) rR(:,:,i) vR(:,:,i) aR(:,:,i) omegaR(:,:,i) gR(:,:,i) FA_R(:,:,i) Tension(:,:,i) rQ(:,i) R_KE(:,:,i) ...
       rK(:,i) vK(:,i) aK(:,i) euler(:,i) omegaK(:,i) gK(:,i) FA_K(:,i) FR_K(:,i) FG_K(:,i)  MA_K(:,i) MR_K(:,i) MG_K(:,i) MMC_K(:,i) alfa_K(i) beta_K(i)... 
       rG(:,:,i) vG(:,:,i) aG(:,:,i) omegaG(:,:,i) gG(:,:,i) FA_G(:,:,i)  FK_G(:,:,i) MA_G(:,:,i) MK_G(:,:,i)  MMC_G(:,:,i) xc_out(:,i) xs_target_out(:,i) Error_C(i)] = Fun_Post_KF(PD,PND,T(i),X(i,:)',Flag_Dim);
          
       %Plot_FlyGen_KF(X(i,:)',Tout(i),rQ(:,i),rR(:,:,i),rK(:,i),rG(:,:,i),R_KE(:,:,i),rR_Edge(:,:,i),PND,Flag_Dim,PD);

       pause(0.01)
       title('')
end
            %[Kite Position, Velocity, Euler, alfa&beta, Tension, Control, Rotor angular velocity, Error_C
Flag_Plot = [1             ,     1   ,   1  ,     1    ,    1   ,    1   ,       1               ,   1  ];
     
Plot_Results_KF(Tout,rR_Edge, rR, vR, aR, omegaR, gR, FA_R, Tension, rQ, R_KE, ...
    rK, vK, aK, euler, omegaK, gK, FA_K, FR_K, FG_K,  MA_K, MR_K, MG_K, MMC_K, alfa_K, beta_K,... 
    rG, vG, aG, omegaG, gG, FA_G,  FK_G, MA_G, MK_G,  MMC_G, xc_out,xs_target_out,Error_C, Flag_Dim, Flag_Plot,PD)

end
'--------------------------------------------------------'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Close Loop Configuration                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Close Loop Results')
% Set control close loop control
PND.Control.Type      = 5;  
% Add the deflection of the control surface to the state vector
u0                    = [u0;delta_a_Eq;0;0];  %[u0;Aileron;Rudder;Elevator]
% Check that the right-hand-side vanishes (u0 is equilibrium state)   
Error                = max(abs(Fun_ODE_Lag_KF(0,u0)))
% Compute Jacobian
Jnum                 = Jacobian('Fun_ODE_Lag_KF',0,u0,PND);
% Compute eigenvalues
Val                  = eig(Jnum)

% Random Perturbation
Pert = [0.231594386708524    0.488897743920167    0.624060088173690    0.679135540865748    0.395515215668593...
        0.367436648544477    0.987982003161633    0.037738866239552    0.885168008202475    0.913286827639239...
        0.796183873585212    0.098712278655574    0.261871183870716    0.335356839962797    0.679727951377338...
        0.136553137355370    0.721227498581740    0.106761861607241    0.653757348668560    0.494173936639270...
        0.779051723231275    0.715037078400694    0.903720560556316]';

% Integrate the equations of motion
[T X] = ode45('Fun_ODE_Lag_KF',[0:0.01:25.0],u0+0.001*Pert);
   

for i=1:1:length(T)
       [Tout(i) rR_Edge(:,:,i) rR(:,:,i) vR(:,:,i) aR(:,:,i) omegaR(:,:,i) gR(:,:,i) FA_R(:,:,i) Tension(:,:,i) rQ(:,i) R_KE(:,:,i) ...
       rK(:,i) vK(:,i) aK(:,i) euler(:,i) omegaK(:,i) gK(:,i) FA_K(:,i) FR_K(:,i) FG_K(:,i)  MA_K(:,i) MR_K(:,i) MG_K(:,i) MMC_K(:,i) alfa_K(i) beta_K(i)... 
       rG(:,:,i) vG(:,:,i) aG(:,:,i) omegaG(:,:,i) gG(:,:,i) FA_G(:,:,i)  FK_G(:,:,i) MA_G(:,:,i) MK_G(:,:,i)  MMC_G(:,:,i) xc_out(:,i) xs_target_out(:,i) Error_C(i)] = Fun_Post_KF(PD,PND,T(i),X(i,:)',Flag_Dim);
          
      % Plot_FlyGen_KF(X(i,:)',Tout(i),rQ(:,i),rR(:,:,i),rK(:,i),rG(:,:,i),R_KE(:,:,i),rR_Edge(:,:,i),PND,Flag_Dim,PD);

       pause(0.01)
       title('')
end
            %[Kite Position, Velocity, Euler, alfa&beta, Tension, Control, Rotor angular velocity, Error_C
Flag_Plot = [1             ,     1   ,   1  ,     1    ,    1   ,    1   ,       1               ,   1  ];
     
Plot_Results_KF(Tout,rR_Edge, rR, vR, aR, omegaR, gR, FA_R, Tension, rQ, R_KE, ...
    rK, vK, aK, euler, omegaK, gK, FA_K, FR_K, FG_K,  MA_K, MR_K, MG_K, MMC_K, alfa_K, beta_K,... 
    rG, vG, aG, omegaG, gG, FA_G,  FK_G, MA_G, MK_G,  MMC_G, xc_out,xs_target_out,Error_C, Flag_Dim, Flag_Plot,PD)

