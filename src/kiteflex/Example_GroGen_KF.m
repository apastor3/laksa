%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Example for a groung-generation system                         %
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
PD       = Fun_PD_GroGen_KF;
% Construct the dimensionless parameters
PND             = Fun_PND_KF(PD);
 
display('Computing Equilibrium State')
[u0  Error Flag PND ]=Equilibrium_GroGen_KF(0,PND);
if Flag==1
    display('Equilibrium state was compputed successfully')
    display('Computing the trajectory')
    [T X] = ode45('Fun_ODE_Lag_KF',[0:0.0025:.1],u0);
    
    for i=1:1:length(T)
       [Tout(i) rR_Edge(:,:,i) rR(:,:,i) vR(:,:,i) aR(:,:,i) omegaR(:,:,i) gR(:,:,i) FA_R(:,:,i) Tension(:,:,i) rQ(:,i) R_KE(:,:,i) ...
       rK(:,i) vK(:,i) aK(:,i) euler(:,i) omegaK(:,i) gK(:,i) FA_K(:,i) FR_K(:,i) FG_K(:,i)  MA_K(:,i) MR_K(:,i) MG_K(:,i) MMC_K(:,i) alfa_K(i) beta_K(i)... 
       rG(:,:,i) vG(:,:,i) aG(:,:,i) omegaG(:,:,i) gG(:,:,i) FA_G(:,:,i)  FK_G(:,:,i) MA_G(:,:,i) MK_G(:,:,i)  MMC_G(:,:,i) xc_out(:,i) xs_target_out(:,i) Error_C(i)] = Fun_Post_KF(PD,PND,T(i),X(i,:)',Flag_Dim);
            
       Plot_GroGen_KF(X(i,:)',Tout(i),rQ(:,i),rR(:,:,i),rK(:,i),rG(:,:,i),R_KE(:,:,i),rR_Edge(:,:,i),PND,Flag_Dim,PD);
       pause(0.01)
       title('')
    end
    %[Kite Position, Velocity, Euler, alfa&beta, Tension, Control, Rotor angular velocity, Error_C
    Flag_Plot = [1             ,     1   ,   1  ,     1    ,    1   ,    1   ,       1               ,   1  ];
              
    Plot_Results_KF(Tout,rR_Edge, rR, vR, aR, omegaR, gR, FA_R, Tension, rQ, R_KE, ...
    rK, vK, aK, euler, omegaK, gK, FA_K, FR_K, FG_K,  MA_K, MR_K, MG_K, MMC_K, alfa_K, beta_K,... 
    rG, vG, aG, omegaG, gG, FA_G,  FK_G, MA_G, MK_G,  MMC_G, xc_out,xs_target_out,Error_C, Flag_Dim, Flag_Plot,PD)

end



   