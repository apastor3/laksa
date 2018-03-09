

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Example                                                        %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %%
% Inputs: No inputs                                             %%
%                                                               %%
% Outputs: Integration Results are placed in the workspace      %% 
%                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Define the dimensionless parameters as global variables
global PND

% Add the path of the common folder
addpath('../Common/')

% Load the physical parameters of KiteAcrobat
PD  = Fun_PD_KA;

% Find the dimensionless parameters of KiteAcrobat
PND = Fun_PND_KA(PD);

% Compute equilibrium state for PND.Env.Type=0 & PND.Ctr.Type=0
[u0 Error Flag]=Equilibrium_KA(0,PND);

% Integrate the Equations of motion
options = odeset('RelTol',PND.Num.RelTol,'AbsTol',PND.Num.AbsTol);

[T u] = ode45('Fun_ODE_Lag_KA',[0:0.01:.4]*pi/PND.Ctr.delta1,u0,options);


Flag_Dim = 1;     % --> 1 Get the results with dimensions
for i=1:1:length(T)
   [Tout(i) RBE(:,:,i) rk(:,i)  vk(:,i)  ak(:,i) euler(:,i) omega(:,i) omega_p(:,i)...
   Lambda(:,i) FAP(:,i) FAM(:,i)   MAP(:,i) MAM(:,i) FA(:,i) MA(:,i)    W(:,i)...
   alfa(i)  beta(i)     LP(i)    LM(i)]                               = Fun_Post_KA(PD,PND,T(i),u(i,:)',Flag_Dim);
   Plot_KA(T(i),rk(:,i),RBE(:,:,i),PND,Flag_Dim,PD)
   pause(0.001)
   title('')
end
Flag_Plot = [1 1 1 1 1 1 ]; %[Position Velocity Euler alfa&beta Tension Lengths ]
Plot_Results_KA(Tout, RBE, rk, vk, ak, euler, omega, omega_p, Lambda, FAP, FAM, MAP, MAM, FA, MA, W, alfa, beta, LP, LM,Flag_Dim,Flag_Plot)


%% Repeat the calculations with the Hamiltonian formulation 
w0    = From_Lag2Ham_KA(0,u0,PND);
[T w] = ode45('Fun_ODE_Ham_KA',T,w0,options);

clear u Tout RBE rk  vk  ak euler omega omega_p Lambda FAP FAM MAP MAM FA MA W alfa beta LP LM 
for i=1:1:length(T)
   u = From_Ham2Lag_KA(T(i),w(i,:)',PND);
   [Tout(i) RBE(:,:,i) rk(:,i)  vk(:,i)  ak(:,i) euler(:,i) omega(:,i) omega_p(:,i)...
   Lambda(:,i) FAP(:,i) FAM(:,i)   MAP(:,i) MAM(:,i) FA(:,i) MA(:,i)    W(:,i)...
   alfa(i)  beta(i)     LP(i)    LM(i)]                               = Fun_Post_KA(PD,PND,T(i),u,Flag_Dim);
end
Flag_Plot = [1 1 1 1 1 1 ]; %[Position Velocity Euler alfa&beta Tension Lengths ]
Plot_Results_KA(Tout, RBE, rk, vk, ak, euler, omega, omega_p, Lambda, FAP, FAM, MAP, MAM, FA, MA, W, alfa, beta, LP, LM,Flag_Dim,Flag_Plot)



   