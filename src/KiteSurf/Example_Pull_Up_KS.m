
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Example program                                                %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%                                                                    %%
% Inputs: No inputs                                                  %%
%                                                                    %%
% Outputs: Integration Results are placed in the workspace           %%
%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Define the dimensionless parameters as global variables
global PND 

% Add the path of the common folder
addpath('../Common/')
% Load the physical parameters of KiteAcrobat
PD          = Fun_PD_KS;
PD.Ctr.Type = 1;
% Find the dimensionless parameters of KiteAcrobat
PND         = Fun_PND_KS(PD);

% Select Dimensions or Dimensionless outputs
Flag_Dim = 1;
         
% Compute Equilibrium
[u0  Error Flag]=Equilibrium_KS(0,PND);

% Integrate the Equations of motion
display('Computing Trajectory')
options = odeset('RelTol',PND.Num.RelTol,'AbsTol',PND.Num.AbsTol);

TF = 15;

[T u]   = ode15s('Fun_ODE_KS',[0:0.001:TF],u0,options);
XF      = u(end,:)';

% Make a Movie and post-process the results
Ind =1;

for i=1:1:length(T)
    
   [T_out(i) RBE(:,:,i) R2E(:,:,i) R3E(:,:,i) rk(:,i) vk(:,i) ak(:,i) euler(:,i) omega(:,i) omega_p(:,i)...
   Lambda(:,i) FAP(:,i) FAM(:,i) MAP(:,i) MAM(:,i) FBP(:,i) FBM(:,i) MBP(:,i) MBM(:,i)...
   FA(:,i) MA(:,i) W(:,i) alfa(i) beta(i) Rp(:,1:3,i)  Rm(:,1:3,i) ...
   Elong_p(:,i) Elong_m(:,i) xc(:,i) Error0(i)] = Fun_Post_KS(PD,T(i),u(i,:)',Flag_Dim,PND);   
 
   if Ind==700 || i==1
       if i>1
          title('')
       end
       Plot_KS(u(i,:)',T(i),rk(:,i),vk(:,i),ak(:,i),RBE(:,:,i),R2E(:,:,i),R3E(:,:,i),Rp(:,:,i),Rm(:,:,i),PND,Flag_Dim,PD)
       Ind = 1;
   else
       Ind = Ind+1;
   end
   pause(0.001)

end

% Plot the results of the simulation
Flag_Plot = [1 1 1 1 1 ...   %[Kite Position, Velocity, Euler, alfa&beta, Tension at inelastic tethers,...
             1 1 1 1 1 ...   % Control, Position of elastic tether, Tension of elastic tethers, Elongation, Moments  ]
             1 1];             % Error Angular Velocity   
           
Plot_Results_KS(T_out, RBE, rk, vk, ak, euler, omega, omega_p, Lambda, FAP, FAM, MAP, MAM, FBP, FBM, MBP, MBM, ...
               FA, MA, W, alfa, beta, Rp, Rm, Elong_p, Elong_m,xc,Error0,Flag_Dim,Flag_Plot)
     
           