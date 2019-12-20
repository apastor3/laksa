

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Example                                                        %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %%
% Inputs: No inputs                                             %%
%                                                               %%
% Outputs: Integration Results are placed in the workspace      %% 
%                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

global PND

% Add the path of the common folder
addpath('../Common/')

% Parameter with and without dimensions
PD = Fun_PD_KE;
PND = Fun_PND_KE(PD);

Flag_Dim = 1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% Compute Equilibrium  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('---------------------------------------------------------')
display('Computing Equilibrium State')
 
X0 = 0;
T  = 0;
[X_Amp Error Flag] = Equilibrium_KE(X0,PND);
% Get the data with dimensions
Flag_Plot =  [0            0         0       0             0           ];
[R_KE rK vK eulerK omegaK F_Tether M_Tether FA MA alfa beta  rM DF T_U T_D] = Fun_ODE_Full_Output_KE(T,X_Amp,PND);
[T R_KE rk vK eulerK omegaK F_Tether  FA  alfa  beta rM T_U T_D]            = Plot_Results_KE(T, R_KE, rK, vK, eulerK, omegaK, F_Tether, FA, alfa, beta,rM,T_U,T_D,Flag_Plot,Flag_Dim,PD,PND);

% Stability Analysis
display('---------------------------------------------------------')
display('Stability Analysis of the Equilibrium State')
Jnum        = Jacobian('Fun_ODE_KE',0,X_Amp,PND);
[Vec Val]   = eig(Jnum);

figure(73)
plot(real(Val),imag(Val),'+')
xlabel('Real(\lambda)','fontsize',12)
ylabel('Imag(\lambda)')
set(gca,'box','on','fontsize',12)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute tension modulus
for i=1:1:PD.Kite.Num 
      Tension(i) = sqrt(squeeze(T_U(:,i))'*squeeze(T_U(:,i)))/2;
end

figure(5)

subplot(3,1,1)
hold on
plot(-rk(1,:),-rk(3,:),'-b')
plot(-rk(1,:),-rk(3,:),'b+')
text(0,500,'(a)','fontsize',12)
xlabel('$-x_E\ (m)$','fontsize',12,'interpreter','latex')
ylabel('$-z_E\ (m)$','fontsize',12,'rotation',0,'interpreter','latex','fontsize',14)
set(gca,'fontsize',12,'box','on')

subplot(3,1,2)
hold on
plot([1:1:PD.Kite.Num  ],alfa,'-b')
plot([1:1:PD.Kite.Num  ],alfa,'b+')
text(0,7,'(b)','fontsize',12)
xlabel('Kite number','fontsize',12)
ylabel('$\alpha_i\ (^o)$','fontsize',12,'rotation',0,'interpreter','latex','fontsize',14)
set(gca,'fontsize',12,'box','on')


subplot(3,1,3)
hold on
plot([1:1:PD.Kite.Num  ],Tension,'-b')
plot([1:1:PD.Kite.Num  ],Tension,'b+')
text(0,0,'(c)','fontsize',12)
xlabel('Kite number','fontsize',12)
ylabel('$|T_U^+|\ (N)$','fontsize',12,'rotation',0,'interpreter','latex','fontsize',14)
set(gca,'fontsize',12,'box','on')
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate 
%%%%%%%%%%%%%%%%%%%%%%%%
display('---------------------------------------------------------')
display('Integrating Equations of motion')
display('Initial Condition = Equilibrium State + Random Perturbation')
    
[T X]          = ode15s('Fun_ODE_KE',[0:0.5:20],X_Amp+0.001*rand(12*PND.Kite.Num+6*PND.Tether.MassT,1));
Ind = 1;
hFig = figure(23);
set(hFig, 'Position', [100 100 1000 600]);
for k=1:1:length(T)

 [R_KE(:,:,:,k) rK(:,:,k) vK(:,:,k) eulerK(:,:,k) omegaK(:,:,k)...
 F_Tether(:,:,k) M_Tether(:,:,k) FA(:,:,k) MA(:,:,k) alfa(:,k) beta(:,k)  rM(:,:,:,k) DF(:,k) T_U(:,:,k) T_D(:,:,k)] = Fun_ODE_Full_Output_KE(T(k),X(k,:)',PND);

 ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
 cla(ax1)
 hold on
 set(get(ax1, 'title'), 'string',['Time = ', num2str(T(k))])
 Plot_KE(T(k),X(k,:)',PND,ax1)
 pause(0.001)
 set(get(ax1, 'title'), 'string','                  ')
 
end
%             Position  Velocity   Euler alfa&beta      Forces      
Flag_Plot =  [1            1         1       1             1           ];
                                             
[T R_KE rK vK eulerK omegaK F_Tether  FA  alfa  beta rM T_U T_D] = Plot_Results_KE(T, R_KE, rK, vK, eulerK, omegaK, F_Tether, FA, alfa, beta,rM,T_U,T_D,Flag_Plot,Flag_Dim,PD,PND);
 

