%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Example for a groung-generation system (steering manoeuver)    %
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

% Choose output with or without dimensions
Flag_Dim = 1;

% Load physical parameters
PD       = Fun_PD_Steering_KF;
% Construct the dimensionless parameters
PND      = Fun_PND_KF(PD);
% Dimensionless Integration Time
T = 5.4;
     
% Initial conditions    
if PND.Num.N==1       
      u0 = [    0.469233524677695  -0.227397130570123   -0.149390230751017   0.992559120972679   1.180911752930972...
               -0.094335796908624  -0.361418859118427   -3.056073040383586   1.298774216003980   1.334445795576000]';
end
 
if PND.Num.N==2
    u0 = [   0.354734839871189    0.357713412408990  -0.239700153459800  -0.288798563027027...
            -0.310056702802329    0.990750105751635   1.356937981100839...
            -0.073235535834771   -0.232493725808031  -0.344769982936826  -0.310273817634688...
            -3.148743279454925    0.683919300138707   1.332228086703546]'; 
end

if PND.Num.N==3
    u0 = [   0.331931859385302    0.337965204354309    0.335658208906848   -0.241297502398784   -0.259296100879425   -0.304680736035203...
        -0.332805676997182    0.984799380755165    1.388071537429363...
        -0.074331179740453   -0.128664390549729   -0.277211313644053   -0.338353508451515   -0.330161175356686   -0.294834376245505...
        -3.150629219939743    0.572909822873071    1.354685327210377]'; 
end

if PND.Num.N==4
     u0 = [ 0.322854561835376    0.327908149934098    0.330554577840026    0.326161751944103   -0.243221015682775  -0.252560239908856  -0.274248742651091  -0.314187671215925...
      -0.344880939472145    0.982061418911141    1.401543329456379...
      -0.077258033540035   -0.103982959128069   -0.171724233883546   -0.303572607165787   -0.335773414683133  -0.331961823316321  -0.319458489830549  -0.285121425688943...
      -3.149791696490080    0.525150088469110   1.364113752575836]'; 
end
       
  

     
[u0,Error,Floquet,flag] = Corrector_Non_Autonomous('Fun_ODE_Lag_KF',u0,T,5,2e-3,PND.Num.IntTol,PND.Num.IntTol,PND);

figure(1)
hold on
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'b')
plot(real(Floquet),imag(Floquet),'+')


    
[T X] = ode45('Fun_ODE_Lag_KF',[0:0.005:T],u0);
    
Error = max(abs(X(end,:)'-u0))
    
    
Ind = 20;
for i=1:1:length(T)
       [Tout(i) rR_Edge(:,:,i) rR(:,:,i) vR(:,:,i) aR(:,:,i) omegaR(:,:,i) gR(:,:,i) FA_R(:,:,i) Tension(:,:,i) rQ(:,i) R_KE(:,:,i) ...
       rK(:,i) vK(:,i) aK(:,i) euler(:,i) omegaK(:,i) gK(:,i) FA_K(:,i) FR_K(:,i) FG_K(:,i)  MA_K(:,i) MR_K(:,i) MG_K(:,i) MMC_K(:,i) alfa_K(i) beta_K(i)... 
       rG(:,:,i) vG(:,:,i) aG(:,:,i) omegaG(:,:,i) gG(:,:,i) FA_G(:,:,i)  FK_G(:,:,i) MA_G(:,:,i) MK_G(:,:,i)  MMC_G(:,:,i) xc_out(:,i) xs_target_out(:,i) Error_C(i)] = Fun_Post_KF(PD,PND,T(i),X(i,:)',Flag_Dim);
       
         
       if Ind==20
           Plot_GroGen_KF(X(i,:)',Tout(i),rQ(:,i),rR(:,:,i),rK(:,i),rG(:,:,i),R_KE(:,:,i),rR_Edge(:,:,i),PND,Flag_Dim,PD);
           pause(0.01)
           title('')
           Ind = 1;
          
       else
           Ind = Ind+1;
       end
       
       eta(i)  =  xc_out(4,i);  % Bridle angle eta
end
   
    
%[Kite Position, Velocity, Euler, alfa&beta, Tension, Control, Rotor angular velocity, Error_C
%Flag_Plot = [1             ,     1   ,   1  ,     1    ,    1   ,    1   ,       1               ,   1  ];
           
%Plot_Results_KF(Tout,rR_Edge, rR, vR, aR, omegaR, gR, FA_R, Tension, rQ, R_KE, ...
%rK, vK, aK, euler, omegaK, gK, FA_K, FR_K, FG_K,  MA_K, MR_K, MG_K, MMC_K, alfa_K, beta_K,... 
%rG, vG, aG, omegaG, gG, FA_G,  FK_G, MA_G, MK_G,  MMC_G, xc_out,xs_target_out,Error_C, Flag_Dim, Flag_Plot,PD)



% Control Law and trajectory

if Flag_Dim == 0
    Unit_Time         = '(\sqrt{L_0/g})$';
    Unit_Length       = '(L_0$)';
    Unit_Velocity     = '\sqrt{gL_0}$';
    Unit_deg          = '(rad)$';
    Unit_Force        = '(mg)$'
    Unit_Moment       = '(mgL_0)$';
    Unit_Ang_Velocity = '(g/L_0)^{1/2}$';    
else
    Unit_Time         = '(s)$';
    Unit_Length       = '(m)$';
    Unit_Velocity     = '(m/s)$';
    Unit_deg          = '(^\circ )$';
    Unit_Force        = '(N)$';
    Unit_Moment       = '(N m )$';
    Unit_Ang_Velocity = '(rpm)$';
end
Width = 1.5;

figure(102)
subplot(2,1,1)
hold on
plot(Tout,xc_out(4,:),'linewidth',Width)
xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
ylabel(['$\eta\ ' Unit_deg],'fontsize',12,'interpreter','latex')
grid on
set(gca,'box','on','fontsize',12)
xlim([0 max(Tout)])    

subplot(2,1,2)
hold on
hold on
plot(rK(2,:),-rK(3,:),'linewidth',Width)
plot(rK(2,end),-rK(3,end),'+')
xlabel(['$Y\ ' Unit_Length],'fontsize',12,'interpreter','latex' )
ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
grid on
set(gca,'box','on','fontsize',12)    

   
%%
 figure(104)
 
 subplot(3,1,1)
 hold on
 plot(Tout,euler(1,:),'linewidth',Width)
 %set(gca,'position',[0.13 0.73 0.8 0.23])
 xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
 ylabel(['$\psi\ ' Unit_deg],'fontsize',12,'interpreter','latex')
 grid on
 set(gca,'box','on','fontsize',12)
 xlim([0 max(Tout)])      
 
 subplot(3,1,2)
 hold on
 plot(Tout,euler(2,:),'linewidth',Width)
 %set(gca,'position',[0.13 0.43 0.8 0.23])
 xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
 ylabel(['$\theta\ ' Unit_deg],'fontsize',12,'interpreter','latex')
 grid on
 set(gca,'box','on','fontsize',12)
 xlim([0 max(Tout)])   
 
 subplot(3,1,3)
 hold on
 plot(Tout, euler(3,:),'linewidth',Width)
 xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
 ylabel(['$\phi\ ' Unit_deg],'fontsize',12,'interpreter','latex')
 grid on
 set(gca,'box','on','fontsize',12)
 xlim([0 max(Tout)])  
 ylim([0 360])  
 %set(gca,'position',[0.13 0.13 0.8 0.23])





   