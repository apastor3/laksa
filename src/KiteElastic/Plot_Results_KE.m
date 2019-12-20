function [T R_KE rK vK eulerK omegaK F_Tether  FA  alfa  beta rM T_U T_D]=Plot_Results_KE(T, R_KE, rK, vK, eulerK, omegaK, F_Tether, FA, alfa, beta,rM,T_U,T_D,Flag_Plot,Flag_Dim,PD,PND)
 
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga                                        %
% Language  : Matlab                                                         %
% Synopsis  : Compute all relevant quantities and plot                       %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:                                                                  %
%           T -> Time                                                        %
%           R_KE      -> Kite-Earth rotation matrices                        %
%           rK        -> Kite position vectors (SE components)               %
%           vK        -> Kite velocity vector (SB components)                %
%           eulerK    -> Euler angles                                        %
%           omegaK    -> Kite angular velocities (SB components)             %
%           F_Tether  -> Tether force upon the kites                         %
%           FA        -> Aerodynamic force upon the kite                     %
%           alfa      -> Angle of attack                                     %
%           beta      -> Sideslip angle                                      %
%           rM        -> Position vectors of the masses (SE components)      %
%           Flag_Plot -> Flag controlling the figures                        %
%           Flag_Dim  -> Flag controlling the figure units (0/1)             %
%           PD        -> Physical parameters                                 %
%           PND       -> Dimensionless parameters                            %
%   Outputs:                                                                 %
%           Figures with the results                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Flag_Dim==0
  Unit_Time       = '(\sqrt{L_0/g})$';
  Unit_Length     = '(L_0)$';
  Unit_Velocity   = '(\sqrt{gL_0})$';
  Unit_deg        = '(rad)$';
  Unit_Force      = '(mg)$';
else
  T        = T*sqrt(PD.Tether.L0/PD.Env.g);
  rK       = rK*PD.Tether.L0;
  vK       = vK*sqrt(PD.Env.g*PD.Tether.L0);
  eulerK   = eulerK*180/pi;
  omegaK   = omegaK*sqrt(PD.Env.g/PD.Tether.L0);
  F_Tether = F_Tether*PD.Kite.m(1)*PD.Env.g;
  T_U      = T_U*PD.Kite.m(1)*PD.Env.g; 
  T_D      = T_D*PD.Kite.m(1)*PD.Env.g;
  FA       = FA*PD.Kite.m(1)*PD.Env.g;
  alfa     = alfa*180/pi;
  beta     = beta*180/pi;
  rM       = rM*PD.Tether.L0;
  
  Unit_Time       = '(s)$';
  Unit_Length     = '(m)$)';
  Unit_Velocity   = '(m/s)$';
  Unit_deg        = '(^\circ)$';
  Unit_Force      = '(N)$';
  
end


if Flag_Plot(1) == 1 % Plot Position               %%
    
    for i=1:1:PND.Kite.Num
        figure(10+i)
        subplot(3,1,1)
        title(['Position of Kite number = ' num2str(i)])
        hold on
        plot(T,squeeze(rK(1,i,:)))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$X\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
        subplot(3,1,2)
        hold on
        plot(T,squeeze(rK(2,i,:)))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$Y\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
        subplot(3,1,3)
        hold on
        plot(T,-squeeze(rK(3,i,:)))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    end
    
    for i=1:1:PND.Tether.Num 
        figure(100+i)
        for j=1:1:PND.Tether.Mass(i)
            subplot(PND.Tether.Mass(i),1,j)
            hold on
            if j==1
               title(['Position of the masses modelling tether number = ' num2str(i)])
            end
            plot(T,squeeze(rM(1,i,j,:)))
            plot(T,squeeze(rM(2,i,j,:)))
            plot(T,-squeeze(rM(3,i,j,:)))
            xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
            Legend(1,1:length(['$X\ ' Unit_Length])) = ['$X\ ' Unit_Length]; 
            Legend(2,1:length(['$Y\ ' Unit_Length])) = ['$Y\ ' Unit_Length]; 
            Legend(3,1:length(['$H\ ' Unit_Length])) = ['$H\ ' Unit_Length]; 
            s = legend(Legend);
            set(s,'fontsize',12,'interpreter','latex')
    
        end
    end
    
end
   
if Flag_Plot(2) == 1 %Plot velocity               %%
     for i=1:1:PND.Kite.Num
        figure(200+i)
        subplot(3,1,1)
        hold on
        title(['SB-components of velocity of kite number = ' num2str(i)])
        plot(T,squeeze(vK(1,i,:)))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$v_x\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    
        subplot(3,1,2)
        hold on
        plot(T,squeeze(vK(2,i,:)))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$v_y\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    
        subplot(3,1,3)
        hold on
        plot(T,-squeeze(vK(3,i,:)))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$v_z\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    end
end

if  Flag_Plot(3) == 1 %Plot Euler Angles           %%
    for i=1:1:PND.Kite.Num
        figure(300+i)
        subplot(3,1,1)
        hold on
        title(['Euler angles of kite number' num2str(i)])
        plot(T,squeeze(eulerK(3,i,:)))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$\psi\ ' Unit_deg],'fontsize',12,'interpreter','latex')

        subplot(3,1,2)
        hold on
        plot(T,squeeze(eulerK(2,i,:)))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$\theta\ ' Unit_deg],'fontsize',12,'interpreter','latex')

        subplot(3,1,3)
        hold on
        plot(T,squeeze(eulerK(1,i,:)))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$\phi\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    end
end

if Flag_Plot(4) == 1 %Plot alfa and beta          %%
    for i=1:1:PND.Kite.Num
        figure(300+i)
        subplot(2,1,1)
        hold on
        title(['Angle of attack and sideslip angle of kite number ' num2str(i)])
        plot(T,alfa(i,:))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$\alpha\ ' Unit_deg],'fontsize',12,'interpreter','latex')

        subplot(2,1,2)
        hold on
        plot(T,beta(i,:))
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        ylabel(['$\beta\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    end
end


if Flag_Plot(5) == 1 %% Plot Leading edge tether tension        %%
    
    %Compute modulus
    for i=1:1:PND.Kite.Num
       Mod_FT = zeros(length(T),1);
       Mod_FA = zeros(length(T),1);
        
       for k=1:1:length(T)    
           Mod_FT(k) = sqrt(squeeze(F_Tether(:,i,k))'*squeeze(F_Tether(:,i,k))); 
           Mod_FA(k) = sqrt(squeeze(FA(:,i,k))'*squeeze(FA(:,i,k)));
       end
       figure(400+i)
       subplot(2,1,1)
       hold on
       title(['Tether and aerodynamic forces upon kite number ' num2str(i)])
       plot(T,Mod_FT)
       xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
       ylabel(['$F_T\ ' Unit_Force],'fontsize',12,'interpreter','latex')

       subplot(2,1,2)
       hold on
       plot(T,Mod_FA)
       xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
       ylabel(['$F_A\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    end  
    
end







end