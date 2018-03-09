function Plot_Results_KA(T, RBE, rk, vk, ak, euler, omega, omega_p, Lambda, FAP, FAM, MAP, MAM, FA, MA, W, alfa, beta, LP, LM,Flag_Dim,Flag_Plot)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Plot the results                                               %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %%
% Inputs:  T                 -> Time                                %%
%          Kinematics quantities of the kite                        %%
%                      RBE   -> Body-Earth rotation matrix          %%
%                      rk    -> Earth components of the position,   %%
%                      vk       velocity and acceleration vectors   %%
%                      ak       of the center of mass of the kite   %% 
%                      euler -> Euler angles                        %%
%                      omega -> Earth components of the kite        %%
%                      omega_p  angular velocity and acceleration   %%
%          Kite forces and moments                                  %%
%                      Lambda -> Tether Tensions along              %%
%                                tether directions                  %%
%                      FAP   -> Earth components of the Forces      %%
%                      FAM      exerted by the tethers linked at    %%
%                               points A plus and A Minus           %%
%                      MAP  ->  Earth components of the Torques     %%
%                      MAM      exerted by the tethers linked at    %%
%                               points A plus and A minus           %%
%                      FA   ->  Earth components of the Aerodynamic %% 
%                      MA       force and torque                    %% 
%                      W    -> Earth components of the kite weight  %%
%           Others                                                  %%
%                      alfa -> Angle of attack                      %%  
%                      beta -> Sideslip angle                       %%
%                      LP   -> Length of the A plus-tether          %%
%                      LM   -> Length of the A minus-tether         %%
%         Flag_Dim  -> Flag controlling the dimension of the inputs %%
%                                                                   %%
%         Flag_Plot -> vector with 0/1 to set plotting options      %%
%                                                                   %%
%                      Flag_Plot(1) = 1 Plot Position               %%
%                      Flag_Plot(2) = 1 Plot velocity               %%
%                      Flag_Plot(3) = 1 Plot Euler Angles           %%
%                      Flag_Plot(4) = l Plot alfa and beta          %%
%                      Flag_Plot(5) = 1 Plot Tether tensions        %%
%                      Flag_Plot(6) = 1 Plot Tether Lengths         %%
% Outputs: Figures with the results                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Flag_Dim == 0
    Unit_Time       = '(\sqrt{L_0/g})$';
    Unit_Length     = '(L_0$)';
    Unit_Velocity   = '\sqrt{gL_0}$';
    Unit_deg        = '(rad)$';
    Unit_Force      = '(mg)$'
else
    Unit_Time       = '(s)$';
    Unit_Length     = '(m)$';
    Unit_Velocity   = 'm/s$';
    Unit_deg        = '^\circ $';
    Unit_Force      = '(N)$';
    
end


if Flag_Plot(1) == 1 % Plot Position               %%
    
    figure(101)
    subplot(3,1,1)
    hold on
    plot(T,rk(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$X\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,2)
     hold on
    plot(T,rk(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$Y\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,3)
    hold on
    plot(T,-rk(3,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
    figure(102)
    subplot(2,1,1)
    hold on
    hold on
    plot(rk(2,:),-rk(3,:))
    plot(rk(2,1),-rk(3,1),'go')
    plot(rk(2,end),-rk(3,end),'+r')
    xlabel(['$Y\ ' Unit_Length],'fontsize',12,'interpreter','latex' )
    ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
    figure(102)
    subplot(2,1,2)
    hold on
    plot(rk(1,:),-rk(3,:))
    plot(rk(1,1),-rk(3,1),'go')
    plot(rk(1,end),-rk(3,end),'+r')
    xlabel(['$X\ ' Unit_Length],'fontsize',12,'interpreter','latex' )
    ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
end
   
if Flag_Plot(2) == 1 %Plot velocity               %%
    figure(103)
    subplot(3,1,1)
    hold on
    plot(T,vk(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$V_x\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,2)
    hold on
    plot(T,vk(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$V_y\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,3)
    hold on
    plot(T, vk(3,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$V_z\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
end

if  Flag_Plot(3) == 1 %Plot Euler Angles           %%
    
    figure(104)
    hold on
    subplot(3,1,1)
    plot(T,euler(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\psi\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,2)
    hold on
    plot(T,euler(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\theta\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,3)
    hold on
    plot(T, euler(3,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\phi\ ' Unit_deg],'fontsize',12,'interpreter','latex')
end

if Flag_Plot(4) == 1 %Plot alfa and beta          %%
    figure(105)
    subplot(2,1,1)
    hold on
    plot(T,alfa)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\alpha\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    
    subplot(2,1,2)
    hold on
    plot(T,beta)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\beta\ ' Unit_deg],'fontsize',12,'interpreter','latex')
end


if Flag_Plot(5) == 1 %% Plot Tether tensions        %%
    figure(106)
    subplot(2,1,1)
    hold on
    plot(T,Lambda(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T^+\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(2,1,2)
    hold on
    plot(T,Lambda(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T^-\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
end


if Flag_Plot(6) == 1 % Plot Tether Lengths         %%
    figure(107)
    subplot(2,1,1)
    hold on
    plot(T,LP)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$L^+\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
    subplot(2,1,2)
    hold on
    plot(T,LM)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$L^-\ ' Unit_Length],'fontsize',12,'interpreter','latex')
end


end