function Plot_Results_KT(T_out,RBE,rk,vk,ak,euler,omega,omega_p,rA,rC,Lambda,TA,TC,FA,MA,W,alfa,beta,Error,Flag_Dim,Flag_Plot,PND)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Jose Antonio Serrano-Iglesias         %
% Language  : Matlab                                                         %
% Synopsis  : Plot the results                                               %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %%
% Inputs:  T_out             -> Time                                %%
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
%                      TA     -> SE-comp Forces exerted by the      %%
%                      TC          tethers linked at A+- and C+-    %%
%                      FA     -> SB-comp of the Aerodynamic force   %% 
%                      MA          and torque                       %% 
%                      W      -> SE-comp kite weights               %%
%           Others                                                  %%
%                      alfa -> Angle of attack                      %%  
%                      beta -> Sideslip angle                       %%                      
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
    Unit_Force      = '(mg)$';
    Unit_Torque     = '(mg\cdot L_0)$';
    Unit_AngVel     = '(rad\cdot \sqrt{L_0/g})$';
else
    Unit_Time       = '(s)$';
    Unit_Length     = '(m)$';
    Unit_Velocity   = 'm/s$';
    Unit_deg        = '^\circ $';
    Unit_Force      = '(N)$'; 
    Unit_Torque     = '(N\cdot m)$';
    Unit_AngVel     = '(rad\cdot s)$';
end

for i = 1:1:PND.Kite.N
    
  if     Flag_Plot(1) == 1 % Plot Position %
    
         figure()
         subplot(3,1,1)
         plot(T_out,squeeze(rk(1,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$X\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
         subplot(3,1,2)
         plot(T_out,squeeze(rk(2,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$Y\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
         subplot(3,1,3)
         plot(T_out,squeeze(-rk(3,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
         figure()
         subplot(2,1,1)
         plot(squeeze(rk(2,i,:))',squeeze(-rk(3,i,:))')
         hold on
         plot(squeeze(rk(2,i,1))',squeeze(-rk(3,i,1))','go')
         plot(squeeze(rk(2,i,end))',squeeze(-rk(3,i,end))','+r')
         xlabel(['$Y\ ' Unit_Length],'fontsize',12,'interpreter','latex' )
         ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
         subplot(2,1,2)
         plot(squeeze(-rk(1,i,:))',squeeze(-rk(3,i,:))')
         hold on
         plot(squeeze(-rk(1,i,1))',squeeze(-rk(3,i,1))','go')
         plot(squeeze(-rk(1,i,end))',squeeze(-rk(3,i,end))','+r')
         xlabel(['$X\ ' Unit_Length],'fontsize',12,'interpreter','latex' )
         ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    
  end
  if Flag_Plot(2) == 1 % Plot velocity %
         figure()
         subplot(3,1,1)
         plot(T_out,squeeze(vk(1,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$V_x\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
     
         subplot(3,1,2)
         plot(T_out,squeeze(vk(2,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$V_y\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    
         subplot(3,1,3)
         plot(T_out,squeeze(vk(3,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$V_z\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')

  end
  if Flag_Plot(3) == 1 % Plot Euler Angles %
    
         figure()
         subplot(3,1,1)
         plot(T_out,squeeze(euler(1,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$\psi\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    
         subplot(3,1,2)
         plot(T_out,squeeze(euler(2,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$\theta\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    
         subplot(3,1,3)
         plot(T_out,squeeze(euler(3,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$\phi\ ' Unit_deg],'fontsize',12,'interpreter','latex')

  end
  if Flag_Plot(4) == 1 % Plot alfa and beta %
    
         figure()
         subplot(2,1,1)
         plot(T_out,alfa(i,:))
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$\alpha\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    
         subplot(2,1,2)
         plot(T_out,beta(i,:))
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$\beta\ ' Unit_deg],'fontsize',12,'interpreter','latex')

  end
  if Flag_Plot(5) == 1 % Plot Tether tensions %
    
         figure()
         subplot(2,1,1)
         plot(T_out,squeeze(Lambda(1,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$T^+\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
         subplot(2,1,2)
         plot(T_out,squeeze(Lambda(2,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$T^-\ ' Unit_Force],'fontsize',12,'interpreter','latex')
         
  end
  if Flag_Plot(10) == 1 % Plot Aerodynamic Moments %
    
         figure()
         subplot(3,1,1)
         plot(T_out,squeeze(MA(1,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$M_x\ ' Unit_Torque],'fontsize',12,'interpreter','latex')
    
         subplot(3,1,2)
         plot(T_out,squeeze(MA(2,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$M_y\ ' Unit_Torque],'fontsize',12,'interpreter','latex')     
    
         subplot(3,1,3)
         plot(T_out,squeeze(MA(3,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$M_z\ ' Unit_Torque],'fontsize',12,'interpreter','latex')  
         
  end
  if Flag_Plot(11) == 1 % Plot Error %
    
         figure()
         plot(T_out,Error)
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$Error\ -$'],'fontsize',12,'interpreter','latex')
    
  end
  if Flag_Plot(12) == 1 % Plot Angular Valocity %
    
         figure()
         subplot(3,1,1)
         plot(T_out,squeeze(omega(1,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$\omega_x\ ' Unit_AngVel],'fontsize',12,'interpreter','latex')
    
         subplot(3,1,2)
         plot(T_out,squeeze(omega(2,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$\omega_y\ ' Unit_AngVel],'fontsize',12,'interpreter','latex')     
    
         subplot(3,1,3)
         plot(T_out,squeeze(omega(3,i,:))')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$\omega_z\ ' Unit_AngVel],'fontsize',12,'interpreter','latex')     

  end
end
