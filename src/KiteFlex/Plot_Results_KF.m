function Plot_Results_KS(T,rR_Edge, rR, vR, aR, omegaR, gR, FA_R, Tension, rQ, R_KE, ...
          rK, vK, aK, euler, omegaK, gK, FA_K, FR_K, FG_K,  MA_K, MR_K, MG_K, MMC_K, alfa_K, beta_K,... 
          rG, vG, aG, omegaG, gG, FA_G,  FK_G, MA_G, MK_G,  MMC_G, xc_out,xs_target_out,Error_C, Flag_Dim, Flag_Plot,PD)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Plot the simulation results                                    %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %%
% Inputs: Dimensionless or with dimensions depending on Flag_Dim          %%
%          T         -> Time                                               %%
%          rR_Edge   -> Position of rods tips  (SE components)             %%
%          rR        -> Position of rods' centers (SE components)          %%
%          vR        -> Velocity of the rods (SE components)               %%
%          aR        -> Acceleration of the rods (SE components)           %%
%          omegaR    -> Angular velocity of the rods (SR components)       %%
%          gR        -> Angular acceleration of the rods (SR components)   %%
%          FA_R      -> Aerodynamic force upon the rods  (SE components)   %%
%          Tension   -> Tension at the rods  (SE components)               %%
%          rQ        -> Position of vector Q  (SE components)              %%
%          R_KE      -> SK-SE rotation matrix                              %%
%          rK        -> Kite position  (SE components)                     %%
%          vK        -> Kite velocity  (SE components)                     %%
%          aK        -> Kite acceleration  (SE components)                 %%
%          euler     -> Kite Euler angles                                  %%
%          omegaK    -> Kite angular velocity  (SK components)             %%
%          gK        -> Kite angular acceleration  (SK components)         %%
%          FA_K      -> Kite aerodynamic force  (SE components)            %%
%          FR_K      -> Rod force upon the kite  (SE components)           %%
%          FG_K      -> Rotor force upon the kite  (SE components)         %%
%          MA_K      -> Kite aerodynamic torque   (SK components)          %%
%          MR_K      -> Rotor torque upon the kite   (SK components)       %%
%          MG_K      -> Rods torque upon the kite  (SK components)         %%
%          MMC_K     -> Motor Cont. torque upon the kite (SK components)   %%
%          alfa_K    -> Kite angle of attack                               %%
%          beta_K    -> Kite sideslip angles                               %%
%          rG        -> Rotor position  (SE components)                    %%
%          vG        -> Rotor velocity (SE components)                     %%
%          aG        -> Rotor acceleration (SE components)                 %%
%          omegaG    -> Rotor angular velocity   (SK components)           %%
%          gG        -> Rotor angular acceleration  (SK components)        %%
%          FA_G      -> Rotor aerodynamic force (SE components)            %%
%          FK_G      -> Kite forceupon the rotors (SE components)          %%
%          MA_G      -> Rotor Aerodynamic torque   (SK components)         %%
%          MK_G      -> Kite torque upon the rotors  (SK components)       %%
%          MMC_G     -> Motor Cont. torque upon the kite   (SK components) %%
%          xc_out    -> Control vector                                     %%  
%          xs_target -> Target State vector if Close-Loop                 %%  
%          Error_C   -> Error compared with classical formulation          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LWidth = 1; % Line Width

if Flag_Dim == 0
    Unit_Time         = '(\sqrt{L_0/g})$';
    Unit_Length       = '(L_0$)';
    Unit_Velocity     = '\sqrt{gL_0}$';
    Unit_deg          = '(rad)$';
    Unit_Force        = '(mg)$';
    Unit_Moment       = '(mgL_0)$';
    Unit_Ang_Velocity = '(g/L_0)^{1/2}$';    
else
    Unit_Time         = '(s)$';
    Unit_Length       = '(m)$';
    Unit_Velocity     = '(m/s$)';
    Unit_deg          = '(^\circ)$';
    Unit_Force        = '(N)$';
    Unit_Moment       = '(N m )$';
    Unit_Ang_Velocity = '(rad/s)$';
end

NR = PD.Num.N ;   % Number of Rods
NG = PD.Gen.Num;  % Number of Generators
NT = 4+NG+3;      % Control variables


if Flag_Plot(1) == 1 % Plot Position               %%
    
    figure(101)
    subplot(3,1,1)
    hold on
    plot(T,rK(1,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$X\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(3,1,2)
    hold on
    plot(T,rK(2,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$Y\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(3,1,3)
    hold on
    plot(T,-rK(3,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    figure(102)
    subplot(2,1,1)
    hold on
    hold on
    plot(rK(2,:),-rK(3,:),'linewidth',LWidth)
    plot(rK(2,1),-rK(3,1),'go')
    plot(rK(2,end),-rK(3,end),'+r')
    xlabel(['$Y\ ' Unit_Length],'fontsize',12,'interpreter','latex' )
    ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    figure(102)
    subplot(2,1,2)
    hold on
    plot(rK(1,:),-rK(3,:),'linewidth',LWidth)
    plot(rK(1,1),-rK(3,1),'go')
    plot(rK(1,end),-rK(3,end),'+r')
    xlabel(['$X\ ' Unit_Length],'fontsize',12,'interpreter','latex' )
    ylabel(['$H\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)    
end
   
if Flag_Plot(2) == 1 %Plot velocity               %%
    figure(103)
    subplot(3,1,1)
    hold on
    plot(T,vK(1,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$V_x\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(3,1,2)
    hold on
    plot(T,vK(2,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$V_y\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(3,1,3)
    hold on
    plot(T, vK(3,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$V_z\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
end

if  Flag_Plot(3) == 1 %Plot Euler Angles           %%
    
    figure(104)
   
    subplot(3,1,1)
    plot(T,euler(1,:),'linewidth',LWidth)
    hold on
    if PD.Control.Type==5 % Close-Loop
       plot(T,xs_target_out(2*NR+2,:),'r--','linewidth',LWidth)
    end
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\psi\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(3,1,2)
    hold on
    plot(T,euler(2,:),'linewidth',LWidth)
    if PD.Control.Type==5 % Close-Loop
       plot(T,xs_target_out(2*NR+1,:),'r--','linewidth',LWidth)
    end
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\theta\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(3,1,3)
    hold on
    plot(T, euler(3,:),'linewidth',LWidth)
    if PD.Control.Type==5 % Close-Loop
       plot(T,xs_target_out(2*NR+3,:),'r--','linewidth',LWidth)
    end
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\phi\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
end

if Flag_Plot(4) == 1 %Plot alfa and beta          %%
    figure(105)
    subplot(2,1,1)
    hold on
    plot(T,alfa_K,'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\alpha\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(2,1,2)
    hold on
    plot(T,beta_K,'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\beta\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
end


if Flag_Plot(5) == 1 %% Plot Tension        %%
    
    %Compute the modulus of the tension
    for i=1:1:length(Tension(1,:,1))    % Bar loop
        for j=1:1:length(Tension(1,1,:)) % Temporal Loop
            Mod(i,j) = sqrt(Tension(:,i,j)'*Tension(:,i,j));
        end
    end
    
    figure(106)
    hold on
    plot(T,Mod,'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)  
end


if Flag_Plot(6) == 1 %% Plot Control        %%
    
    figure(107)
    subplot(4,1,1)
    hold on
    plot(T,xc_out(1,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$l_R\ ' Unit_Length] ,'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(4,1,2)
    hold on
    plot(T,xc_out(2,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$l_B\ ' Unit_Length],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(4,1,3)
    hold on
    plot(T,xc_out(3,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\delta\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(4,1,4)
    hold on
    plot(T,xc_out(4,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\eta\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    if PD.Gen.Num>0
        figure(108)
        for i=1:1:PD.Gen.Num
            subplot(PD.Gen.Num/2,2,i)
            hold on
            plot(T,xc_out(4+i,:),'linewidth',LWidth)
            xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
            ylabel(['$m_j\ ' Unit_Moment] ,'fontsize',12,'interpreter','latex')
        end
    end
    
    figure(109)
    subplot(3,1,1)
    hold on
    plot(T,xc_out(4+PD.Gen.Num+1,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\delta_a\ ' Unit_deg] ,'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(3,1,2)
    hold on
    plot(T,xc_out(4+PD.Gen.Num+2,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\delta_r\ ' Unit_deg] ,'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
    
    subplot(3,1,3)
    hold on
    plot(T,xc_out(4+PD.Gen.Num+3,:),'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\delta_e\ ' Unit_deg] ,'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
end


if Flag_Plot(7) == 1 %% Absolute Angular velocity of the rotor 
    
    for i=1:1:PD.Gen.Num
         figure(110)
         subplot(PD.Gen.Num/2,2,i)
         hold on
         plot(T,squeeze(omegaG(1,i,:)),'linewidth',LWidth)
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
         ylabel(['$\omega_{xK}\ '  Unit_Ang_Velocity] ,'fontsize',12,'interpreter','latex')
         grid on
         set(gca,'box','on','fontsize',12)
    end    
    
   
    
end


if Flag_Plot(8) == 1 %% Check consistency with classical mechanics   %%
    
    figure(111)
    hold on
    plot(T,Error_C,'linewidth',LWidth)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$Error$ '] ,'fontsize',12,'interpreter','latex')
    grid on
    set(gca,'box','on','fontsize',12)
end



  

end