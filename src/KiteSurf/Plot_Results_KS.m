function Plot_Results_KS(T, RBE, rk, vk, ak, euler, omega, omega_p, Lambda, FAP, FAM, MAP, MAM, FBP, FBM, MBP, MBM,...
                         FA, MA, W, alfa, beta,Rp, Rm, Elong_p, Elong_m,xc,Error0 ,Flag_Dim,Flag_Plot)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Author    : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Plot all the relevant quantities                               %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
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
%                               points  plus and A minus           %%
%                      FBP   -> Earth components of the Forces      %%
%                      FBM      exerted by the tethers linked at    %%
%                               points B plus and A Minus           %%
%                      MBP  ->  Earth components of the Torques     %%
%                      MBM      exerted by the tethers linked at    %%
%                               points B plus and A minus           %%

%                      FA   ->  Earth components of the Aerodynamic %% 
%                      MA       force and torque                    %% 
%                      W    -> Earth components of the kite weight  %%
%           Others                                                  %%
%                      alfa -> Angle of attack                      %%  
%                      beta -> Sideslip angle                       %%
%           Elastic tethers                                         %%
%                      Rp,Rm-> position vector of elastic tethers   %%
%                      Elong_p, Elong_m-> elongation of elastic     %%
%                                tethers                            %%
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
    Unit_Velocity   = '(\sqrt{gL_0})$';
    Unit_deg        = '(rad)$';
    Unit_Force      = '(mg)$';
    Unit_Omega      = '(\sqrt{g/L0})$';
    Unit_Acce       = '(g)$';
    Unit_Alp        = '(g/L_0)$';
else
    Unit_Time       = '(s)$';
    Unit_Length     = '(m)$';
    Unit_Velocity   = '(m/s)$';
    Unit_deg        = '(^\circ )$';
    Unit_Force      = '(N)$';
    Unit_Omega      = '(rad/s)$';
    Unit_Acce       = '(m/s^2)$';
    Unit_Alp        = '(rad/s^2)$';
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
    for i=1:1:length(T)
        vkb(:,i) =  squeeze(RBE(:,:,i))*vk(:,i);   % 
        akb(:,i) =  squeeze(RBE(:,:,i))*ak(:,i);   % 
        Alp(:,i) =  squeeze(RBE(:,:,i))*omega_p(:,i);
    end
    
    figure(103)
    subplot(3,1,1)
    hold on
    plot(T,vkb(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$V_{Bx}\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,2)
    hold on
    plot(T,vkb(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$V_{By}\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,3)
    hold on
    plot(T, vkb(3,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$V_{Bz}\ ' Unit_Velocity],'fontsize',12,'interpreter','latex')
    
    
    figure(1003)
    subplot(3,1,1)
    hold on
    plot(T,akb(1,:),'b-')
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$a_{Bx}\ ' Unit_Acce],'fontsize',12,'interpreter','latex')

    subplot(3,1,2)
    hold on
    plot(T,akb(2,:),'b-')
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$a_{By}\ ' Unit_Acce],'fontsize',12,'interpreter','latex')

    subplot(3,1,3)
    hold on
    plot(T, akb(3,:),'b-')
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$a_{Bz}\ ' Unit_Acce],'fontsize',12,'interpreter','latex')
    
    figure(8003)
    subplot(3,1,1)
    hold on
    plot(T,Alp(1,:),'b-')
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\alpha_{Bx}\ ' Unit_Alp],'fontsize',12,'interpreter','latex')

    subplot(3,1,2)
    hold on
    plot(T, Alp(2,:),'b-')
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\alpha_{By}\ ' Unit_Alp],'fontsize',12,'interpreter','latex')

    subplot(3,1,3)
    hold on
    plot(T,Alp(3,:),'b-')
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\alpha_{Bz}\ ' Unit_Alp],'fontsize',12,'interpreter','latex')
    
end



if  Flag_Plot(3) == 1 %Plot Euler Angles           %%
    
    figure(104)
    subplot(3,1,1)
    hold on
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


if Flag_Plot(5) == 1 %% Plot inelastic Tether tensions        %%
    for i=1:1:length(T)  
       % Tensions applied at the kite points A+-
       T_Ap(:,i) =  squeeze(RBE(:,:,i))*FAP(:,i);   % Tension at B+ (components in the Earth frame)
       T_Am(:,i) =  squeeze(RBE(:,:,i))*FAM(:,i);   % Tension at B- (components in the Earth frame)
    end    
    
    figure(1061)
    subplot(3,1,1)
    hold on
    plot(T,T_Ap(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{A^+x}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,2)
    hold on
    plot(T,T_Ap(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{A^+y}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,3)
    hold on
    plot(T,T_Ap(3,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{A^+z}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    figure(1062)
    subplot(3,1,1)
    hold on
    plot(T,T_Am(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{A^-x}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,2)
    hold on
    plot(T,T_Am(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{A^-y}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,3)
    hold on
    plot(T,T_Am(3,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{A^-z}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
end


if Flag_Plot(6) == 1 %% Plot Control        %%
    figure(107)
    subplot(2,1,1)
    hold on
    plot(T,xc(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel('$u_p$','fontsize',12,'interpreter','latex')
    
    subplot(2,1,2)
    hold on
    plot(T,xc(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\nu\ ' Unit_deg],'fontsize',12,'interpreter','latex')
    
end


if Flag_Plot(8) == 1 %% Tension of elastic tether   %%
    
     for i=1:1:length(T)  
       % Tensions applied at the kite points B+-
       T_Bp(:,i) =  squeeze(RBE(:,:,i))*FBP(:,i);   % Tension at B+ (SB components)
       T_Bm(:,i) =  squeeze(RBE(:,:,i))*FBM(:,i);   % Tension at B- (SB components)
    end    
    figure(1081)
    subplot(4,1,1)
    hold on
    plot(T,T_Bp(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{B^+x}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(4,1,2)
    hold on
    plot(T,T_Bp(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{B^+y}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(4,1,3)
    hold on
    plot(T,T_Bp(3,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{B^+z}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(4,1,4)
    hold on
    plot(T,sqrt(T_Bp(1,:).^2+T_Bp(2,:).^2+T_Bp(3,:).^2))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{B^+}\ ' Unit_Force],'fontsize',12,'interpreter','latex')

    
    figure(1082)
    subplot(4,1,1)
    hold on
    plot(T,T_Bm(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{B^-x}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(4,1,2)
    hold on
    plot(T,T_Bm(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{B^-y}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(4,1,3)
    hold on
    plot(T,T_Bm(3,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{B^-z}\ ' Unit_Force],'fontsize',12,'interpreter','latex')
    
    subplot(4,1,4)
    hold on
    plot(T,sqrt(T_Bm(1,:).^2+T_Bm(2,:).^2+T_Bm(3,:).^2))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$T_{B^-}\ ' Unit_Force],'fontsize',12,'interpreter','latex')

    
end


if Flag_Plot(9) == 1 %% Elongation   %%
    
    figure(109)
    subplot(2,1,1)
    hold on
    plot(T,Elong_p)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel('$Elongation_{B^+}$','fontsize',12,'interpreter','latex')
    
    subplot(2,1,2)
    hold on
    plot(T,Elong_m)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel('$Elogantion_{B^-}$','fontsize',12,'interpreter','latex')
    
end


if Flag_Plot(10) == 1 %% Moments   %%
    
    for i=1:1:length(T)  
        MAP(:,i)  = squeeze(RBE(:,:,i))*MAP(:,i);
        MAM(:,i)  = squeeze(RBE(:,:,i))*MAM(:,i);
        MBP(:,i)  = squeeze(RBE(:,:,i))*MBP(:,i);
        MBM(:,i)  = squeeze(RBE(:,:,i))*MBM(:,i);
        MA(:,i)   = squeeze(RBE(:,:,i))*MA(:,i);
    end    
    
    for i=1:1:3
        figure(1110) % Aerodynamic Torque
        subplot(3,1,i)
        if i==1
            title('Aerodynamic Torque')
        end
        hold on
        plot(T,MA(i,:),'b')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        if i==1
           ylabel('$M_{Bx}$','fontsize',12,'interpreter','latex')
        else
            if i==2
               ylabel('$M_{By}$','fontsize',12,'interpreter','latex')
            else
               ylabel('$M_{Bz}$','fontsize',12,'interpreter','latex')
            end
        end
         
        figure(1111) % Inelastic Tether Torque 
        subplot(3,1,i)
        if i==1
            title('Inelastic Tether Torque')
        end
        hold on
        plot(T,MAP(i,:)+MAM(i,:),'b')
         xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        if i==1
           ylabel('$M_{Bx}$','fontsize',12,'interpreter','latex')
        else
            if i==2
               ylabel('$M_{By}$','fontsize',12,'interpreter','latex')
            else
               ylabel('$M_{Bz}$','fontsize',12,'interpreter','latex')
            end
        end
        
        figure(1112) % Elastic Tether Torque 
        subplot(3,1,i)
        if i==1
            title('Elastic Tether Torque')
        end
        hold on
        plot(T,MBP(i,:)+MBM(i,:),'b')
       
       
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        if i==1
           ylabel('$M_{Bx}$','fontsize',12,'interpreter','latex')
        else
            if i==2
               ylabel('$M_{By}$','fontsize',12,'interpreter','latex')
            else
               ylabel('$M_{Bz}$','fontsize',12,'interpreter','latex')
            end
        end
        
        
        figure(1113) % Total Toruqe 
        subplot(3,1,i)
        if i==1
            title('Total Torque')
        end
        hold on
        plot(T,MBP(i,:)+MBM(i,:)+MAP(i,:)+MAM(i,:)+MA(i,:),'b')
       
       
        xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
        if i==1
           ylabel('$M_{Bx}$','fontsize',12,'interpreter','latex')
        else
            if i==2
               ylabel('$M_{By}$','fontsize',12,'interpreter','latex')
            else
               ylabel('$M_{Bz}$','fontsize',12,'interpreter','latex')
            end
        end
                
    end 
end

if Flag_Plot(11) == 1 %% Error -> comparison with F = ma and dL/dt = M   %%
    figure(111)
    hold on
    semilogy(T,Error0)
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel('$Error$','fontsize',12,'interpreter','latex')
end

if Flag_Plot(12) == 1 %% SB-Components of the angular velocity
    
    for i=1:1:length(T)  
        omega_B(:,i)  = squeeze(RBE(:,:,i))*omega(:,i);
    end   
    
    figure(112)
    subplot(3,1,1)
    hold on
    plot(T, omega_B(1,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\omega_{Bx}\ ' Unit_Omega],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,2)
    hold on
    plot(T, omega_B(2,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\omega_{By}\ ' Unit_Omega],'fontsize',12,'interpreter','latex')
    
    subplot(3,1,3)
    hold on
    plot(T,  omega_B(3,:))
    xlabel(['$Time\ ' Unit_Time],'fontsize',12,'interpreter','latex' )
    ylabel(['$\omega_{Bz}\ ' Unit_Omega],'fontsize',12,'interpreter','latex')
end



  
  

end