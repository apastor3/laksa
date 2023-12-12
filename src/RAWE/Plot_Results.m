function Plot_Results(Fig,Ax,Flag_Dim,Flag_Movie,NP,PP,PD,tau,xs)

N_fe = NP(1);
N_nc = NP(4);
Nu   = 3*(N_fe -1 ) + N_fe + 1;

[s_fe s_nc phi phi_s phi_ss  phi_sss varphi varphi_s varphi_ss Shape Shape_s Shape_ss Shape_sss B_alpha_inv] = Fun_Global(NP(1),NP(2));       
[P0_Inv Q0_Inv p0 P4 Q1 P00 Q00] = Compute_Matrix(NP(1),100);



if Flag_Dim==1   
        TLab      = '$Time\ (s)$';
        xLab      = '$x\ (m)$';
        yLab      = '$y\ (m)$';
        zLab      = '$z\ (m)$';
        vxLab     = '$v_x\ (m/s)$';
        vyLab     = '$v_y\ (m/s)$';
        vzLab     = '$v_z\ (m/s)$';
        alpLab    = '$\alpha\ (deg)$';
        alp_tLab  = '$\dot{\alpha}\ (rpm)$';
        MLab      = '$M_0\ (N m)$';
        FLab      = '$T (N)$';
else
        TLab      = '$\tau$';
        xLab      = '$x $';
        yLab      = '$y $';
        zLab      = '$z $';
        vxLab     = '$v_x $';
        vyLab     = '$v_y $';
        vzLab     = '$v_z $';
        alpLab    = '$\alpha$';
        alp_tLab  = '$\alpha_t$';
        MLab      = '$M (N m)$';
        FLab      = '$T$ ';
end

 % Make a movie  
    r       = zeros(length(tau),NP(1)+1,3);
    r_t     = zeros(length(tau),NP(1)+1,3);
    alpha   = zeros(length(tau),NP(1)+1);
    alpha_t = zeros(length(tau),NP(1)+1);
    
    Ind     = 1;
    for i=1:1:length(tau)  
        [r(i,:,:) alpha(i,:) r_t(i,:,:) alpha_t(i,:) ] = From_xs_to_Physical(tau(i),xs(i,:)',NP,PP);
        
        [r0 rN alpha_s_N r0_t rN_t rN_tt m_N_tt] = Compute_BC(tau(i),PP);
        if i<=2 
            E_m_0_tt = 0;
        else
            alpha_tt0 = alpha_tt(i-1) + (alpha_tt(i-1) - alpha_tt(i-2))/(tau(i-1) - tau(i-2))*( tau(i) - tau(i-1));                
            alpha_ttt0 = alpha_ttt(i-1) + (alpha_ttt(i-1) - alpha_ttt(i-2))/(tau(i-1) - tau(i-2))*( tau(i) - tau(i-1));                
            
            E_m_0_tt = PP(10)*alpha_tt0 + PP(11)*alpha_ttt0;
        end
        m_N(i)     = alpha_s_N; 
        m_0(i)     = m_N(i) + xs(i,2*Nu+1);
        m_0_tt     = m_N_tt + E_m_0_tt; 
        
        [V_alpha V_alpha_s V_alpha_ss V_alpha_tt] = Fun_Global_Valpha(NP, m_0(i), m_N(i),B_alpha_inv,Shape,Shape_s,Shape_ss,m_0_tt,m_N_tt );

        [DF1 DG1 r_N1 alpha_N1 r_t1 alpha_t1 m_01 kappa1   tau11 Ten1(:,i) F_NL1 RHS_r1 RHS_alpha1] = Fun_ODE_RAWE(tau(i),xs(i,:)',E_m_0_tt);

        alpha_t(i)   = xs(i,Nu + 3*(N_fe-1)+1);
        alpha_tt(i)  = DF1(Nu + 3*(N_fe-1)+1,1);
        if i==1
            alpha_ttt(1) = 0;
        else
            alpha_ttt(i) = (alpha_tt(i) - alpha_tt(i-1))/(tau(i)-tau(i-1));
        end
  
        if Ind == 100 && Flag_Movie==1
            Plot_Rotary(Fig,Ax,NP,squeeze(r(i,:,:)),squeeze(alpha(i,:)),phi,phi_s,phi_ss,phi_sss,varphi,varphi_s,varphi_ss,V_alpha,V_alpha_s,V_alpha_ss);
            Ind = 1;
         
        else
            Ind = Ind+1;
        end
    end
    
    if Flag_Dim==1
        tau       = tau*PD(2);
        r         = PD(1)*r;
        r_t       = PD(1)/PD(2)*r_t;
        alpha_t   = (1/PD(2))*(60/(2*acos(-1)))*alpha_t;
        M_0       = m_0*PD(5)/PD(1);
        M_N       = m_N*PD(5)/PD(1);
        Ten1      = Ten1*PD(3)*PD(1)/PD(2)^2;
    else
        M_0       = m_0;
        M_N       = m_N;
    end
    
    figure(12)
    subplot(3,1,1)
    hold on
    plot(tau,r(:,:,1),'linewidth',1.5)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel(xLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
    
    DeltaY = 0.1*max([abs(min(min(r(:,:,1)))) abs(max(max(r(:,:,1))))] );
    %axis([min(tau) max(tau) min(min(r(:,:,1))) - DeltaY  max(max(r(:,:,1))) + DeltaY ]);
    
    
    subplot(3,1,2)
      hold on
    plot(tau,r(:,:,2),'linewidth',1.5)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel(yLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
    DeltaY = 0.1*max([abs(min(min(r(:,:,2)))) abs(max(max(r(:,:,2))))] );
    %axis([min(tau) max(tau) min(min(r(:,:,2))) - DeltaY  max(max(r(:,:,2))) + DeltaY ]);

    
    
    subplot(3,1,3)
      hold on
    plot(tau,r(:,:,3),'linewidth',1.5)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel(zLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
    DeltaY = 0.1*max([abs(min(min(r(:,:,3)))) abs(max(max(r(:,:,3))))] );
    %axis([min(tau) max(tau) min(min(r(:,:,3))) - DeltaY  max(max(r(:,:,3))) + DeltaY ]);

    
    figure(13)
    subplot(3,1,1)
      hold on
    plot(tau,r_t(:,:,1),'linewidth',1.5)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel(vxLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
    DeltaY = 0.1*max([abs(min(min(r_t(:,:,1)))) abs(max(max(r_t(:,:,1))))] );
    %axis([min(tau) max(tau) min(min(r_t(:,:,1))) - DeltaY  max(max(r_t(:,:,1))) + DeltaY ]);

    
    subplot(3,1,2)
      hold on
    plot(tau,r_t(:,:,2),'linewidth',1.5)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel(vyLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
    DeltaY = 0.1*max([abs(min(min(r_t(:,:,2)))) abs(max(max(r_t(:,:,2))))] );
    %axis([min(tau) max(tau) min(min(r_t(:,:,2))) - DeltaY  max(max(r_t(:,:,2))) + DeltaY ]);

    
    subplot(3,1,3)
      hold on
    plot(tau,r_t(:,:,3),'linewidth',1.5)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel(vzLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
    DeltaY = 0.1*max([abs(min(min(r_t(:,:,3)))) abs(max(max(r_t(:,:,3))))] );
    %axis([min(tau) max(tau) min(min(r_t(:,:,3))) - DeltaY  max(max(r_t(:,:,3))) + DeltaY ]);
   
    
    figure(14)
    subplot(3,1,1)
    hold on
    plot(tau,alpha/(2*pi),'linewidth',1.5)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel('Turn No.','fontsize',12)
    set(gca,'fontsize',12,'box','on')
    DeltaY = 0.1*max([abs(min(min(alpha/(2*pi)))) abs(max(max(alpha/(2*pi))))] );
    %axis([min(tau) max(tau) min(min(alpha/(2*pi))) - DeltaY  max(max(alpha/(2*pi))) + DeltaY ]);

    
    subplot(3,1,2)
    hold on
    plot(tau,alpha_t,'linewidth',1.5)
    plot([tau(1) tau(end)],PD(4)*[1 1]*60/(2*pi),'r-','linewidth',2)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel(alp_tLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
    DeltaY = 0.1*max([abs(min(min(alpha_t))) abs(max(max(alpha_t)))] );
    %axis([min(tau) max(tau) min(min(alpha_t)) - DeltaY  max(max(alpha_t)) + DeltaY ]);
  
    
    subplot(3,1,3)
    hold on
    plot(tau,M_N,'-b')
    plot(tau,M_0,'-- r')
    plot([tau(1) tau(end)]*PD(2),PD(6)*[1 1],'r-','linewidth',2)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel(MLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
    legend('M_N (Nm)','M_0 (Nm)')
    %DeltaY = 0.1*max([abs(min(min(xs(:,end)*PD(5)/PD(1)))) abs(max(max(xs(:,end)*PD(5)/PD(1))))] );
    %axis([min(tau) max(tau) min(min(xs(:,end)*PD(5)/PD(1))) - DeltaY  max(max(xs(:,end)*PD(5)/PD(1))) + DeltaY ]);
  
    
   
    figure(15)
    subplot(3,1,1)
    hold on
    plot(tau,alpha_t,'linewidth',1.5)
    plot([tau(1) tau(end)],PD(4)*[1 1]*60/(2*pi),'r-','linewidth',1.5)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel(alp_tLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
    %ylim([0 120])
   
    
    subplot(3,1,2)
    hold on
    plot(tau,Ten1(1:NP(2):end,:),'linewidth',1.5) 
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    ylabel(FLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
   for i=1:1:length(Ten1(1:NP(2):end,1))
       text = ['Node = ' num2str(i)];
       leg(i,1:length(text)) = text;
   end
   legend(leg)
   %ylim([195 205])
    
    subplot(3,1,3)
    hold on
    plot(tau,M_N,'-b')
    plot(tau,M_0,'-- r')
    %plot([tau(1) tau(end)],PD(6)*[1 1],'r-','linewidth',2)
    xlabel(TLab,'fontsize',12,'interpreter','latex')
    %ylabel(MLab,'fontsize',12,'interpreter','latex')
    set(gca,'fontsize',12,'box','on')
    legend('M_N (Nm)','M_0 (Nm)')
   %ylim([-5 25])
    
   
 

end