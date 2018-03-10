function Plot_KS(X,t,rk,vk,ak,RBE,R2E,R3E,Rp,Rm,PND,Flag_Dim,PD)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Author    : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Plot the mechanical system                                     %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %%
% Inputs: t   -> dimensionless time                                    %%
%         rk  -> center of mass position (SE components)               %%
%         RBE -> Body-Earth rotation matrix                            %%
%         PND -> Dimensionless parameter                               %%
%         Flag_Dim  -> Flag controlling the units of the inputs        %%
%         PD  -> Physical Parameters                                   %%
%                                                                      %%
% Outputs: Plot the KiteAcrobat System                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontsi = 12;
Xmax = 1.0;
Ymax = 1.0;
Zmax = 1.1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover control variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xc xc_p xc_pp] = Fun_Control_KS(t,X,PND,rk,vk,ak);
PR              = xc(1,1);
lambda          = xc(2,1);

d_cb            = (1-PR)*(PND.Bar.Ls-PND.Bar.Lds-PND.Bar.Lps);


hFig = figure(23);
set(hFig, 'Position', [100 100 1000 600]);
ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
axes(ax1)
title(ax1,'')
cla(ax1)
hold on

if Flag_Dim == 1
    b   = PD.Tether.Ll*PND.Kite.b;
    c   = PD.Tether.Ll*PND.Kite.c;
    h   = PD.Tether.Ll*PND.Kite.h;
    hg  = PD.Tether.Ll*PND.Kite.hg;
    xLE = PD.Tether.Ll*abs(PND.Tether.XA);
    xTE = PD.Tether.Ll*abs(PND.Tether.XB);
    yLE1 = abs(PD.Tether.YA);
    zLE1 = PD.Tether.ZA;
else
    b   = PND.Kite.b;
    c   = PND.Kite.c;
    h   = PND.Kite.h;
    hg  = PND.Kite.hg;
    xLE = abs(PND.Tether.XA);
    xTE = abs(PND.Tether.XB);
    yLE1 = abs(PND.Tether.YA);
    zLE1 = PND.Tether.ZA;
end

% Plot the Kite
Plot_Kite_KS(rk,RBE,b,xLE,xTE,yLE1,zLE1,h,hg);
h1 = gca;
% Plot the tethers
if  Flag_Dim == 1
    GAP    =  RBE'*PD.Tether.Ll*[PND.Tether.XA    PND.Tether.YA   PND.Tether.ZA]';   % Earth Components
    GAM    =  RBE'*PD.Tether.Ll*[PND.Tether.XA   -PND.Tether.YA   PND.Tether.ZA]';   % Earth Components
    
    OEF_OE = -R3E'*PD.Tether.Ll*[0 0 PND.Bar.Ls]';               % Earth Components
    OEF_C0 = -R3E'*PD.Tether.Ll*[0 0 PND.Bar.Lps+d_cb]';         % Earth Components
    
    OE_O2  = -R2E'*PD.Tether.Ll*[0 0 PND.Tether.l]';
    
    C0_CP  =  R3E'*PD.Tether.Ll*PND.Bar.Lc/2*[0 cos(lambda) sin(lambda)]';    % Earth Components
    C0_CM  = -R3E'*PD.Tether.Ll*PND.Bar.Lc/2*[0 cos(lambda) sin(lambda)]';    % Earth Components
    OEF_PS = -R3E'*PD.Tether.Ll*[0 0 PND.Bar.Lps]';              % Earth Components
    OEF_DS = -R3E'*PD.Tether.Ll*[0 0 PND.Bar.Ls-PND.Bar.Lds]';   % Earth Components
else
    GAP    =  RBE'*[PND.Tether.XA    PND.Tether.YA   PND.Tether.ZA]';  % Earth Components
    GAM    =  RBE'*[PND.Tether.XA   -PND.Tether.YA   PND.Tether.ZA]';  % Earth Components
    
    OEF_OE = -R3E'*[0 0 PND.Bar.Ls]';               % Earth Components
    OEF_C0 = -R3E'*[0 0 PND.Bar.Lps+d_cb]';         % Earth Components
  
    OE_O2  = -R2E'*[0 0 PND.Tether.l]';
  
    C0_CP  =  R3E'*PND.Bar.Lc/2*[0 cos(lambda) sin(lambda)]';    % Earth Components
    C0_CM  = -R3E'*PND.Bar.Lc/2*[0 cos(lambda) sin(lambda)]';    % Earth Components
    OEF_PS = -R3E'*[0 0 PND.Bar.Lps]';              % Earth Components
    OEF_DS = -R3E'*[0 0 PND.Bar.Ls-PND.Bar.Lds]';   % Earth Components
    
end
 
OEF_AP  =  rk + GAP;
OEF_AM  =  rk + GAM;

OEF_O2  = OEF_OE+OE_O2;

plot3([-OEF_OE(1) -OEF_AP(1)],[OEF_OE(2) OEF_AP(2)],[-OEF_OE(3) -OEF_AP(3)],'k')
plot3([-OEF_OE(1) -OEF_AM(1)],[OEF_OE(2) OEF_AM(2)],[-OEF_OE(3) -OEF_AM(3)],'k')

%% Plot the Bar
plot3([0 -OEF_OE(1)],[0 OEF_OE(2)],[0  -OEF_OE(3)],'b','linewidth',2 )
plot3(-OEF_PS(1),OEF_PS(2),-OEF_PS(3),'or','MarkerFaceColor',[1 0 0],'MarkerSize',4)
plot3(-OEF_DS(1),OEF_DS(2),-OEF_DS(3),'or','MarkerFaceColor',[1 0 0],'MarkerSize',4)
plot3([-OEF_C0(1) -(OEF_C0(1)+C0_CP(1))],[OEF_C0(2) (OEF_C0(2)+C0_CP(2))],[-OEF_C0(3) -(OEF_C0(3)+C0_CP(3))],'g','linewidth',2)
plot3([-OEF_C0(1) -(OEF_C0(1)+C0_CM(1))],[OEF_C0(2) (OEF_C0(2)+C0_CM(2))],[-OEF_C0(3) -(OEF_C0(3)+C0_CM(3))],'g','linewidth',2)

% Plot the Flexible tethers
plot3(-Rp(:,1),Rp(:,2),-Rp(:,3),'k')
plot3(-Rm(:,1),Rm(:,2),-Rm(:,3),'k')


grid on
view([-45 20])
if  Flag_Dim == 1
    xlabel('x (m)','fontsize',fontsi)
    ylabel('y (m)','fontsize',fontsi)
    zlabel('z (m)','fontsize',fontsi)
    title(ax1,[' t  = ' num2str(t*sqrt(PD.Tether.Ll/PD.Env.g)) ' s'],'fontsize',fontsi)  
    axis([0 Xmax -Ymax Ymax 0 Zmax]*PD.Tether.Ll)
else
    xlabel('x/L_l','fontsize',fontsi)
    ylabel('y/L_l','fontsize',fontsi)
    zlabel('z/L_l','fontsize',fontsi)
    title([' \tau  = ' num2str(t) ])
    axis([0 Xmax -Ymax Ymax 0 Zmax])
end
set(gca,'fontsize',fontsi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot an inset with the kite  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2 = axes('Position',[0.65 0.7 0.28 0.28]);
cla(ax2)
hold on
Plot_Kite_KS(rk,RBE,b,xLE,xTE,yLE1,zLE1,h,hg);
% Plot the inelastic tethers
plot3([-OEF_OE(1) -OEF_AP(1)],[OEF_OE(2) OEF_AP(2)],[-OEF_OE(3) -OEF_AP(3)],'k')
plot3([-OEF_OE(1) -OEF_AM(1)],[OEF_OE(2) OEF_AM(2)],[-OEF_OE(3) -OEF_AM(3)],'k')
% Plot the Flexible tethers
plot3(-Rp(:,1),Rp(:,2),-Rp(:,3),'k')
plot3(-Rm(:,1),Rm(:,2),-Rm(:,3),'k')


grid on
view([-45 20])
axis equal
axis([-rk(1)-1.1*b  -rk(1)+1.1*b  rk(2)-1.5*b rk(2)+1.5*b  -rk(3)-h  -rk(3)+h]  )
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot an inset with the control bar  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3 = axes('Position',[0.05 0.78 0.25 0.2]);
cla(ax3)
hold on
%% Plot the Bar in S3 axes
OEF_OE = R2E*OEF_OE; 
OEF_PS = R2E*OEF_PS; 
OEF_DS = R2E*OEF_DS; 
OEF_C0 = R2E*OEF_C0;
C0_CP  = R2E*C0_CP;
C0_CM  = R2E*C0_CM;
OEF_AP = R2E*OEF_AP;
OEF_AM = R2E*OEF_AM;
OEF_O2 = R2E*OEF_O2;


plot(-[0 OEF_OE(2)],[0  -OEF_OE(3)],'b','linewidth',2 )
plot(-OEF_PS(2),-OEF_PS(3),'or','MarkerFaceColor',[1 0 0],'MarkerSize',4)
plot(-OEF_DS(2),-OEF_DS(3),'or','MarkerFaceColor',[1 0 0],'MarkerSize',4)
plot(-[OEF_C0(2) (OEF_C0(2)+C0_CP(2))],[-OEF_C0(3) -(OEF_C0(3)+C0_CP(3))],'g','linewidth',2)
plot(-[OEF_C0(2) (OEF_C0(2)+C0_CM(2))],[-OEF_C0(3) -(OEF_C0(3)+C0_CM(3))],'g','linewidth',2)

plot(-(OEF_C0(2)+C0_CP(2)),-(OEF_C0(3)+C0_CP(3)),'ob','MarkerFaceColor',[0 0 1],'MarkerSize',4)
plot(-(OEF_C0(2)+C0_CM(2)), -(OEF_C0(3)+C0_CM(3)),'ok','MarkerFaceColor',[0 0 0],'MarkerSize',4)


plot(-[OEF_C0(2)++C0_CP(2)  OEF_C0(2)+C0_CM(2)],[-OEF_C0(3) -OEF_C0(3)],'k','linewidth',.25)


% Plot the inelastic tethers
plot(-[OEF_OE(2) OEF_AP(2)],-[OEF_OE(3) OEF_AP(3)],'k')
plot(-[OEF_OE(2) OEF_AM(2)],-[OEF_OE(3) OEF_AM(3)],'k')
plot(-[OEF_OE(2) OEF_O2(2)],-[OEF_OE(3) OEF_O2(3)],'--k')


xticklabels({''})
yticklabels({''})
grid on
set(gca,'box','on','fontsize',12)

if Flag_Dim==0
    axis([-1.3*PND.Bar.Lc/2 1.3*PND.Bar.Lc/2 0 2*PND.Bar.Ls]);
    text(-0.5*PND.Bar.Lc/2,0.2*PND.Bar.Ls,'\Pi')
else
    axis([-1.3*PD.Bar.Lc/2  1.3*PD.Bar.Lc/2 0 2*PD.Bar.Ls]);
    text(-0.5*PD.Bar.Lc/2,0.2*PD.Bar.Ls,'\Pi')
end


% Set ax1 as current axes
axes(ax1)

end