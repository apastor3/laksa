function Plot_GroGen_KF(xs_amp,t,rQ,rR,rK,rG,R_KE,rR_Edge,PND,Flag_Dim,PD)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Plot a Ground-Generation system                                %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%     t        -> dimensionless time                                      %%
%     rQ       -> Point Q vector position  (SE components)                %%                                                                         %%
%     rR       -> Rod position vectors (SE components)                    %%
%     rK       -> Kite position vectors (SE components)                   %%                                                                      %%
%     rG       -> Rotor position vectors (SE components)                  %%
%     R_KE     -> SK-SE rotation matrix                                   %%
%     PND      -> dimensionless parameters                                %%
%     Flag_Dim -> 0: inputs and plots are dimensionless                   %%
%                 1: inputs and plot have dimensions                      %%
%     PD       -> Physical parameters                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recover control variables
xc_amp          = Fun_Control_KF(t,xs_amp,PND);


if Flag_Dim == 1
    b   = PD.Tether.L*PND.Kite.b;
    c   = PD.Tether.L*PND.Kite.c;
    h   = PD.Tether.L*PND.Kite.h;
    lt  = PD.AerC.lt;
    hv  = PD.AerC.hv;
    ya  = PD.AerC.ya;
    b_e = sqrt(PD.AerC.St/0.75);
    c_e = 0.75*b_e;
    
    hT_v = 2*PD.AerC.hv;
    c_v  = PD.AerC.Sv/hT_v;
    
    L_blade = PD.Gen.L; 
   
    h0  = PND.Kite.h*PD.Tether.L;
    h0g = h0/2;
else
    b   = PND.Kite.b;
    c   = PND.Kite.c;
    h   = PND.Kite.h;
    lt  = PD.AerC.lt/PD.Tether.L;
    hv  = PD.AerC.hv/PD.Tether.L;
    ya  = PD.AerC.ya/PD.Tether.L;
    
    b_e = sqrt(PD.AerC.St/0.75)/PD.Tether.L;
    c_e = 0.75*b_e;
    
    hT_v = 2*PD.AerC.hv/PD.Tether.L; 
    c_v  = PD.AerC.Sv/(PD.Tether.L^2*hT_v);
    
    L_blade = PD.Gen.L/PD.Tether.L;
    
    h0  = PND.Kite.h;
    h0g = h/2;
end



% Bridle
Q_Plus  = [c/2  b/2 h0g];
Q_Minus = [c/2 -b/2 h0g];
Q_Tail  = [-c/2 0 -(h0-h0g)];
Q_Plus  = rK + R_KE'*Q_Plus';
Q_Minus = rK + R_KE'*Q_Minus';
Q_Tail  = rK + R_KE'*Q_Tail';


% Body frame 
rf1 = rK+R_KE'*[1.1*b/2 0 0]';
rf2 = rK+R_KE'*[0 1.1*b/2 0]';
rf3 = rK+R_KE'*[0 0 1.1*b/2]';

hFig = figure(23);
set(hFig, 'Position', [100 100 1000 600]);
ax1 = axes('Position',[0.08 0.1 0.65 0.7]);
title(ax1,'')
cla(ax1)
hold on

Fac = 5;
% Plot the Kite
Plot_Kite(rK,R_KE,Fac*b,Fac*c,Fac*h0,Fac*h0g)
% Plot point Q, G, rods, and the center of mass of the generators
plot3(-rQ(1),rQ(2),-rQ(3),'*k')
plot3(-rK(1),rK(2),-rK(3),'*k')
plot3(-rR_Edge(1,:),rR_Edge(2,:),-rR_Edge(3,:),'k-')
plot3(-rR(1,:),rR(2,:),-rR(3,:),'ko')
plot3(-[rQ(1) Q_Plus(1)],[rQ(2) Q_Plus(2)],-[rQ(3) Q_Plus(3)],'k-')
plot3(-[rQ(1) Q_Minus(1)],[rQ(2) Q_Minus(2)],-[rQ(3) Q_Minus(3)],'k-')
plot3(-[rQ(1) Q_Tail(1)],[rQ(2) Q_Tail(2)],-[rQ(3) Q_Tail(3)],'k-')

% Plot local horizontal axes
h = plot3([-rK(1)  -rK(1)-1.1*b/2],[rK(2)  rK(2)],[-rK(3)  -rK(3)],'color','red','LineWidth',0.1);
h = plot3([-rK(1)  -rK(1)]    ,[rK(2)  rK(2)+1.1*b/2],[-rK(3)  -rK(3)],'color','red','LineWidth',0.1);
h = plot3([-rK(1)  -rK(1)]    ,[rK(2)  rK(2)],[-rK(3)  -rK(3)-1.1*b/2],'color','red','LineWidth',0.1);

h = plot3([-rK(1)  -rf1(1)],[rK(2)  rf1(2)],[-rK(3)  -rf1(3)],'color','blue','LineWidth',0.1);
h = plot3([-rK(1)  -rf2(1)],[rK(2)  rf2(2)],[-rK(3)  -rf2(3)],'color','blue','LineWidth',0.1);
h = plot3([-rK(1)  -rf3(1)],[rK(2)  rf3(2)],[-rK(3)  -rf3(3)],'color','blue','LineWidth',0.1);


grid on
%view([-45 20])
view([-60 20])
%view([-90 10])
if  Flag_Dim == 1
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    title(ax1,[' t  = ' num2str(t) ' s'])
    axis([0 1.1 -0.75 0.75 0 1.1]*PD.Tether.L)
else
    xlabel('x/L_0')
    ylabel('y/L_0')
    zlabel('z/L_0')
    title([' \tau  = ' num2str(t) ])
    axis([0 1.1 -0.75 0.75 0 1.1])   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot an inset with the Drone %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2 = axes('Position',[0.65 0.65 0.34 0.34]);
cla(ax2)
hold on


% Plot point Q, G, rods, and the center of mass of the generators
plot3(-rQ(1),rQ(2),-rQ(3),'*k')
plot3(-rK(1),rK(2),-rK(3),'*k')
plot3(-[rQ(1) Q_Plus(1)],[rQ(2) Q_Plus(2)],-[rQ(3) Q_Plus(3)],'k-')
plot3(-[rQ(1) Q_Minus(1)],[rQ(2) Q_Minus(2)],-[rQ(3) Q_Minus(3)],'k-')
plot3(-[rQ(1) Q_Tail(1)],[rQ(2) Q_Tail(2)],-[rQ(3) Q_Tail(3)],'k-')

% Plot the Kite
Plot_Kite(rK,R_KE,b,c,h0,h0g)

% Plot local horizontal axes
h = plot3([-rK(1)  -rK(1)-1.1*b/2],[rK(2)  rK(2)],[-rK(3)  -rK(3)],'color','red','LineWidth',0.1);
h = plot3([-rK(1)  -rK(1)]    ,[rK(2)  rK(2)+1.1*b/2],[-rK(3)  -rK(3)],'color','red','LineWidth',0.1);
h = plot3([-rK(1)  -rK(1)]    ,[rK(2)  rK(2)],[-rK(3)  -rK(3)-1.1*b/2],'color','red','LineWidth',0.1);

h = plot3([-rK(1)  -rf1(1)],[rK(2)  rf1(2)],[-rK(3)  -rf1(3)],'color','blue','LineWidth',0.1);
h = plot3([-rK(1)  -rf2(1)],[rK(2)  rf2(2)],[-rK(3)  -rf2(3)],'color','blue','LineWidth',0.1);
h = plot3([-rK(1)  -rf3(1)],[rK(2)  rf3(2)],[-rK(3)  -rf3(3)],'color','blue','LineWidth',0.1);



grid off
view([-25 20])
grid on
%axis equal
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','')
set(gca,'xtick',[],'ytick',[],'ztick',[])

axis([-rK(1)-b -rK(1)+b rK(2)-b rK(2)+b  -rK(3)-b -rK(3)+b]) 

% Set ax1 as current axes
axes(ax1)


end