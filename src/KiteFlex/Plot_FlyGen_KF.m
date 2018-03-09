function Plot_FlyGen_KF(xs_amp,t,rQ,rR,rK,rG,R_KE,rR_Edge,PND,Flag_Dim,PD)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Plot a Fly-Generation System                                   %
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
    for i=1:1:PD.Gen.Num
        xG(i) = PD.Gen.x(i);
        yG(i) = PD.Gen.y(i);
        zG(i) = PD.Gen.z(i);
    end
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
    for i=1:1:PD.Gen.Num
        xG(i) = PD.Gen.x(i)/PD.Tether.L;
        yG(i) = PD.Gen.y(i)/PD.Tether.L;
        zG(i) = PD.Gen.z(i)/PD.Tether.L;
    end
end

% Do some calculations

% Bridle
Q_Plus  = [0 b/4 0];
Q_Minus = [0 -b/4 0];
Q_Tail  = [-lt 0 0];
Q_Plus  = rK + R_KE'*Q_Plus';
Q_Minus = rK + R_KE'*Q_Minus';
Q_Tail  = rK + R_KE'*Q_Tail';


% Rotors
angle = [0:30:360 0]*pi/180;
for j=1:1:length(angle)
    Blade_Body(j,:) = L_blade*[0 cos(angle(j)) sin(angle(j))];
end
for i=1:1:PND.Gen.Num 
    Support_Body(1,1:3,i) = [ xG(i) yG(i) zG(i)];
    Support_Body(2,1:3,i) = [ xG(i) yG(i) 0];
    Support_Body(3,1:3,i) = [ 0     yG(i) 0];
end

% Plot the wing of the drone
x_wing   = c/2*[-1  1];
y_wing   = b/2*[-1  1];

Wing_Body = [-c/2 -b/2 0;...
              c/2 -b/2 0;...
              c/2  b/2 0;...
             -c/2  b/2 0];
         
Elevator_Body = [-lt-c_e/2 -b_e/2 0;...
                 -lt+c_e/2 -b_e/2 0;...
                 -lt+c_e/2  b_e/2 0;...
                 -lt-c_e/2  b_e/2 0];         

Fuselage_Body = [0.1*lt+c/2  0  0.1*hv;...
                 -lt     0  0.1*hv;...
                 -lt     0 -0.1*hv;...
                 0.1*lt+c/2  0 -0.1*hv;];
             
Tail_Body        = [-lt+c_v/2  0  0.1*hT_v;...
                    -lt-c_v/2  0  0.1*hT_v;...
                    -lt-c_v/2  0 -0.1*hT_v;...
                    -lt+c_v/2  0 -0.1*hT_v;];             
                
for i = 1:1:4
    Wing_Earth(i,:)     = rK + R_KE'*Wing_Body(i,:)';
    Fuselage_Earth(i,:) = rK + R_KE'*Fuselage_Body(i,:)';
    Elevator_Earth(i,:) = rK + R_KE'*Elevator_Body(i,:)';
    Tail_Earth(i,:)     = rK + R_KE'*Tail_Body(i,:)';
end
for i=1:1:3
    for j=1:1:PND.Gen.Num 
        Support_Earth(i,:,j) = rK + R_KE'*Support_Body(i,:,j)';
    end
end
for i=1:1:length(angle)
    for j=1:1:PND.Gen.Num 
        Blade_Earth(i,:,j) = rG(:,j) + R_KE'*Blade_Body(i,:)';
    end
end

% Body axis
rf1 = rK+R_KE'*[1.1*b/2 0 0]';
rf2 = rK+R_KE'*[0 1.1*b/2 0]';
rf3 = rK+R_KE'*[0 0 1.1*b/2]';

hFig = figure(23);
set(hFig, 'Position', [100 100 1000 600]);
ax1 = axes('Position',[0.08 0.1 0.65 0.7]);
title(ax1,'')
cla(ax1)
hold on

% Plot point Q, G, rods, and the center of mass of the generators
plot3(-rQ(1),rQ(2),-rQ(3),'*k')
plot3(-rK(1),rK(2),-rK(3),'*k')
plot3(-rR_Edge(1,:),rR_Edge(2,:),-rR_Edge(3,:),'k-')
plot3(-rR(1,:),rR(2,:),-rR(3,:),'ko')
plot3(-[rQ(1) Q_Plus(1)],[rQ(2) Q_Plus(2)],-[rQ(3) Q_Plus(3)],'k-')
plot3(-[rQ(1) Q_Minus(1)],[rQ(2) Q_Minus(2)],-[rQ(3) Q_Minus(3)],'k-')
plot3(-[rQ(1) Q_Tail(1)],[rQ(2) Q_Tail(2)],-[rQ(3) Q_Tail(3)],'k-')

for j=1:1:PND.Gen.Num 
  patch(-Support_Earth(:,1,j),Support_Earth(:,2,j),-Support_Earth(:,3,j),'.r')
  plot3(-Blade_Earth(:,1,j), Blade_Earth(:,2,j),-Blade_Earth(:,3,j),'k' )
end

patch(-Wing_Earth(:,1),Wing_Earth(:,2),-Wing_Earth(:,3),'.b')
patch(-Fuselage_Earth(:,1),Fuselage_Earth(:,2),-Fuselage_Earth(:,3),'.b')
patch(-Elevator_Earth(:,1),Elevator_Earth(:,2),-Elevator_Earth(:,3),'.r')
patch(-Tail_Earth(:,1),Tail_Earth(:,2),-Tail_Earth(:,3),'.r')

% Plot local horizontal axes
h = plot3([-rK(1)  -rK(1)-1.1*b/2],[rK(2)  rK(2)],[-rK(3)  -rK(3)],'color','red','LineWidth',0.1);
h = plot3([-rK(1)  -rK(1)]    ,[rK(2)  rK(2)+1.1*b/2],[-rK(3)  -rK(3)],'color','red','LineWidth',0.1);
h = plot3([-rK(1)  -rK(1)]    ,[rK(2)  rK(2)],[-rK(3)  -rK(3)-1.1*b/2],'color','red','LineWidth',0.1);

h = plot3([-rK(1)  -rf1(1)],[rK(2)  rf1(2)],[-rK(3)  -rf1(3)],'color','blue','LineWidth',0.1);
h = plot3([-rK(1)  -rf2(1)],[rK(2)  rf2(2)],[-rK(3)  -rf2(3)],'color','blue','LineWidth',0.1);
h = plot3([-rK(1)  -rf3(1)],[rK(2)  rf3(2)],[-rK(3)  -rf3(3)],'color','blue','LineWidth',0.1);


grid on
%view([-45 20])
view([-90 0])

if  Flag_Dim == 1
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    title(ax1,[' t  = ' num2str(t) ' s'])
    axis([0 1.2 -0.6 0.6 0 1.2]*PD.Tether.L)
else
    xlabel('x/L_0')
    ylabel('y/L_0')
    zlabel('z/L_0')
    title([' \tau  = ' num2str(t) ])
    axis([0 1.2 -0.6 0.6 0 1.2])   
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

for j=1:1:PND.Gen.Num 
  patch(-Support_Earth(:,1,j),Support_Earth(:,2,j),-Support_Earth(:,3,j),'.r')
  plot3(-Blade_Earth(:,1,j), Blade_Earth(:,2,j),-Blade_Earth(:,3,j),'k') 
end

patch(-Wing_Earth(:,1),Wing_Earth(:,2),-Wing_Earth(:,3),'.b')
patch(-Fuselage_Earth(:,1),Fuselage_Earth(:,2),-Fuselage_Earth(:,3),'.b')
patch(-Elevator_Earth(:,1),Elevator_Earth(:,2),-Elevator_Earth(:,3),'.r')
patch(-Tail_Earth(:,1),Tail_Earth(:,2),-Tail_Earth(:,3),'.r')

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