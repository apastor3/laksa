function Plot_KA(t,rk,RBE,PND,Flag_Dim,PD)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Plot the mechanical system                                     %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %%
% Inputs: t   -> dimensionless time                             %%
%         rk  -> center of mass position (SE components)        %%
%         RBE -> Body-Earth rotation matrix                     %%
%         PND -> Dimensionless parameter                        %%
%         Flag_Dim  -> Flag controlling the units of the inputs %%
%         PD  -> Physical Parameters                            %%
%                                                               %%
% Outputs: Plot the KiteAcrobat System                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hFig = figure(23);
delete(hFig.Children);
set(hFig, 'Position', [100 100 1000 600]);
ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
title(ax1,'')
cla(ax1)
hold on

if Flag_Dim == 1
    b  = PD.Tether.L0*PND.Kite.b;
    c  = PD.Tether.L0*PND.Kite.c;
    h  = PD.Tether.L0*PND.Kite.h;
    hg = PD.Tether.L0*PND.Kite.hg;
else
    b  = PND.Kite.b;
    c  = PND.Kite.c;
    h  = PND.Kite.h;
    hg = PND.Kite.hg;
end

% Plot the Kite

Plot_Kite(rk,RBE,b,c,h,hg);
h1 = gca;
% Plot the tethers
if  Flag_Dim == 1
    GAP = PD.Tether.L0*[PND.Tether.XA    PND.Tether.YA   PND.Tether.ZA]';
    GAM = PD.Tether.L0*[PND.Tether.XA   -PND.Tether.YA   PND.Tether.ZA]';
else
    GAP = [PND.Tether.XA    PND.Tether.YA   PND.Tether.ZA]';
    GAM = [PND.Tether.XA   -PND.Tether.YA   PND.Tether.ZA]';
end
OEAP = rk+RBE'*GAP;
OEAM = rk+RBE'*GAM;


plot3([0 -OEAP(1)],[0 OEAP(2)],[0 -OEAP(3)],'k')
plot3([0 -OEAM(1)],[0 OEAM(2)],[0 -OEAM(3)],'k')

%axis equal

grid on
view([-45 20])
if  Flag_Dim == 1
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    title(ax1,[' \tau  = ' num2str(t*sqrt(PD.Tether.L0/PD.Env.g)) ])
    axis([0 1. -0.5 0.5 0 1.1]*PD.Tether.L0)
else
    xlabel('x/L_0')
    ylabel('y/L_0')
    zlabel('z/L_0')
    title([' \tau  = ' num2str(t) ])
     axis([0 1. -0.5 0.5 0 1.1])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot an inset with the kite  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2 = axes('Position',[0.65 0.7 0.28 0.28]);
cla(ax2)
hold on
Plot_Kite(rk,RBE,b,c,h,hg);
% Plot the tether
plot3([0 -OEAP(1)],[0 OEAP(2)],[0 -OEAP(3)],'k')
plot3([0 -OEAM(1)],[0 OEAM(2)],[0 -OEAM(3)],'k')
grid on
view([-45 20])
axis equal
axis([-rk(1)-1.1*b  -rk(1)+1.1*b  rk(2)-1.5*b rk(2)+1.5*b  -rk(3)-h  -rk(3)+h]  )
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','')

% Set ax1 as current axes
axes(ax1)

end
