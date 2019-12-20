function Plot_KT(T_out,rk,RBE,PND,PD,Flag_Dim,ax1)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Jose Antonio Serrano-Iglesias         %
% Language  : Matlab                                                         %
% Synopsis  : Plot the mechanical system                                     %
% Copyright : Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %%
% Inputs: T_out     -> dimensionless time                       %%
%         rk        -> center of mass position (SE components)  %%
%         RBE       -> Body-Earth rotation matrix               %%
%         PND       -> Dimensionless parameter                  %%
%         Flag_Dim  -> Flag controlling the units of the inputs %%
%         PD        -> Physical Parameters                      %%
%                                                               %%
% Outputs: Plot the KiteTrain System                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:1:PND.Kite.N 
 if Flag_Dim == 1
    b  = PD.Tether.L(i)*PND.Kite.b(i);
    c  = PD.Tether.L(i)*PND.Kite.c(i);
    h  = PD.Tether.L(i)*PND.Kite.h(i);
    hg = PD.Tether.L(i)*PND.Kite.hg(i);
 else
    b  = PND.Kite.b(i);
    c  = PND.Kite.c(i);
    h  = PND.Kite.h(i);
    hg = PND.Kite.hg(i);
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%  Plot the kite  %%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 Plot_Kite(rk(:,i),RBE(:,:,i),b,c,h,hg);
 h1 = gca;
 
 if  Flag_Dim == 1
   GAP = PD.Tether.L(i)*[PND.Tether.XA(i)    PND.Tether.YA(i)   PND.Tether.ZA(i)]';
   GAM = PD.Tether.L(i)*[PND.Tether.XA(i)   -PND.Tether.YA(i)   PND.Tether.ZA(i)]';
   
   GCP = PD.Tether.L(i)*[PND.Tether.XC(i)    PND.Tether.YC(i)   PND.Tether.ZC(i)]';
   GCM = PD.Tether.L(i)*[PND.Tether.XC(i)   -PND.Tether.YC(i)   PND.Tether.ZC(i)]';
 else
   GAP = [PND.Tether.XA(i)    PND.Tether.YA(i)   PND.Tether.ZA(i)]';
   GAM = [PND.Tether.XA(i)   -PND.Tether.YA(i)   PND.Tether.ZA(i)]';
   
   GCP = [PND.Tether.XC(i)    PND.Tether.YC(i)   PND.Tether.ZC(i)]';
   GCM = [PND.Tether.XC(i)   -PND.Tether.YC(i)   PND.Tether.ZC(i)]';
 end
 
  OEAP(:,i) = rk(:,i)+RBE(:,:,i)'*GAP;
  OEAM(:,i) = rk(:,i)+RBE(:,:,i)'*GAM;
   
  OECP(:,i) = rk(:,i)+RBE(:,:,i)'*GCP;
  OECM(:,i) = rk(:,i)+RBE(:,:,i)'*GCM;
   
  if i == 1
   % Plot the tethers of the first kite %
   plot3([0 -OEAP(1,i)],[0 OEAP(2,i)],[0 -OEAP(3,i)],'k')
   plot3([0 -OEAM(1,i)],[0 OEAM(2,i)],[0 -OEAM(3,i)],'k')

  else
   % Plot the tethers of the rest of kites %
   plot3([-OECP(1,i-1) -OEAP(1,i)],[OECP(2,i-1) OEAP(2,i)],[-OECP(3,i-1) -OEAP(3,i)],'k')
   plot3([-OECM(1,i-1) -OEAM(1,i)],[OECM(2,i-1) OEAM(2,i)],[-OECM(3,i-1) -OEAM(3,i)],'k')
 
  end
 
 
end



 grid on
 view([-45 20])
 axis equal
 
 if Flag_Dim == 1
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    title(ax1,[' time  = ' num2str(T_out*sqrt(PD.Tether.L(1)/PD.Env.g)) ' s '])
    axis([0 1. -1.0 1.0 0 1.1]*PD.Tether.L(1)*PD.Kite.N)
   
 else
    xlabel('x/L_1')
    ylabel('y/L_1')
    zlabel('z/L_1')
    title([' \tau  = ' num2str(T_out) ])
    axis([0 1. -0.5 0.5 0 1.1])
 end


end


