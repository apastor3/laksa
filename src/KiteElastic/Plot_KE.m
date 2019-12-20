function Plot_KE(t,X,PND,ax1)

%%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute all relevant quantities                                %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %%
% Inputs: t   -> dimensionless time                                         %%
%         X   -> extended state vector                                      %%
%         PND -> Dimensionless parameter                                    %%
% Outputs: Plot the system                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:PND.Kite.Num
    % Kinemtic quantities
    rk(:,i)     = X(12*(i-1)+1:12*(i-1)+3,1);   % Position vector (SE components)
    eulerK(:,i) = X(12*(i-1)+7:12*(i-1)+9,1);   % Euler angles  
    % Rotation matrix
    R_KE(:,:,i) = Matrix_R_KE_KE(eulerK(:,i));
end
% Recover Mass positions
Ind  = 12*PND.Kite.Num;
for i=1:1:PND.Tether.Num 
    for j=1:1:PND.Tether.Mass(i)
        rM(:,i,j)   = X(Ind+1:Ind+3,1); % Masses position vectors
        Ind         = Ind+3;
    end
end



% Plot the Kites
for i=1:1:PND.Kite.Num
    b    = PND.Kite.b(i);
    c    = PND.Kite.c(i);
    h    = PND.Kite.h(i);
    hg   = PND.Kite.hg(i);
   % xLE  = PND.Kite.c(i)/2;
   % xTE  = PND.Kite.c(i)/2(i);
   % yLE1 = PND.Kite.b(i)/2;
   % zLE1 = PND.Kite.h(i)-PND.Kite.hg(i);
    
     Plot_Kite(rk(:,i),R_KE(:,:,i),b,c,h,hg);
end
h1 = gca;

for i=1:1:PND.Tether.Num  % Do all the tethers
    %Position vector of the Lower point  (SE components)
    r = [];
    Kite_Down = PND.Tether.Down(i);
    if Kite_Down==0
        r = [0 0 0]';
    else
       GD  = [PND.Tether.Dx(i)  PND.Tether.Dy(i)   PND.Tether.Dz(i)]';    % SB components
       r   = rk(:,Kite_Down) + squeeze(R_KE(:,:,Kite_Down))'*GD;          % SE components
    end
    for j=1:1:PND.Tether.Mass(i)
        r = [r rM(:,i,j)];
    end  
    %Position vector of the Upper point  (SE components)
    Kite_Up   = PND.Tether.Up(i); 
    GU        = [PND.Tether.Ux(i)  PND.Tether.Uy(i)   PND.Tether.Uz(i)]';    % SB components
    rU        = rk(:,Kite_Up) + squeeze(R_KE(:,:,Kite_Up))'*GU;              % SE components
    r         = [r rU];
    plot3(-r(1,:),r(2,:),-r(3,:),'k-')
    
    plot3(-r(1,:),r(2,:),-r(3,:),'k.','markersize',10)
end



%title(['Time = ' num2str(t)])
grid on
view([-45 20])
%axis equal


% Set ax1 as current axes
axes(ax1)
 pause(0.001)
end