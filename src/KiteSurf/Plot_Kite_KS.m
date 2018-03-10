function Plot_Kite_KS(rk,RBE,b,xLE,xTE,yLE1,zLE1,h,hg)
     
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Author    : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Plot the kite                                                  %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%----------------------------------------------------------------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %%
% Inputs: rk   -> center of mass position (SE components)       %%
%         RBE  -> Body-Earth rotation matrix                    %%
%         b    -> Kite span                                     %%
%         xLE  -> x-coordinate of the leading edge              %%
%         xTE  -> x-coordinate of the trailing edge             %%
%         yLE1 -> y-coordinate of A                             %%
%         zLE1 -> z-coordinate of A                             %% 
%         h    -> Kite height                                   %%
%         hg   -> Height of the center of mass with respect to  %%
%                kite tips                                      %% 
% Outputs: Plot the kite, center of mass, horizontal and body   %% 
%          axes in the current figure                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Center of mass position
plot3(-rk(1), rk(2),-rk(3),'o','LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);

% Find the ellipse in 2D with simple axis (chi-eta)
chi = [-20:1:20]*b/(2*20);
eta1 = real(sqrt((h)^2.*(1-(2*chi./b).^2)));

eta2 = interp1([-b/2 -yLE1 0 yLE1 b/2],[0 hg-zLE1 h hg-zLE1 0 ],chi,'spline');

% Find the ellipse in body axes (displace h_g distance and create a third component)
x_ellip_body_LE =  xLE*ones(1,length(chi));
y_ellip_body_LE = -chi;
z_ellip_body_LE = -(eta2-hg);
ellipse_1_Body = [x_ellip_body_LE' y_ellip_body_LE' z_ellip_body_LE'];

x_ellip_body_TE = -xTE*ones(1,length(chi));
y_ellip_body_TE = -chi;
z_ellip_body_TE = -(eta1-hg);
ellipse_2_Body = [x_ellip_body_TE' y_ellip_body_TE' z_ellip_body_TE'];


% Find the ellipse in Earth Component
for i = 1:length(ellipse_1_Body)
    
    ellipse_1_Earth(i,:) = rk + RBE'*ellipse_1_Body(i,:)';
    ellipse_2_Earth(i,:) = rk + RBE'*ellipse_2_Body(i,:)';
    % Plot of the lines joining both ellipses
    plot3(-[ellipse_2_Earth(i,1) ellipse_1_Earth(i,1)],[ellipse_2_Earth(i,2) ellipse_1_Earth(i,2)],-[ellipse_2_Earth(i,3) ellipse_1_Earth(i,3)],'b','linewidth',1.5)
end

% Plot the surface of the Kite
for i = 2:length(ellipse_1_Earth)
    
    X = [ellipse_1_Earth(i,1) ellipse_1_Earth(i-1,1) ellipse_2_Earth(i-1,1) ellipse_2_Earth(i,1)];
    Y =  [ellipse_1_Earth(i,2) ellipse_1_Earth(i-1,2) ellipse_2_Earth(i-1,2) ellipse_2_Earth(i,2)];
    Z = [ellipse_1_Earth(i,3) ellipse_1_Earth(i-1,3) ellipse_2_Earth(i-1,3) ellipse_2_Earth(i,3)];
    patch(-X,Y,-Z,'c')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Plot The Kite
%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot3(-rk(1),rk(2),-rk(3),'o','LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);

% Plot local horizontal axes
plot3([-rk(1)  -rk(1)-1.1*b/2],[rk(2)  rk(2)],[-rk(3)  -rk(3)],'color','red','LineWidth',0.1);
plot3([-rk(1)  -rk(1)]    ,[rk(2)  rk(2)+1.1*b/2],[-rk(3)  -rk(3)],'color','red','LineWidth',0.1);
plot3([-rk(1)  -rk(1)]    ,[rk(2)  rk(2)],[-rk(3)  -rk(3)-1.1*b/2],'color','red','LineWidth',0.1);

% Plot Body axis
rf1 = rk+RBE'*[1.1*b/2 0 0]';
rf2 = rk+RBE'*[0 1.1*b/2 0]';
rf3 = rk+RBE'*[0 0 1.1*b/2]';

plot3([-rk(1)  -rf1(1)],[rk(2)  rf1(2)],[-rk(3)  -rf1(3)],'color','blue','LineWidth',0.1);
plot3([-rk(1)  -rf2(1)],[rk(2)  rf2(2)],[-rk(3)  -rf2(3)],'color','blue','LineWidth',0.1);
plot3([-rk(1)  -rf3(1)],[rk(2)  rf3(2)],[-rk(3)  -rf3(3)],'color','blue','LineWidth',0.1);


end