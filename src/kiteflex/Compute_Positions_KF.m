function [rQ rR rR_Edge rK rG ]     = Compute_Positions_KF(xs,xc,R_KE,PND)
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute position vector of the rigid bodies                    %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %%
%     State variables  -> xs -> gamma,varphi,theta,psi,phi                %%
%     Control variable -> xc -> lt,lb,delta,eta                           %%
%     Rotation Matrix  -> R_KE                                            %%
%     PND              -> Dimensionless parameters                        %%
% Output: All vector components are in the Earth frame                    %%
%    rQ = position of point Q                                             %%
%    rR = position vector of the center of mass of the rods               %%
%    rK = position of the center of mass of the kite                      %%
%    rG = position vector of the generators                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      Point Q    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
rQ = zeros(3,1);
for j=1:1:NR
  rQ(1,1) = rQ(1,1)-xc.lt*cos(xs.gamma(j))*cos(xs.varphi(j));
  rQ(2,1) = rQ(2,1)-xc.lt*cos(xs.gamma(j))*sin(xs.varphi(j));
  rQ(3,1) = rQ(3,1)-xc.lt*sin(xs.gamma(j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      Tether     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
rR      = zeros(3,NR);
rR_Edge = zeros(3,NR+1); 
for i=1:1:NR
    for j=1:1:i
        if j==i
            pj = 0.5;
        else
            pj = 1;
        end
        rR(1,i) = rR(1,i)-pj*xc.lt*cos(xs.gamma(j))*cos(xs.varphi(j));
        rR(2,i) = rR(2,i)-pj*xc.lt*cos(xs.gamma(j))*sin(xs.varphi(j));
        rR(3,i) = rR(3,i)-pj*xc.lt*sin(xs.gamma(j));
        
        rR_Edge(1,i+1) = rR_Edge(1,i)-xc.lt*cos(xs.gamma(j))*cos(xs.varphi(j));
        rR_Edge(2,i+1) = rR_Edge(2,i)-xc.lt*cos(xs.gamma(j))*sin(xs.varphi(j));
        rR_Edge(3,i+1) = rR_Edge(3,i)-xc.lt*sin(xs.gamma(j));
    end
end


%%%%%%%%%%%%%%%%%%%%
%%%%    Kite    %%%%
%%%%%%%%%%%%%%%%%%%%
% QG components in the body frame
QG = -xc.lb*[cos(xc.delta)*cos(xc.eta) cos(xc.delta)*sin(xc.eta) sin(xc.delta)]';
% QG components in the Earth frame
QG = R_KE'*QG;
% Vector position of the center of mass of the Kite
rK = rQ+QG;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Generators   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
rG = zeros(3,NG);
for i=1:1:NG
    GOG     = [ PND.Gen.x(i)  PND.Gen.y(i)  PND.Gen.z(i)]';
    rG(:,i) = rK + R_KE'*GOG; 
end