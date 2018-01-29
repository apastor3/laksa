function U_xs = Grad_U_KF(xs,xc,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Gradient of the potential energy                               %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%      xs and xc  - > State and control vectors                           %%
%      Dimensionless Parameters -> PND                                    %%
% Outputs                                                                 %%
%      Potential gradient -> U_xs                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tether Parameter
sigmaR    = PND.Tether.Sigma;
% Generator Parameters

NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators

if NG>0
    sigmaG    = PND.Gen.Sigma;
else
    sigmaG    = 0;
end

Nvar   = 2*NR+3;
Nvar_p = 2*NR+3+NG;

% Derivative with respect to gamma
U_xs = zeros(Nvar,1);

for j=1:1:NR
   U_xs(j,1) = (1+NG*sigmaG)*xc.lt*cos(xs.gamma(j));
   for k=1:1:j
       if k==j
           g = 1/2;
       else
           g = 1;
       end
       U_xs(k,1) =  U_xs(k,1) +sigmaR*xc.lt^2*g*cos(xs.gamma(k)); 
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contribution from QG  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative with respect to theta
U_xs(2*NR+1,1) = -(1+NG*sigmaG)*xc.lb*(sin(xs.theta)*(cos(xc.delta)*sin(xc.eta)*sin(xs.phi)+sin(xc.delta)*cos(xs.phi))+cos(xs.theta)*cos(xc.delta)*cos(xc.eta));
% Derivative with respect to phi
U_xs(2*NR+3,1) =  (1+NG*sigmaG)*xc.lb*cos(xs.theta)*(cos(xc.delta)*sin(xc.eta)*cos(xs.phi)-sin(xc.delta)*sin(xs.phi));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Contribution from GOG %
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:NG
    % Derivative with respect to theta
    U_xs(2*NR+1,1) = U_xs(2*NR+1,1) + sigmaG*(PND.Gen.x(i)*cos(xs.theta)+PND.Gen.y(i)*sin(xs.theta)*sin(xs.phi)+PND.Gen.z(i)*sin(xs.theta)*cos(xs.phi));
    % Derivative with respect to phi
    U_xs(2*NR+3,1) = U_xs(2*NR+3,1) - sigmaG*(PND.Gen.y(i)*cos(xs.theta)*cos(xs.phi)-PND.Gen.z(i)*cos(xs.theta)*sin(xs.phi));
end

end
