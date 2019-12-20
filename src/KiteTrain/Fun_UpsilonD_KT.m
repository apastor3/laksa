function [UpsD UpsD_xs]= Fun_UpsilonD_KT(xs,Grad_ZX,Grad2_ZX,Kite,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : SB-S2 rotation matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %%
% Inputs:  xs       -> state vector of the system                    %%
%          Grad_d   -> Gradient of d                                 %% 
%          Grad2_d  -> xs-derivatives of gradient of d               %%
%          Kite     -> Kite Number
%          PND      -> Dimensionless parameters                      %%
%                                                                    %%
% Outputs: UpsD    -> Matrix Upsilon^D                               %%
%          UpsD_xs -> xs-Derivative of UpsD                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the matrix
UpsD    = zeros(3,4*PND.Kite.N);
UpsD_xs = zeros(3,4*PND.Kite.N,4*PND.Kite.N);
Vec1    = zeros(3,1);
Vec1_xs = zeros(3,4);
Vec2    = zeros(3,1);
Vec2_xs = zeros(3,4);
% Recover the variables
varphi = xs(1,1);
gamma  = xs(2,1);
eta    = xs(3,1);

Grad_zeta  = squeeze(Grad_ZX(1,:));
Grad2_zeta = squeeze(Grad2_ZX(1,:,:));

Grad_xi    = squeeze(Grad_ZX(2,:));
Grad2_xi   = squeeze(Grad2_ZX(2,:,:));

%% Compute UpsD

% Terms proportional to Gradient of xi 
Vec1(1,1)    = cos(eta)*sin(gamma)*cos(varphi)+sin(eta)*sin(varphi);
Vec1(2,1)    = cos(eta)*sin(gamma)*sin(varphi)-sin(eta)*cos(varphi);
Vec1(3,1)    = cos(eta)*cos(gamma);

Vec1_xs(1,1) =-cos(eta)*sin(gamma)*sin(varphi)+sin(eta)*cos(varphi);
Vec1_xs(1,2) = cos(eta)*cos(gamma)*cos(varphi);
Vec1_xs(1,3) =-sin(eta)*sin(gamma)*cos(varphi)+cos(eta)*sin(varphi);

Vec1_xs(2,1) = cos(eta)*sin(gamma)*cos(varphi)+sin(eta)*sin(varphi);
Vec1_xs(2,2) = cos(eta)*cos(gamma)*sin(varphi);
Vec1_xs(2,3) =-sin(eta)*sin(gamma)*sin(varphi)-cos(eta)*cos(varphi);

Vec1_xs(3,2) = -cos(eta)*sin(gamma);
Vec1_xs(3,3) = -sin(eta)*cos(gamma);


% Terms proportional to Gradient of zeta 
Vec2(1,1)    = sin(eta)*sin(gamma)*cos(varphi)-cos(eta)*sin(varphi);
Vec2(2,1)    = cos(eta)*cos(varphi)+sin(eta)*sin(gamma)*sin(varphi);
Vec2(3,1)    = sin(eta)*cos(gamma);

Vec2_xs(1,1) = -sin(eta)*sin(gamma)*sin(varphi)-cos(eta)*cos(varphi);
Vec2_xs(1,2) =  sin(eta)*cos(gamma)*cos(varphi);
Vec2_xs(1,3) =  cos(eta)*sin(gamma)*cos(varphi)+sin(eta)*sin(varphi);

Vec2_xs(2,1) =-cos(eta)*sin(varphi)+sin(eta)*sin(gamma)*cos(varphi);
Vec2_xs(2,2) = sin(eta)*cos(gamma)*sin(varphi);
Vec2_xs(2,3) =-sin(eta)*cos(varphi)+cos(eta)*sin(gamma)*sin(varphi);

Vec2_xs(3,2) = -sin(eta)*sin(gamma);
Vec2_xs(3,3) =  cos(eta)*cos(gamma);


UpsD      = -Vec1*Grad_xi - Vec2*Grad_zeta;

for k=1:1:4*PND.Kite.N   
   UpsD_xs(:,:,k)  = -Vec1*Grad2_xi(:,k)' - Vec2*Grad2_zeta(:,k)';
end
for k=1:1:3
  Ind = 4*(Kite-1)+k;
  UpsD_xs(:,:,Ind) = UpsD_xs(:,:,Ind) - Vec1_xs(:,k)*Grad_xi- Vec2_xs(:,k)*Grad_zeta;
end

end