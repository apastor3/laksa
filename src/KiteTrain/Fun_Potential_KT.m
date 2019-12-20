function [U U_xs] = Fun_Potential_KT(xs,ZX,Grad_ZX,PND)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : Potential Energy and its gradient                              %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  xs    -> state vector                                   %%
%          xc    -> control vector                                 %%
%          PND   -> dimensionless parameters of the system         %%
%                                                                  %%
% Outputs: U    -> Dimensionless potential                         %%
%          U_xs -> partial derivatives of U with respect to        %%
%                  the state vector components                     %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U         = 0;
U_xs      = zeros(4*PND.Kite.N,1);
ZE        = zeros(PND.Kite.N,1);
ZE_xs     = zeros(PND.Kite.N,4*PND.Kite.N);

Grad_zeta = squeeze(Grad_ZX(1,:,:));
Grad_xi   = squeeze(Grad_ZX(2,:,:));
   


for i=1:1:PND.Kite.N
   
   gamma   = xs(4*(i-1)+2,1);
   eta     = xs(4*(i-1)+3,1);
   theta   = xs(4*(i-1)+4,1);
   
    
   xA      = PND.Tether.XA(i);     % XA Attachment Point Coordinates
   zA      = PND.Tether.ZA(i);     % ZA Attachment Point Coordinates
   sigma   = PND.Kite.sigma(i);
   
   zeta      = ZX(1,i);
   xi        = ZX(2,i);
   
   ZE(i,1) = (xA*cos(theta)+zA*sin(theta))*sin(gamma)-(xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*cos(gamma)-zeta*sin(eta)*cos(gamma);
   if i>1
      ZE(i,1) = ZE(i-1,1) + ZE(i,1);
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%
   % Compute the Potential %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%
   U       = U - sigma*ZE(i,1);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Compute the gradient of the potential %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Ind          = 4*(i-1)+2; % Derivative with respect to gamma
   ZE_xs(i,Ind) = (xA*cos(theta)+zA*sin(theta))*cos(gamma)+(xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*sin(gamma)+zeta*sin(eta)*sin(gamma);
   
   Ind          = 4*(i-1)+3; % Derivative with respect to eta
   ZE_xs(i,Ind) = (xi-xA*sin(theta)+zA*cos(theta))*sin(eta)*cos(gamma)-zeta*cos(eta)*cos(gamma);
   
   Ind          = 4*(i-1)+4; % Derivative with respect to theta
   ZE_xs(i,Ind) = -(xA*sin(theta)-zA*cos(theta))*sin(gamma)+(xA*cos(theta)+zA*sin(theta))*cos(eta)*cos(gamma);
   
    
   ZE_xs(i,:)   =  ZE_xs(i,:) - cos(eta)*cos(gamma)*Grad_xi(:,i)' - sin(eta)*cos(gamma)*Grad_zeta(:,i)';
   
   if i>1
        ZE_xs(i,:)  =  ZE_xs(i-1,:) +  ZE_xs(i,:); 
   end
   
   U_xs = U_xs - sigma*ZE_xs(i,:)'; 

end



end