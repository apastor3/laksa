function [Ups Ups_xs]  = Fun_Upsilon_KT(xs,ZX,Grad_ZX,Grad2_ZX,PND)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : SB-S2 rotation matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %%
% Inputs:  xs       -> state vector                                 %%
%          ZX       -> Zeta and xi                                  %%
%          Grad_ZX  -> Gradient of ZX                               %% 
%          Grad2_ZX -> Gradient of grad_ZX                          %%
%          PND      -> Dimensionless parameters                     %%
%                                                                   %%
% Outputs: Ups      -> matrix for the computation of the components %%
%                  in SE of the dimensionless velocity of the kites %%
%                  with respect to SE:                              %%
%                  vk  = Ups*xs_p                                   %%
%          Ups_xs   -> tensor with the partial derivatives of Ups   %%
%                  with respect to the state vector components      %%
%                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
Ups       = zeros(3,4*PND.Kite.N,PND.Kite.N);
Ups_xs    = zeros(3,4*PND.Kite.N,4*PND.Kite.N,PND.Kite.N);
UpsC_xi   = zeros(3,4);
UpsC_zeta = zeros(3,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            Compute matrix Ups                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First kite
[UpsC UpsC_xs]      = Fun_UpsilonC_KT(xs(1:4,1),PND.Tether.XA(1),PND.Tether.ZA(1),ZX(:,1),PND);
Ups(:,1:4,1)        = UpsC;
Ups_xs(:,1:4,1:4,1) = UpsC_xs; 

Grad_zeta  = squeeze(Grad_ZX(1,:,:));
Grad_xi    = squeeze(Grad_ZX(2,:,:));

for i=2:1:PND.Kite.N
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Recover the variables %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    varphi              = xs(4*(i-1)+1,1);
    gamma               = xs(4*(i-1)+2,1);
    eta                 = xs(4*(i-1)+3,1);
   
    %%%%%%%%%%%%%%%%%%%%
    % Previous Kite   %% 
    %%%%%%%%%%%%%%%%%%%%
    Ups(:,:,i)          = Ups(:,:,i-1);      % Previous kite
    Ups_xs(:,:,:,i)     = Ups_xs(:,:,:,i-1); % Previous kite  
    
    %%%%%%%%%%%%%%%%%%%%%%
    % UpsC Contribution %%
    %%%%%%%%%%%%%%%%%%%%%%
    Ind0                = 4*(i-1)+1;
    IndF                = 4*i;
    
    % UpsC
    [UpsC UpsC_xs]      = Fun_UpsilonC_KT(xs(Ind0:IndF),PND.Tether.XA(i),PND.Tether.ZA(i),ZX(:,i),PND);
    Ups(:,Ind0:IndF,i)  = Ups(:,Ind0:IndF,i) + UpsC;    
    
    % Gradient of UpsC
    Ups_xs(:,Ind0:IndF,Ind0:IndF,i)= Ups_xs(:,Ind0:IndF,Ind0:IndF,i) + UpsC_xs;
    
    UpsC_xi(1,1) =  cos(eta)*sin(gamma)*sin(varphi)-sin(eta)*cos(varphi);
    UpsC_xi(1,2) = -cos(eta)*cos(gamma)*cos(varphi);
    UpsC_xi(1,3) =  sin(eta)*sin(gamma)*cos(varphi)-cos(eta)*sin(varphi);
    
    UpsC_xi(2,1) = -(cos(eta)*sin(gamma)*cos(varphi)+sin(eta)*sin(varphi));
    UpsC_xi(2,2) = -cos(eta)*cos(gamma)*sin(varphi);
    UpsC_xi(2,3) =  (sin(eta)*sin(gamma)*sin(varphi)+cos(eta)*cos(varphi));
   
    UpsC_xi(3,2) = cos(eta)*sin(gamma);
    UpsC_xi(3,3) = sin(eta)*cos(gamma);
    
    %
    UpsC_zeta(1,1) =  cos(eta)*cos(varphi)+sin(eta)*sin(gamma)*sin(varphi);
    UpsC_zeta(1,2) =  -sin(eta)*cos(gamma)*cos(varphi);
    UpsC_zeta(1,3) =  -sin(eta)*sin(varphi)-cos(eta)*sin(gamma)*cos(varphi);
    
    UpsC_zeta(2,1) =  cos(eta)*sin(varphi)-sin(eta)*sin(gamma)*cos(varphi);
    UpsC_zeta(2,2) = -sin(eta)*cos(gamma)*sin(varphi);
    UpsC_zeta(2,3) =  sin(eta)*cos(varphi)-cos(eta)*sin(gamma)*sin(varphi);
   
    UpsC_zeta(3,2) = sin(eta)*sin(gamma);
    UpsC_zeta(3,3) =-cos(eta)*cos(gamma);
    
    for k=1:1:4*PND.Kite.N
        Ups_xs(:,Ind0:IndF,k,i)= Ups_xs(:,Ind0:IndF,k,i) + UpsC_xi*Grad_xi(k,i) + UpsC_zeta*Grad_zeta(k,i);
    end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % UpsD Contribution   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    [UpsD UpsD_xs]       = Fun_UpsilonD_KT(xs(Ind0:IndF),Grad_ZX(:,:,i),Grad2_ZX(:,:,:,i),i,PND);
    
    Ups(:,:,i)           = Ups(:,:,i) +  UpsD; 
    Ups_xs(:,:,:,i)      = Ups_xs(:,:,:,i) + UpsD_xs;
    
end


    


end