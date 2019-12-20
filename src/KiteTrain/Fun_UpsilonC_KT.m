function [UpsC UpsC_xs]= Fun_UpsilonC_KT(xs,xA,zA,ZX,PND)

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
% Inputs:  xs      -> state vector of a kite                        %%
%          xA,zA   -> SB coordinates of point A                     %% 
%          ZX      -> zeta and xi                                   %%
%          d_xs    -> x_s-Gradient of d                             %%
%                                                                   %%
% Outputs: UpsC    -> matrix for the computation of the components  %%
%                  in SE of the dimensionless velocity of the kites %%
%                  with respect to SE:                              %%
%                  v^n  = v^{n-1} + UpsC*xs_n_p + UpsD*xs_p         %%
%           UpsC_xs -> Derivative of UpsC with respect to xs,       %% 
%                      except the terms involving d                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the matrix
UpsC    = zeros(3,4);
UpsC_xs = zeros(3,4,4);

% Recover the variables
varphi = xs(1,1);
gamma  = xs(2,1);
eta    = xs(3,1);
theta  = xs(4,1);

% Recover auxiliary variables
zeta   = ZX(1,1);
xi     = ZX(2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute UpsC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UpsC(1,1) = (xA*cos(theta)+zA*sin(theta))*cos(gamma)*sin(varphi)+...
            (xi-xA*sin(theta)+zA*cos(theta))*(cos(eta)*sin(gamma)*sin(varphi)-sin(eta)*cos(varphi))...
            +zeta*(sin(eta)*sin(gamma)*sin(varphi) + cos(eta)*cos(varphi));

UpsC(1,2) = (xA*cos(theta)+zA*sin(theta))*sin(gamma)*cos(varphi)-...
            (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*cos(gamma)*cos(varphi)...
            -zeta*sin(eta)*cos(gamma)*cos(varphi);
        
UpsC(1,3) = (xi-xA*sin(theta)+zA*cos(theta))*(sin(eta)*sin(gamma)*cos(varphi)-cos(eta)*sin(varphi))...
            -zeta*(cos(eta)*sin(gamma)*cos(varphi)+sin(eta)*sin(varphi));

UpsC(1,4) = (xA*sin(theta)-zA*cos(theta))*cos(gamma)*cos(varphi)+...
            (xA*cos(theta)+zA*sin(theta))*(cos(eta)*sin(gamma)*cos(varphi) + sin(eta)*sin(varphi));
%%
UpsC(2,1) = -(xA*cos(theta)+zA*sin(theta))*cos(gamma)*cos(varphi)-...
             (xi-xA*sin(theta)+zA*cos(theta))*(cos(eta)*sin(gamma)*cos(varphi)+sin(eta)*sin(varphi))...
             +zeta*(cos(eta)*sin(varphi)-sin(eta)*sin(gamma)*cos(varphi));

UpsC(2,2) =  (xA*cos(theta)+zA*sin(theta))*sin(gamma)*sin(varphi)-...
             (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*cos(gamma)*sin(varphi)...
             -zeta*sin(eta)*cos(gamma)*sin(varphi);

UpsC(2,3) = (xi-xA*sin(theta)+zA*cos(theta))*(sin(eta)*sin(gamma)*sin(varphi)+cos(eta)*cos(varphi))...
             +zeta*(sin(eta)*cos(varphi)-cos(eta)*sin(gamma)*sin(varphi));

UpsC(2,4) = (xA*sin(theta)-zA*cos(theta))*cos(gamma)*sin(varphi)+...
            (xA*cos(theta)+zA*sin(theta))*(cos(eta)*sin(gamma)*sin(varphi) - sin(eta)*cos(varphi));
%%
UpsC(3,2) = (xA*cos(theta)+zA*sin(theta))*cos(gamma)+...
             (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*sin(gamma)+zeta*sin(eta)*sin(gamma);
        
UpsC(3,3) = (xi-xA*sin(theta)+zA*cos(theta))*sin(eta)*cos(gamma)-zeta*cos(eta)*cos(gamma);

UpsC(3,4) = -(xA*sin(theta)-zA*cos(theta))*sin(gamma)+...
            (xA*cos(theta)+zA*sin(theta))*cos(eta)*cos(gamma);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient Compute UpsC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UpsC(1,1)
UpsC_xs(1,1,1) =  (xA*cos(theta)+zA*sin(theta))*cos(gamma)*cos(varphi)+...
                  (xi-xA*sin(theta)+zA*cos(theta))*(cos(eta)*sin(gamma)*cos(varphi)+sin(eta)*sin(varphi))...
                  -zeta*(cos(eta)*sin(varphi)-sin(eta)*sin(gamma)*cos(varphi));
UpsC_xs(1,1,2) = -(xA*cos(theta)+zA*sin(theta))*sin(gamma)*sin(varphi)+...
                  (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*cos(gamma)*sin(varphi)...
                 +zeta*sin(eta)*cos(gamma)*sin(varphi);
UpsC_xs(1,1,3) =-(xi-xA*sin(theta)+zA*cos(theta))*(sin(eta)*sin(gamma)*sin(varphi)+cos(eta)*cos(varphi))...
                 -zeta*(sin(eta)*cos(varphi)-cos(eta)*sin(gamma)*sin(varphi));
UpsC_xs(1,1,4) =-(xA*sin(theta)-zA*cos(theta))*cos(gamma)*sin(varphi)-...
                 (xA*cos(theta)+zA*sin(theta))*(cos(eta)*sin(gamma)*sin(varphi)-sin(eta)*cos(varphi));
             
        
        
% UpsC(1,2)        
UpsC_xs(1,2,1) = -(xA*cos(theta)+zA*sin(theta))*sin(gamma)*sin(varphi)+...
                  (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*cos(gamma)*sin(varphi)...
                 +zeta*sin(eta)*cos(gamma)*sin(varphi);
UpsC_xs(1,2,2) = (xA*cos(theta)+zA*sin(theta))*cos(gamma)*cos(varphi)+...
                 (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*sin(gamma)*cos(varphi)...
                 +zeta*sin(eta)*sin(gamma)*cos(varphi);
UpsC_xs(1,2,3) = (xi-xA*sin(theta)+zA*cos(theta))*sin(eta)*cos(gamma)*cos(varphi)...
                 -zeta*cos(eta)*cos(gamma)*cos(varphi);
UpsC_xs(1,2,4) = -(xA*sin(theta)-zA*cos(theta))*sin(gamma)*cos(varphi)+...
                  (xA*cos(theta)+zA*sin(theta))*cos(eta)*cos(gamma)*cos(varphi);
     
        
% UpsC(1,3)
UpsC_xs(1,3,1) = -(xi-xA*sin(theta)+zA*cos(theta))*(sin(eta)*sin(gamma)*sin(varphi)+cos(eta)*cos(varphi))...
                  -zeta*(sin(eta)*cos(varphi)-cos(eta)*sin(gamma)*sin(varphi));;
UpsC_xs(1,3,2) =  (xi-xA*sin(theta)+zA*cos(theta))*(sin(eta)*cos(gamma)*cos(varphi))...
                  -zeta*cos(eta)*cos(gamma)*cos(varphi); ;
UpsC_xs(1,3,3) =  (xi-xA*sin(theta)+zA*cos(theta))*(cos(eta)*sin(gamma)*cos(varphi)+sin(eta)*sin(varphi))...
                  -zeta*(cos(eta)*sin(varphi)-sin(eta)*sin(gamma)*cos(varphi)); ;
UpsC_xs(1,3,4) = -(xA*cos(theta)+zA*sin(theta))*(sin(eta)*sin(gamma)*cos(varphi)-cos(eta)*sin(varphi));

        
% UpsC(1,4)
UpsC_xs(1,4,1) = -(xA*sin(theta)-zA*cos(theta))*cos(gamma)*sin(varphi)-...
                  (xA*cos(theta)+zA*sin(theta))*(cos(eta)*sin(gamma)*sin(varphi) - sin(eta)*cos(varphi));
UpsC_xs(1,4,2) = -(xA*sin(theta)-zA*cos(theta))*sin(gamma)*cos(varphi)+...
                  (xA*cos(theta)+zA*sin(theta))*cos(eta)*cos(gamma)*cos(varphi);
UpsC_xs(1,4,3) = -(xA*cos(theta)+zA*sin(theta))*(sin(eta)*sin(gamma)*cos(varphi) - cos(eta)*sin(varphi));
UpsC_xs(1,4,4) =  (xA*cos(theta)+zA*sin(theta))*cos(gamma)*cos(varphi)-...
                  (xA*sin(theta)-zA*cos(theta))*(cos(eta)*sin(gamma)*cos(varphi) + sin(eta)*sin(varphi));
       
%UpsC(2,1)
UpsC_xs(2,1,1) =  (xA*cos(theta)+zA*sin(theta))*cos(gamma)*sin(varphi)+...
                  (xi-xA*sin(theta)+zA*cos(theta))*(cos(eta)*sin(gamma)*sin(varphi)-sin(eta)*cos(varphi))...
                   +zeta*(cos(eta)*cos(varphi)+sin(eta)*sin(gamma)*sin(varphi));
UpsC_xs(2,1,2) =  (xA*cos(theta)+zA*sin(theta))*sin(gamma)*cos(varphi)-...
                  (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*cos(gamma)*cos(varphi)...
                  +zeta*(-sin(eta)*cos(gamma)*cos(varphi));
UpsC_xs(2,1,3) =  (xi-xA*sin(theta)+zA*cos(theta))*(sin(eta)*sin(gamma)*cos(varphi)-cos(eta)*sin(varphi))...
                  +zeta*(-sin(eta)*sin(varphi)-cos(eta)*sin(gamma)*cos(varphi));
UpsC_xs(2,1,4) =  (xA*sin(theta)-zA*cos(theta))*cos(gamma)*cos(varphi)+...
                  (xA*cos(theta)+zA*sin(theta))*(cos(eta)*sin(gamma)*cos(varphi)+sin(eta)*sin(varphi));

         
% UpsC(2,2)         
UpsC_xs(2,2,1) =  (xA*cos(theta)+zA*sin(theta))*sin(gamma)*cos(varphi)-...
                  (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*cos(gamma)*cos(varphi)...
                 -zeta*sin(eta)*cos(gamma)*cos(varphi);
UpsC_xs(2,2,2) =  (xA*cos(theta)+zA*sin(theta))*cos(gamma)*sin(varphi)+...
                  (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*sin(gamma)*sin(varphi)...
                 +zeta*sin(eta)*sin(gamma)*sin(varphi);
UpsC_xs(2,2,3) =  (xi-xA*sin(theta)+zA*cos(theta))*sin(eta)*cos(gamma)*sin(varphi)...
                  -zeta*cos(eta)*cos(gamma)*sin(varphi);
UpsC_xs(2,2,4) =  (-xA*sin(theta)+zA*cos(theta))*sin(gamma)*sin(varphi)+...
                  (xA*cos(theta)+zA*sin(theta))*cos(eta)*cos(gamma)*sin(varphi);
         
% UpsC(2,3)         
UpsC_xs(2,3,1) = (xi-xA*sin(theta)+zA*cos(theta))*(sin(eta)*sin(gamma)*cos(varphi)-cos(eta)*sin(varphi))...
                  -zeta*(sin(eta)*sin(varphi)+cos(eta)*sin(gamma)*cos(varphi));
UpsC_xs(2,3,2) = (xi-xA*sin(theta)+zA*cos(theta))*sin(eta)*cos(gamma)*sin(varphi)...
                  -zeta*cos(eta)*cos(gamma)*sin(varphi);
UpsC_xs(2,3,3) = (xi-xA*sin(theta)+zA*cos(theta))*(cos(eta)*sin(gamma)*sin(varphi)-sin(eta)*cos(varphi))...
                  +zeta*(cos(eta)*cos(varphi)+sin(eta)*sin(gamma)*sin(varphi));
UpsC_xs(2,3,4) =-(xA*cos(theta)+zA*sin(theta))*(sin(eta)*sin(gamma)*sin(varphi)+cos(eta)*cos(varphi));
   
% UpsC(2,4)
UpsC_xs(2,4,1) =  (xA*sin(theta)-zA*cos(theta))*cos(gamma)*cos(varphi)+...
                  (xA*cos(theta)+zA*sin(theta))*(cos(eta)*sin(gamma)*cos(varphi) + sin(eta)*sin(varphi));
UpsC_xs(2,4,2) = -(xA*sin(theta)-zA*cos(theta))*sin(gamma)*sin(varphi)+...
                  (xA*cos(theta)+zA*sin(theta))*cos(eta)*cos(gamma)*sin(varphi) ;
UpsC_xs(2,4,3) = -(xA*cos(theta)+zA*sin(theta))*(sin(eta)*sin(gamma)*sin(varphi) + cos(eta)*cos(varphi));
UpsC_xs(2,4,4) =  (xA*cos(theta)+zA*sin(theta))*cos(gamma)*sin(varphi)-...
                  (xA*sin(theta)-zA*cos(theta))*(cos(eta)*sin(gamma)*sin(varphi) - sin(eta)*cos(varphi));
    
% UpsC(3,2)
UpsC_xs(3,2,2) = -(xA*cos(theta)+zA*sin(theta))*sin(gamma)+...
                  (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*cos(gamma)+zeta*sin(eta)*cos(gamma);
UpsC_xs(3,2,3) = -(xi-xA*sin(theta)+zA*cos(theta))*sin(eta)*sin(gamma)+zeta*cos(eta)*sin(gamma);
UpsC_xs(3,2,4) = -(xA*sin(theta)-zA*cos(theta))*cos(gamma)-...
                  (xA*cos(theta)+zA*sin(theta))*cos(eta)*sin(gamma);
                
% UpsC(3,3)         
UpsC_xs(3,3,2) = -(xi-xA*sin(theta)+zA*cos(theta))*sin(eta)*sin(gamma)+zeta*cos(eta)*sin(gamma);
UpsC_xs(3,3,3) =  (xi-xA*sin(theta)+zA*cos(theta))*cos(eta)*cos(gamma)+zeta*sin(eta)*cos(gamma);
UpsC_xs(3,3,4) = -(xA*cos(theta)+zA*sin(theta))*sin(eta)*cos(gamma);


% UpsC(3,4)    
UpsC_xs(3,4,2) = -(xA*sin(theta)-zA*cos(theta))*cos(gamma)-...
                  (xA*cos(theta)+zA*sin(theta))*cos(eta)*sin(gamma);        
UpsC_xs(3,4,3) = -(xA*cos(theta)+zA*sin(theta))*sin(eta)*cos(gamma);        
UpsC_xs(3,4,4) = -(xA*cos(theta)+zA*sin(theta))*sin(gamma)-...
                  (xA*sin(theta)-zA*cos(theta))*cos(eta)*cos(gamma);        

end