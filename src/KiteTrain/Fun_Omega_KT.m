function [Phi Phi_xs] = Fun_Omega_KT(xs,PND)


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
% Inputs:  xs     -> state vector                                    %%
%          PND    -> Dimensionless parameters                        %%
% Outputs: Phi    -> matrix for the computation of the components    %%
%                  in SB of the dimensionless angular velocities of  %%
%                  of the kites with respect to SE:                  %%
%                  omega = Phi*xs_p                                  %%
%          Phi_xs -> tensor with the partial derivatives of Phi      %%
%                  with respect to the state vector components       %%
%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
Phi     = zeros(3,4*PND.Kite.N,PND.Kite.N);
Phi_xs  = zeros(3,4*PND.Kite.N,4*PND.Kite.N,PND.Kite.N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            Compute matrix phi                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:PND.Kite.N
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   % Recover the variables of the state Vector
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   varphi = xs(4*(i-1)+1,1);
   gamma  = xs(4*(i-1)+2,1);
   eta    = xs(4*(i-1)+3,1);
   theta  = xs(4*(i-1)+4,1);
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Compute Phi 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Terms with varphi_i  
   Phi(1,4*(i-1)+1,i) = -cos(gamma)*cos(eta)*sin(theta)-sin(gamma)*cos(theta);   
   Phi(2,4*(i-1)+1,i) =  cos(gamma)*sin(eta);                                    
   Phi(3,4*(i-1)+1,i) =  cos(gamma)*cos(eta)*cos(theta)-sin(gamma)*sin(theta);   

   % Terms with gamma_i
   Phi(1,4*(i-1)+2,i) =  sin(eta)*sin(theta);  
   Phi(2,4*(i-1)+2,i) =  cos(eta);            
   Phi(3,4*(i-1)+2,i) = -sin(eta)*cos(theta);  

   % Terms with eta_i
   Phi(1,4*(i-1)+3,i) = cos(theta);          
   Phi(3,4*(i-1)+3,i) = sin(theta); 
   
   % Terms with Theta_i
   Phi(2,4*(i-1)+4,i) = 1;
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Compute the gradient of Phi
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Derivative with respect to varphi is equal to zero
   %% Derivative with respect to gamma is:
    Phi_xs(1,4*(i-1)+1,4*(i-1)+2,i) =  sin(gamma)*cos(eta)*sin(theta)-cos(gamma)*cos(theta);   
    Phi_xs(2,4*(i-1)+1,4*(i-1)+2,i) = -sin(gamma)*sin(eta);                                    
    Phi_xs(3,4*(i-1)+1,4*(i-1)+2,i) = -sin(gamma)*cos(eta)*cos(theta)-cos(gamma)*sin(theta);   
   %% Derivative with respect to eta is:
    Phi_xs(1,4*(i-1)+1,4*(i-1)+3,i) =  cos(gamma)*sin(eta)*sin(theta);   
    Phi_xs(2,4*(i-1)+1,4*(i-1)+3,i) =  cos(gamma)*cos(eta);                       
    Phi_xs(3,4*(i-1)+1,4*(i-1)+3,i) = -cos(gamma)*sin(eta)*cos(theta);   
   
    Phi_xs(1,4*(i-1)+2,4*(i-1)+3,i) =  cos(eta)*sin(theta);  
    Phi_xs(2,4*(i-1)+2,4*(i-1)+3,i) = -sin(eta); 
    Phi_xs(3,4*(i-1)+2,4*(i-1)+3,i) = -cos(eta)*cos(theta);
    
   %% Derivative with respect to theta is:
    Phi_xs(1,4*(i-1)+1,4*(i-1)+4,i) = -cos(gamma)*cos(eta)*cos(theta)+sin(gamma)*sin(theta);             
    Phi_xs(3,4*(i-1)+1,4*(i-1)+4,i) = -cos(gamma)*cos(eta)*sin(theta)-sin(gamma)*cos(theta);   
   
    Phi_xs(1,4*(i-1)+2,4*(i-1)+4,i) =  sin(eta)*cos(theta);  
    Phi_xs(3,4*(i-1)+2,4*(i-1)+4,i) =  sin(eta)*sin(theta); 
    
    Phi_xs(1,4*(i-1)+3,4*(i-1)+4,i) = -sin(theta);
    Phi_xs(3,4*(i-1)+3,4*(i-1)+4,i) =  cos(theta);
end



end