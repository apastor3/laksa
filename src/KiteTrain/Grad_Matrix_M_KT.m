function M_xs = Grad_Matrix_M_KT(Ups,Ups_xs,Phi,Phi_xs,PND)
                                                                
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : Gradient of Matrix M                                           %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:                                                                   %%
%        Ups, Ups_xs        -> Matrices Ups and their gradients             %%
%        Phi, Phi_xs        -> Matrices Phi and their gradients             %%
%        PND                -> Dimensionless parameter                      %%
% Outputs                                                                   %%
%         M_xs              -> Gradient of M with respect to xs             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gradients with respect to xs  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_xs   = zeros(4*PND.Kite.N,4*PND.Kite.N,4*PND.Kite.N);

for j=1:1:PND.Kite.N % Kite Loops
     sigma   = PND.Kite.sigma(j);
     for i=1:1:4*PND.Kite.N % Variable Loop 
           M_xs(:,:,i) = M_xs(:,:,i)+...
                         sigma*squeeze(Ups_xs(:,:,i,j))'*squeeze(Ups(:,:,j)) + ...
                         sigma*squeeze(Ups(:,:,j))'*squeeze(Ups_xs(:,:,i,j)) + ... 
                         squeeze(Phi_xs(:,:,i,j))'*squeeze(PND.Kite.iota(:,:,j))*squeeze(Phi(:,:,j)) + ...
                         squeeze(Phi(:,:,j))'*squeeze(PND.Kite.iota(:,:,j))*squeeze(Phi_xs(:,:,i,j));  
     end
end



end