function [Ms Ms_s ] = Fun_Matrix_M_KS(Phi,Phi_xs,Ups_s, Ups_s_xs,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Kinetic energy matrix and its gradient                         %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  Phi, Phi_xs        -> matrices (see Fun_Matrix_Omega)   %%
%          Ups_s,Ups_s_xs,    -> matrices (see Fun_Matrix_Upsilon) %%
%          PND                -> dimensionless parameters          %%
%                                                                  %%
% Outputs: Ms                 -> matrix  for the computation of    %%
%          the kinetic energy                                      %%
%          ek = 0.5*(xs_p'*Ms*xs_p)                                %%
%                                                                  %%
%          Ms_s               -> tensors with the partial          %%
%          derivatives of Ms with respect to the components        %%
%          of xs.                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Recover parameters
iG   = PND.Kite.iota ;
% Compute matrices Ms 
Ms   = Ups_s'*Ups_s + Phi'*iG*Phi;

% Compute tensors Ms_s
Ms_s = zeros(5,5,5);
for i=1:1:5
   Ms_s(:,:,i)  = Ups_s_xs(:,:,i)'*Ups_s + Ups_s'*Ups_s_xs(:,:,i) + Phi_xs(:,:,i)'*iG*Phi + Phi'*iG*Phi_xs(:,:,i);
end


end