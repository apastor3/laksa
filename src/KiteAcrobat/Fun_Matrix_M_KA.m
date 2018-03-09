function [Ms Ms_s Ms_c Msc Msc_s Msc_c ] = Fun_Matrix_M_KA(Phi,Phi_xs,Ups_s, Ups_c, Ups_s_xs, Ups_s_xc, Ups_c_xs, Ups_c_xc,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Matrix M and its gradients                                     %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  Phi, Phi_xs        -> matrices (see Fun_Matrix_Omega)   %%
%          Ups_s, Ups_c,                                           %%
%          Ups_s_xs, Ups_s_xc                                      %%
%          Ups_c_xs, Ups_c_xc -> matrices (see Fun_Matrix_Upsilon) %%
%          PND                -> dimensionless parameters          %%
%                                                                  %%
% Outputs: Ms, Msc            -> matrices for the computation of   %%
%          the kinetic energy  (Mc is not need for Lagrange Eqs.)  %%
%          ek = 0.5*(xs_p'*Ms*xs_p+2xs_p'Msc*xc_p+xc_p'*Mc*xc_p)   %%
%                                                                  %%
%          Ms_s, Ms_c         -> tensors with the partial          %%
%          derivatives of Ms with respect to the components        %%
%          of xs and xc.                                           %%
%                                                                  %%
%          Msc_s, Msc_c       -> tensors with the partial          %%
%          derivatives of Msc with respect to the components       %%
%          of xs and xc.                                           %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Recover parameters
iG      = PND.Kite.iota ;

% Compute matrices Ms and Msc
Ms   = Ups_s'*Ups_s + Phi'*iG*Phi;
Msc  = Ups_s'*Ups_c;

% Compute tensors Ms_s, Ms_c, Msc_s, Msc_c
Ms_s = zeros(4,4,4);
Msc_s = zeros(4,2,4);
for i=1:1:4
   Ms_s(:,:,i)  = Ups_s_xs(:,:,i)'*Ups_s + Ups_s'*Ups_s_xs(:,:,i) + Phi_xs(:,:,i)'*iG*Phi + Phi'*iG*Phi_xs(:,:,i);
   Msc_s(:,:,i) = Ups_s_xs(:,:,i)'*Ups_c + Ups_s'*Ups_c_xs(:,:,i);
end
Ms_c = zeros(4,4,2);
Msc_c = zeros(4,2,2);
for i=1:1:2
   Ms_c(:,:,i)  = Ups_s_xc(:,:,i)'*Ups_s + Ups_s'*Ups_s_xc(:,:,i);
   Msc_c(:,:,i) = Ups_s_xc(:,:,i)'*Ups_c + Ups_s'*Ups_c_xc(:,:,i);
end

end