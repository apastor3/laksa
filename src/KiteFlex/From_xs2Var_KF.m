function [xs xs_p] = From_xs2Var_KF(xs_amp,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : State vector in structure form                                 %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %%
% Input: State vector, first derivatives                                  %%
%        xs_amp = [gamma   varphi   theta   psi   phi ...                 %% 
%                  gamma_p varphi_p theta_p psi_p phi_p...                %%
%                  lambda_1_p... lambda_Ng_p]                             %%
% Output: State vector and its first derivate in structure form           %%

% Recover some dimensions
NR             = PND.Num.N;   % Number of Bars
NG             = PND.Gen.Num; % Number of Generators 
% State Vector
xs.gamma       = xs_amp(1:NR,:);
xs.varphi      = xs_amp(NR+1:2*NR,:);
xs.theta       = xs_amp(2*NR+1,:);
xs.psi         = xs_amp(2*NR+2,:);
xs.phi         = xs_amp(2*NR+3,:);

% Derivative of the state vectors
N0             = 2*NR+3;
xs_p.gamma     = xs_amp(N0+1:N0+NR,:);
xs_p.varphi    = xs_amp(N0+NR+1:N0+2*NR,:);
xs_p.theta     = xs_amp(N0+2*NR+1,:);
xs_p.psi       = xs_amp(N0+2*NR+2,:);
xs_p.phi       = xs_amp(N0+2*NR+3,:);
if NG>0
   xs_p.lambda = xs_amp(N0+2*NR+3+1:N0+2*NR+3+NG,:);
else
   xs_p.lambda = [];
end



end