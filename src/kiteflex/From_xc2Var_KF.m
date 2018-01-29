function [xc xc_p xc_pp] = From_xc2Var_KF(xc_amp,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Control vector in structure form                               %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%----------------------------------------------------------------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %%
% Input:                                                                  %%
% xc_amp = [lt   lb   delta   eta    xi   delta_a   delta_r  delta_e .... %%
%       lt_p lb_p delta_p eta_p  xi_p delta_ap  delta_rp delta_ep...]     %%
%       lt_pp lb_pp delta_pp eta_pp xi_pp delta_app  delta_rpp delta_epp] %%
%                                                                         %%
% Output: control vector, first and second derivativesin structure form   %%                                                                                                                                         %%
%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NG            = PND.Gen.Num; % Number of Generators 

% Control variables
xc.lt         = xc_amp(1,:);
xc.lb         = xc_amp(2,:);
xc.delta      = xc_amp(3,:);
xc.eta        = xc_amp(4,:);
xc.xi         = xc_amp(5:4+NG,:);
xc.delta_a    = xc_amp(4+NG+1,:);
xc.delta_r    = xc_amp(4+NG+2,:);
xc.delta_e    = xc_amp(4+NG+3,:);

% First derivatives
xc_p.lt       = xc_amp(8+NG,:);
xc_p.lb       = xc_amp(9+NG,:);
xc_p.delta    = xc_amp(10+NG,:);
xc_p.eta      = xc_amp(11+NG,:);
xc_p.xi       = xc_amp(12+NG:11+2*NG,:);
xc_p.delta_a  = xc_amp(12+2*NG,:);
xc_p.delta_r  = xc_amp(13+2*NG,:);
xc_p.delta_e  = xc_amp(14+2*NG,:);

% Second derivatives
xc_pp.lt      = xc_amp(15 +2*NG,:);
xc_pp.lb      = xc_amp(16+2*NG,:);
xc_pp.delta   = xc_amp(17+2*NG,:);
xc_pp.eta     = xc_amp(18+2*NG,:);
xc_pp.xi      = xc_amp(19+2*NG:18+3*NG,:);
xc_pp.delta_a = xc_amp(19+3*NG,:);
xc_pp.delta_r = xc_amp(20+3*NG,:);
xc_pp.delta_e = xc_amp(21+3*NG,:);


end