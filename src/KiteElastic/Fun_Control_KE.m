function xc  = Fun_Control_KE(t,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Control inputs                                                 %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  t               -> dimensionless time                   %%
%                                                                  %%
% Outputs: xc, xc_p, xc,pp -> control vector and its derivatives   %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xc = zeros(3*PND.Kite.Num+2*PND.Tether.Num );

    for i=1:1:PND.Kite.Num
           xc(3*(i-1)+1,1)    = PND.Ctr.delta_a0(i) + PND.Ctr.delta_a1(i)*cos(PND.Ctr.Om_delta_a(i)*t); % delta_a 
           xc(3*(i-1)+2,1)    = PND.Ctr.delta_r0(i) + PND.Ctr.delta_r1(i)*cos(PND.Ctr.Om_delta_r(i)*t); % delta_r
           xc(3*(i-1)+3,1)    = PND.Ctr.delta_e0(i) + PND.Ctr.delta_e1(i)*cos(PND.Ctr.Om_delta_e(i)*t); % delta_e
    end
    for i=1:1:PND.Tether.Num
         xc(3*PND.Kite.Num+2*(i-1)+1,1) =                                                         PND.Tether.seda0(i)/(1+PND.Ctr.eps_L(i)*cos(PND.Ctr.Om_L(i)*t)); 
         xc(3*PND.Kite.Num+2*(i-1)+2,1) =-PND.Ctr.eps_L(i)*PND.Ctr.Om_L(i)*sin(PND.Ctr.Om_L(i)*t)*PND.Tether.seda0(i)/(1+PND.Ctr.eps_L(i)*cos(PND.Ctr.Om_L(i)*t))^2;
    end
    




end