function xc  = Fun_Control_KT(t,PND)

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


xc = zeros(3,PND.Kite.N);
if PND.Ctr.Type == 0 
       for i=1:1:PND.Kite.N 
           xc(1,i)    = PND.Ctr.delta_a;                                           % delta_a 
           xc(2,i)    = PND.Ctr.delta_r;                                           % delta_r
           xc(3,i)    = PND.Ctr.delta_e;                                           % delta_e
       end
end

if PND.Ctr.Type == 1 
       for i=1:1:PND.Kite.N 
           xc(1,i)    = PND.Ctr.delta_a*cos(10*t);                                  % delta_a 
           xc(2,i)    = (-1)^i*PND.Ctr.delta_r*cos(4*t);                            % delta_r
           xc(3,i)    = PND.Ctr.delta_e;                                            % delta_e
       end
end



end