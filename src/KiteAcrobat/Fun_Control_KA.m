function [xc xc_p xc_pp] = Fun_Control_KA(t,PND)

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


switch PND.Ctr.Type
    case 0 % No control 
       xc(1,1)    = sqrt(1-PND.Tether.YA ^2);                                  % l
       xc(2,1)    = 0;                                                         % delta
    
       xc_p(1,1)  = 0;                                                         % l_p
       xc_p(2,1)  = 0;                                                         % delta_p
    
       xc_pp(1,1) = 0;                                                         % l_pp
       xc_pp(2,1) = 0;                                                         % delta_pp
       
    case 1 %l constant and delta modulated
       xc(1,1)     = sqrt(1-PND.Tether.YA^2);                                  % l
       xc(2,1)     = PND.Ctr.delta1*sin(PND.Ctr.Om_delta*t);                   % delta

       xc_p(1,1)   = 0;                                                        % l_p
       xc_p(2,1)   = PND.Ctr.delta1*PND.Ctr.Om_delta*cos(PND.Ctr.Om_delta*t);  % delta_p

       xc_pp(1,1)  = 0;                                                        % l_pp
       xc_pp(2,1)  = -PND.Ctr.Om_delta^2*xc(2,1);                              % delta_pp
        
    case 2  % l and delta are modulated      
        
        xc(1,1)     = sqrt(1-PND.Tether.YA^2) + PND.Ctr.l1*sin(PND.Ctr.Om_l*t); % l
        xc(2,1)     = PND.Ctr.delta1*sin(PND.Ctr.Om_delta*t);                   % delta

        xc_p(1,1)   = PND.Ctr.l1*PND.Ctr.Om_l*cos(PND.Ctr.Om_l*t);              % l_p
        xc_p(2,1)   = PND.Ctr.delta1*PND.Ctr.Om_delta*cos(PND.Ctr.Om_delta*t);  % delta_p

        xc_pp(1,1)  = -PND.Ctr.l1*PND.Ctr.Om_l^2*sin(PND.Ctr.Om_l*t);           % l_pp
        xc_pp(2,1)  = -PND.Ctr.Om_delta^2*xc(2,1);                              % delta_pp    
end

end