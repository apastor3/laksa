function [xc xc_p xc_pp] = Fun_Control_KS(t,X,PND,rk,vk,ak)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Define Control variables                                       %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %%
% Inputs:  t      -> dimensionless time                                   %%
%          PND    -> dimensionless parameters                             %%  
%                                                                         %%
% Outputs: xc  = [PR lambda] -> Control variables                         %%
%          xc_p   ->  Time derivative of the control vector               %%
%          xc_pp  -> Second time derivative of the control vector         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initialize the vector
xc      = zeros(2,1);
xc_p    = zeros(2,1);
xc_pp   = zeros(2,1);

xc(1,1) = PND.Ctr.PR0;


if PND.Ctr.Type==1  %  Pull Up Maneuver
    PR0          = PND.Ctr.PR0;
    PRf          = PND.Ctr.PRf;
    Tr           = PND.Ctr.Tr;
    Ts           = PND.Ctr.Ts;
    Xi           = (t-Ts)/Tr;
    
    xc(1,1)    = PR0+(PRf-PR0)*(1+tanh(Xi))/2;                 % PR
    xc_p(1,1)  = (PRf-PR0)/(2*Tr)*(1-(tanh(Xi))^2);            % PR_p
    xc_pp(1,1) = -(PRf-PR0)/(Tr^2)*tanh(Xi)*(1-(tanh(Xi))^2);  % PR_pp
  
end


if PND.Ctr.Type==2  %  Sinusoidal lateral variation   
    xc(2,1)    =  PND.Ctr.Lam*sin(PND.Ctr.OmLam*t);                   % nu
    xc_p(2,1)  =  PND.Ctr.OmLam*PND.Ctr.Lam*cos(PND.Ctr.OmLam*t);     % nu_p
    xc_pp(2,1) = -PND.Ctr.OmLam^2*PND.Ctr.Lam*sin(PND.Ctr.OmLam*t);   % nu_pp 
end


if PND.Ctr.Type==3 % Close-Loop Make Figure of Eight
    
    % Target Trajectory
    y0         = PND.Ctr.Y0;
    omega_star = PND.Ctr.Om;
    
    y_star     =  y0*sin(omega_star*t);
    y_star_p   =  y0*omega_star*cos(omega_star*t);
    y_star_pp  = -y0*omega_star^2*sin(omega_star*t);
    
    E_star     = rk(2)-y_star;
    E_star_p   = vk(2)-y_star_p;
  
    
    xc(2,1)    = X(11);                                 % nu
    if isempty(ak)
       xc_p(2,1)  = PND.Ctr.GI*E_star+PND.Ctr.GP*E_star_p; % nu_p
    else  
       E_star_pp  = ak(2) - y_star_pp;
       xc_p(2,1)  = PND.Ctr.GI*E_star+PND.Ctr.GP*E_star_p+PND.Ctr.GD*E_star_pp; %nu_p
    end
end


end