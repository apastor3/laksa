function [Tension TMod Elong Elong_t] = Tension_KS(xk,xk_p,xc,xc_p,sign,R2E,BarD_Kin_2,PND)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Author    : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Tension due to the springs                                     %
% Copyright:  Universidad Carlos II de Madrid, 2017. All rights reserved     %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %%
% Inputs: xk         -> kite state vector                            %%
%         xk_p       -> time derivative of xk                        %%
%         xc         -> control vector                               %%
%         xc_p       -> time derivative of the control vector        %%
%         sign       -> +1 compute the tension at C+B+ tether        %%  
%                       -1 compute the tension at C-B- tether        %%
%         R2E        -> S2-SE rotation matrix                        %% 
%         BarD_Kin_2 -> [d d_p ]                                     %%
%         PND        -> dimensionless parameters                     %%
%                                                                    %%
% Outputs: Tension -> Tension acting on the kite (SE-components)     %%
%          TMod    -> Tension Modulus                                %%
%          Elong   -> Tether Elongation                              %%
%          Elong_t -> Time derivative of tether Elongation           %%
%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recover the state vector an its first time derivative
varphi   = xk(1,1);
gamma    = xk(2,1);
eta      = xk(3,1);
theta    = xk(4,1); 
chi      = xk(5,1);

varphi_p = xk_p(1,1);
gamma_p  = xk_p(2,1);
eta_p    = xk_p(3,1);
theta_p  = xk_p(4,1); 
chi_p    = xk_p(5,1);
% Recover the control vector an its first time derivative
PR       = xc(1,1);
nu       = xc(2,1);
PR_p     = xc_p(1,1);
nu_p     = xc_p(2,1);

% Recover kinematic quantities related with the distance between OE and C0
d        = BarD_Kin_2(1,1);
d_p      = BarD_Kin_2(2,1);

% Recover Parameter
xA       = PND.Tether.XA;    yA       = PND.Tether.YA;    zA       = PND.Tether.ZA;
xB       = PND.Tether.XB;    yB       = PND.Tether.YB;    zB       = PND.Tether.ZB;

eps_cb   = PND.Bar.Lc;

l0       = PND.Tether.l;
ft       = PND.Tether.ft;
seda     = PND.Tether.seda;
kappa    = PND.Tether.kappa;

% Define auxiliary variables
xBA      = xB-xA;
zBA      = zB-zA;
lambda   = nu+chi;
lambda_p = nu_p+chi_p;

% Vector CB (S2 components )
CB(1,1) = xBA*cos(theta)+zBA*sin(theta);
CB(2,1) = sign*yB-sign*eps_cb/2*cos(lambda)+d*sin(chi);
CB(3,1) = -(d*cos(chi)+l0+sign*eps_cb/2*sin(lambda)+xBA*sin(theta)-zBA*cos(theta));
% Vector CB (SE components)
CB      = R2E'*CB;

% unit vector along BC, and dimensionless elongation
L            = sqrt(CB'*CB);
ut           = CB/L;
Elong(1,1)   = seda*L - 1;   
% Compute L^2 -> only for checking purposes  
%Mod2 = xBA^2+zBA^2+yB^2+l0^2+(eps_cb/2)^2+d^2+sign*l0*eps_cb*sin(lambda)...
%        -yB*(eps_cb*cos(lambda)-2*sign*d*sin(chi))+...
%        2*d*(l0*cos(chi)+sign*eps_cb/2*sin(nu))+...
%        2*(d*cos(chi)+l0+sign*eps_cb/2*sin(lambda))*(xBA*sin(theta)-zBA*cos(theta));
%Error = L-sqrt(Mod2)
   
% Compute time derivative of the elongation  
Aux1    =  d*d_p...
          +sign*eps_cb/2*l0*lambda_p*cos(lambda)...
          +yB*(eps_cb/2*lambda_p*sin(lambda)+sign*d*chi_p*cos(chi)+sign*d_p*sin(chi) )...
          +d_p*(l0*cos(chi)+sign*eps_cb/2*sin(nu))+...
          +d*(sign*eps_cb/2*nu_p*cos(nu)-l0*chi_p*sin(chi))...    
          +(d_p*cos(chi)-d*chi_p*sin(chi)+sign*eps_cb/2*lambda_p*cos(lambda))*(xBA*sin(theta)-zBA*cos(theta))...
          +theta_p*(d*cos(chi)+l0+sign*eps_cb/2*sin(lambda))*(xBA*cos(theta)+zBA*sin(theta)); 
   
Elong_t(1,1) = seda*Aux1/L;      
    
% Compute the tensions at the spring

if Elong(1,1)>0
   TMod          =   kappa*(Elong + ft*Elong_t);
   Tension(:,1)  = - TMod*ut;
else
   TMod     = zeros(1,1);
   Tension  = zeros(3,1);
end


end
