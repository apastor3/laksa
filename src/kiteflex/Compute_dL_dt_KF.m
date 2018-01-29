function Lt  = Compute_dL_dt_KF(Msc,Mc,Ms_xc,Msc_xc,Mc_xc,xs_amp,xc_amp,PND)
      
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Time derivative of the Lagrangian function                     %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
% Input                                                                      %
% Msc,Mc,Ms_xc,Msc_xc,Mc_xc - > Kinematics matrices and gradients            %
% xs_amp, xc_amp            - > State and control vectors                    %
% Output                                                                     %
% Lt          - > partial L/ partial t                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kite Parameters
iK        = PND.Kite.ik;
% Tether Parameter
sigmaR    = PND.Tether.Sigma;
iR        = PND.Tether.Ups;
% Rotor Parameters
if PND.Gen.Num>0
    iG        = PND.Gen.iota; % Generator tensor of inertia (divided by Sigma*lG^2)
    sigmaG    = PND.Gen.Sigma;
    lG        = PND.Gen.l;
else
    iG        = zeros(3,3);
    sigmaG    = 0;
    lG        = 0;
end
% Important Dimensions
Nvar                              = 2*PND.Num.N+3;              % Number of variables 
Nvar_p                            = 2*PND.Num.N+3+PND.Gen.Num;  % Number of variables (dot) 
Nc                                = 4+PND.Gen.Num+3;            % Number of controlparameters
NR                                = PND.Num.N;                  % Number of Bars
NG                                = PND.Gen.Num;
Nc0                               = 4;                          % Number of control variables affecting the kinematics


% Recover the derivatives of the state and control vectors
xs_p   = xs_amp(Nvar+1:Nvar+Nvar_p,1);
xc_p   = xc_amp(1*Nc+1:2*Nc,1);
xc_pp  = xc_amp(2*Nc+1:3*Nc,1);

% State and control vectors in structured form 
[xs_str xs_p_str]           = From_xs2Var_KF(xs_amp,PND);
[xc_str xc_p_str xc_pp_str] = From_xc2Var_KF(xc_amp,PND);

 
% Partial T/partial t
AUX1 = zeros(1,Nc0);
AUX2 = zeros(1,Nc0);
AUX3 = zeros(1,Nc0);

for i=1:1:Nc0
    AUX1(1,i) = xs_p'*squeeze(Ms_xc(:,:,i))*xs_p;
    AUX2(1,i) = xs_p'*squeeze(Msc_xc(:,:,i))*xc_p(1:Nc0,1);
    AUX3(1,i) = xc_p(1:Nc0,1)'*squeeze(Mc_xc(:,:,i))*xc_p(1:Nc0,1);
end

T_t = 1/2*(AUX1+2*AUX2+AUX3)*xc_p(1:4,1) + xs_p'*Msc*xc_pp(1:Nc0,1) + ... 
     (1/2)*xc_p(1:Nc0,1)'*Mc*xc_pp(1:Nc0,1) + (1/2)*xc_pp(1:Nc0,1)'*Mc*xc_p(1:Nc0,1);

% Partial U/Partial t
Ep_xc = zeros(1,Nc);
for j=1:1:NR
    Ep_xc(1,1) = Ep_xc(1,1)+  (1+NG*sigmaG)*sin(xs_str.gamma(j));
    for k=1:1:j
        if k==j
            g = 1/2;
        else
            g = 1;
        end
        Ep_xc(1,1) = Ep_xc(1,1) + 2*sigmaR*xc_str.lt*g*sin(xs_str.gamma(k));
    end
end
Ep_xc(1,2) = (1+NG*sigmaG)*(cos(xs_str.theta)*(cos(xc_str.delta)*sin(xc_str.eta)*sin(xs_str.phi)+sin(xc_str.delta)*cos(xs_str.phi))-sin(xs_str.theta)*cos(xc_str.delta)*cos(xc_str.eta));
Ep_xc(1,3) = (1+NG*sigmaG)*xc_str.lb*(cos(xs_str.theta)*(cos(xc_str.delta)*cos(xs_str.phi)-sin(xc_str.delta)*sin(xc_str.eta)*sin(xs_str.phi))+sin(xs_str.theta)*sin(xc_str.delta)*cos(xc_str.eta));
Ep_xc(1,4) = (1+NG*sigmaG)*xc_str.lb*(cos(xs_str.theta)*cos(xc_str.delta)*cos(xc_str.eta)*sin(xs_str.phi)+sin(xs_str.theta)*cos(xc_str.delta)*sin(xc_str.eta));

Ep_t       = Ep_xc*xc_p;

% Partial Time Derivative of the Lagrangian function 
Lt         = T_t-Ep_t;



end