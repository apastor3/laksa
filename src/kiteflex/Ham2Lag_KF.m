function   X0_Lag = Ham2Lag_KF(t,X_Ham0,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : From Hamiltonian to Lagrangian state vector                    %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------

% Input
% t       - > Time
% X_Ham0  - > Hamiltonian variables X_Ham = [xs p]
% PND     - > Dimensionless variables
% Output
% X0_Lan  - > Lagrangian variables X0_Lan = [xs xs_p]

% Important dimensions
N                            = PND.Num.N; % Number of Bars
Nvar                         = 2*PND.Num.N+3;              % Number of variables 
Nvar_p                       = 2*PND.Num.N+3+PND.Gen.Num;  % Number of variables (dot) 
Nc                           = 4+PND.Gen.Num+3;              % Number of control parameters
Nc0                          = 4;           % Number of control variables affecting the kinematics
NR                           = PND.Num.N; % Number of Bars
NG                           = PND.Gen.Num;
% Check if dh/dt is computed
if length(X_Ham0)>Nvar+Nvar_p
    X_Ham                 = X_Ham0(1:Nvar+Nvar_p,1);
    X0_Lag(Nvar+Nvar_p+1) = X_Ham0(end); 
else
    X_Ham                 = X_Ham0;
end

xc_amp                       = Fun_Control_KF(t,X_Ham,PND);
[xc_str xc_p_str xc_pp_str]  = From_xc2Var_KF(xc_amp,PND);
xc_p                         = xc_amp(1*Nc+1:2*Nc,1);


xs_amp(1:Nvar,1)             = X_Ham(1:Nvar,1);
p                            = X_Ham(Nvar+1:Nvar+Nvar_p,1);
% Recover state vector in structured form
xs_str.gamma  = xs_amp(1:N,:);
xs_str.varphi = xs_amp(N+1:2*N,:); 
xs_str.theta  = xs_amp(2*N+1,:);
xs_str.psi    = xs_amp(2*N+2,:);
xs_str.phi    = xs_amp(2*N+3,:);

% Compute basic matrices
[R_KE Grad_R_KE ]                = Matrix_R_KE_KF(xs_str);
[SR CR OmR SK CK OmK SG CG OmG]  = Matrix_V_KF(xs_str,xc_str,R_KE,PND);
[Ms Msc Mc]                      = Matrix_M_KF(SK,CK,OmK,SR,CR,OmR,SG,CG,xs_str,xc_str,PND);


Hcp                       = Msc'/Ms;
xs_p                      = Ms\(p-Msc*xc_p(1:Nc0,1));

X0_Lag                   = [xs_amp;xs_p];

end