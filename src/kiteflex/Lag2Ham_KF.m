function X_Ham = Lag2Ham_KF(t,X_Lag0,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : From Lagrangian to Hamiltonian state vector                    %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
% Input                                                                       %
% t      - > Time                                                             %
% X_Lag  - > Lagrangian variables X_Lag = [xs xs_p ]                          %
% PND    - > Dimensionless parameters                                         %
% Output                                                                      %
% X_Ham  - > Hamiltonian variables X_Ham = [xs p]                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nvar                         = 2*PND.Num.N+3;              % Number of Variables 
Nvar_p                       = 2*PND.Num.N+3+PND.Gen.Num;  % Number of Variables (dot) 
Nc                           = 4+PND.Gen.Num+3;            % Number of control parameters
Nc0                          = 4;                          % Number of control variables affecting the kinematics


if length(X_Lag0)>Nvar+Nvar_p
    X_Lag                = X_Lag0(1:Nvar+Nvar_p,1);
    X_Ham(Nvar+Nvar_p+1) = X_Lag0(end); 
else
     X_Lag               = X_Lag0;
end

xc_amp                       = Fun_Control_KF(t,X_Lag0,PND);

[xs_str xs_p_str]            = From_xs2Var_KF(X_Lag,PND);
[xc_str xc_p_str xc_pp_str]  = From_xc2Var_KF(xc_amp,PND);

[R_KE Grad_R_KE ]            = Matrix_R_KE_KF(xs_str);
[SR CR OmR SK CK OmK SG CG OmG] = Matrix_V_KF(xs_str,xc_str,R_KE,PND);

[Ms Msc Mc]                  = Matrix_M_KF(SK,CK,OmK,SR,CR,OmR,SG,CG,xs_str,xc_str,PND);
     
xs_p                         = X_Lag(Nvar+1:Nvar+Nvar_p,1);
xc_p                         = xc_amp(1*Nc+1:2*Nc,1);
p                            = Ms*xs_p+Msc*xc_p(1:Nc0,1); 
X_Ham                        = [X_Lag(1:Nvar); p]; 
  
 end