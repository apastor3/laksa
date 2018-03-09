function RHS = Fun_ODE_Lag_KA(t,u)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : RHS of the equations (Lagrangian formulation)                  %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  t               -> dimensionless time                   %%
%          u = [xs xs_dot] -> extended state vector                %%
% Outputs: RHS             -> Right-Hand-Side of Lagrange Eqs.     %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define PND as global variables
global PND
% Recover the state vector
xs              = u(1:4,1);
xs_p            = u(5:8,1);

% Compute the control vector
[xc xc_p xc_pp] = Fun_Control_KA(t,PND);
 
% Rotation matrices
R1E    = Fun_R1E_KA(xs);  
R21    = Fun_R21_KA(xs);
RB2    = Fun_RB2_KA(xs);
RBE    = RB2*R21*R1E;
% Angular velocity in the body frame
[Phi Phi_xs ] = Fun_Matrix_Omega_KA(xs);
omega         = Phi*xs_p;
% Components of the dimensionless velocity of the center of mass in the body frame
[Ups_s Ups_c Ups_s_xs Ups_s_xc Ups_c_xs Ups_c_xc] = Fun_Matrix_Upsilon_KA(xs,xc,PND);
vg                                                = Ups_s*xs_p+Ups_c*xc_p;

% Kinematics Matrices
[Ms Ms_s Ms_c Msc Msc_s Msc_c ] = Fun_Matrix_M_KA(Phi,Phi_xs,Ups_s, Ups_c, Ups_s_xs, Ups_s_xc, Ups_c_xs, Ups_c_xc,PND);

% Componentns of the kite position vector in Earth frame 
OEO2            =  R1E'*[0,  xc(1,1)*sin(xs(3)+xc(2,1)), -xc(1,1)*cos(xs(3)+xc(2,1))]'; 
O2G             = -RBE'*[PND.Tether.XA 0 PND.Tether.ZA]';
rg              =  OEO2+O2G;
% Wind speed components in the Earth Frame
vw              = Fun_Wind(t,rg,PND);
% Wind velocity in the body frame
vw              = RBE*vw;

% Generalized forces

[Q f m alfa beta] = Fun_Q_KA(vw,vg,omega,Ups_s,Phi,PND);

% Gravitational force
[U U_xs] = Fun_Potential_KA(xs,xc,PND);

% Auxiliary operations
Aux1 = zeros(4,4); % (partial Ms/partial xs)*xs_p
Aux3 = zeros(4,2); % (partial Msc/partial xs)*xs_p
Aux5 = zeros(4,1); 
Aux6 = zeros(4,1);
for k=1:1:length(xs)
    Aux1(:,:) = Aux1 + Ms_s(:,:,k)*xs_p(k);
    Aux3(:,:) = Aux3 + Msc_s(:,:,k)*xs_p(k);
    Aux5(k,1) =   xs_p'*Ms_s(:,:,k)*xs_p;
    Aux6(k,1) = 2*xs_p'*Msc_s(:,:,k)*xc_p;
end
Aux2 = zeros(4,4); %(partial Ms/partial xc)*xs
Aux4 = zeros(4,2); %(partial Msc/partial xc)*xc
for k=1:1:length(xc)
    Aux2(:,:) = Aux2 +  Ms_c(:,:,k)*xc_p(k);
    Aux4(:,:) = Aux4 + Msc_c(:,:,k)*xc_p(k); 
end

% Construct the Righ-Hand Side of the equation
RHS0       = Q-U_xs + 0.5*(Aux5+Aux6)-(Aux1+Aux2)*xs_p-(Aux3+Aux4)*xc_p-Msc*xc_pp;
RHS(1:4,1) = xs_p; 
RHS(5:8,1) = Ms\RHS0;

end