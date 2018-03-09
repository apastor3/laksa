function RHS = Fun_ODE_Ham_KA(t,w)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : RHS of the equations (Hamiltonian formulation)                 %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  t               -> dimensionless time                   %%
%          w = [xs p]      -> hamiltonian extended state vector    %%
% Outputs: RHS             -> Right-Hand-Side of Hamilton Eqs.     %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define PND as global variables
global PND
% Recover the state vector
xs              = w(1:4,1);
p               = w(5:8,1);

% Compute the control vector
[xc xc_p xc_pp] = Fun_Control_KA(t,PND);
 
% Rotation matrices
R1E    = Fun_R1E_KA(xs);  
R21    = Fun_R21_KA(xs);
RB2    = Fun_RB2_KA(xs);
RBE    = RB2*R21*R1E;
% Velocity matrices 
[Ups_s Ups_c Ups_s_xs Ups_s_xc Ups_c_xs Ups_c_xc] = Fun_Matrix_Upsilon_KA(xs,xc,PND);
[Phi Phi_xs ] = Fun_Matrix_Omega_KA(xs);
% Kinematics Matrices
[Ms Ms_s Ms_c Msc Msc_s Msc_c ] = Fun_Matrix_M_KA(Phi,Phi_xs,Ups_s, Ups_c, Ups_s_xs, Ups_s_xc, Ups_c_xs, Ups_c_xc,PND);

% Hamiltonian matrices
Hp      = inv(Ms);
Hcp     = Msc'*Hp;
%Hc     = Hcp*Msc-Mc; % This is not need for the equations of motion
xs_p    = Hp*(p-Msc*xc_p);
% Components of the dimensionless velocity and angular velocity of the kite in the body frame
omega   = Phi*xs_p;
vg      = Ups_s*xs_p+Ups_c*xc_p;
% Componentns of the kite position vector in Earth frame 
OEO2    =  R1E'*[0,  xc(1,1)*sin(xs(3)+xc(2,1)), -xc(1,1)*cos(xs(3)+xc(2,1))]'; 
O2G     = -RBE'*[PND.Tether.XA 0 PND.Tether.ZA]';
rg      =  OEO2+O2G;
% Wind speed components in the Earth Frame
vw      = Fun_Wind(t,rg,PND);
% Wind velocity in the body frame
vw      = RBE*vw;

% Generalized forces
[Q f m alfa beta] = Fun_Q_KA(vw,vg,omega,Ups_s,Phi,PND);

% Gravitational force
[U U_xs] = Fun_Potential_KA(xs,xc,PND);

% Compute tensors
Hp_s  = zeros(4,4,4);
Hcp_s = zeros(2,4,4);
Hc_s  = zeros(2,2,4);
for i=1:1:4
    Hp_s(:,:,i)  = -Hp*Ms_s(:,:,i)*Hp;
    Hcp_s(:,:,i) =  squeeze(Msc_s(:,:,i))'*Hp+Msc'*Hp_s(:,:,i);
    Hc_s(:,:,i)  =  Hcp_s(:,:,i)*Msc+Hcp*Msc_s(:,:,i);
end

% Compute the RHS
RHS0 = zeros(4,1);
for i=1:1:4
   RHS0(i,1) = p'*Hp_s(:,:,i)*p-2*xc_p'*Hcp_s(:,:,i)*p+xc_p'*Hc_s(:,:,i)*xc_p; 
end
% Ordinary diferential Equations
RHS(1:4,1) = Hp*p-Hcp'*xc_p;
RHS(5:8,1) = -0.5*RHS0-U_xs+Q;


end