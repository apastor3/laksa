function   u = From_Ham2Lag_KA(t,w,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : From Hamiltonian to Lagrangian state vector                    %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  t               -> dimensionless time                   %%
%          w = [xs p]      -> hamiltonian extended state vector    %%       
%          PND             -> dimensionless parameters             %%
%                                                                  %%
% Outputs: u = [xs xs_dot] -> lagrangian extended state vector     %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recover the state vector
xs              = w(1:4,1);
p               = w(5:8,1);

% Compute the control vector
[xc xc_p xc_pp] = Fun_Control_KA(t,PND);
% Velocity Matrices
[Phi Phi_xs ]   = Fun_Matrix_Omega_KA(xs);
[Ups_s Ups_c Ups_s_xs Ups_s_xc Ups_c_xs Ups_c_xc] = Fun_Matrix_Upsilon_KA(xs,xc,PND);
% Kinematics Matrices
[Ms Ms_s Ms_c Msc Msc_s Msc_c ] = Fun_Matrix_M_KA(Phi,Phi_xs,Ups_s, Ups_c, Ups_s_xs, Ups_s_xc, Ups_c_xs, Ups_c_xc,PND);

% Hamiltonian matrices
Hp      = inv(Ms);
Hcp     = Msc'*Hp;
xs_p    = Hp*(p-Msc*xc_p);
% Lagrangian extended state vector
u       = [xs;xs_p];

end