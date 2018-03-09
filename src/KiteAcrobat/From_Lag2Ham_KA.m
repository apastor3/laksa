function w = From_Lag2Ham_KA(t,u,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : From Lagrangian to Hamiltonian state vector                    %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  t               -> dimensionless time                   %%
%          u = [xs xs_dot] -> lagrangian extended state vector     %%
%          PND             -> dimensionless parameters             %%
%                                                                  %%
% Outputs: w = [xs p]      -> hamiltonian extended state vector    %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recover the state vector
xs              = u(1:4,1);
xs_p            = u(5:8,1);
% Compute the control vector
[xc xc_p xc_pp] = Fun_Control_KA(t,PND);
% Velocity Matrices
[Phi Phi_xs ]   = Fun_Matrix_Omega_KA(xs);
[Ups_s Ups_c Ups_s_xs Ups_s_xc Ups_c_xs Ups_c_xc] = Fun_Matrix_Upsilon_KA(xs,xc,PND);
% Kinematics Matrices
[Ms Ms_s Ms_c Msc Msc_s Msc_c ] = Fun_Matrix_M_KA(Phi,Phi_xs,Ups_s, Ups_c, Ups_s_xs, Ups_s_xc, Ups_c_xs, Ups_c_xc,PND);
% Compute the extended state vector of Hamilton equations
p               = Ms*xs_p+Msc*xc_p; 
w               = [xs;p]; 

  
 end