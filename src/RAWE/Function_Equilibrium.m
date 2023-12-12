function F = Function_Equilibrium(t,x)

global NP

% Load useful Parameters
N_fe = NP(1);

X    = zeros(2*(4*N_fe-2),1);
% Compute the full state vector
X(0*(N_fe-1)+1:1:1*(N_fe-1),1)    = x(0*(N_fe-1)+1:1*(N_fe-1),1);
X(2*(N_fe-1)+1:1:3*(N_fe-1),1)    = x(1*(N_fe-1)+1:2*(N_fe-1),1);

DF = Fun_ODE_RT_NC(0,X);

N0 = 3*(N_fe-1) +  N_fe+ 1; 

F(0*(N_fe-1)+1:1*(N_fe-1),1) = DF(N0 + 0*(N_fe-1)+1: N0 + 1*(N_fe-1),1);
F(1*(N_fe-1)+1:2*(N_fe-1),1) = DF(N0 + 2*(N_fe-1)+1: N0 + 3*(N_fe-1),1);


end