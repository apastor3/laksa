function F = Fun_Zero_Catenaria(s,XX)

global PP NP

options = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Compute boundary conditions
[r0 rN alpha_s_N r0_t rN_t] = Compute_BC(0,PP);

% Prepare initial guess
x0     = [0 0 XX(1,1) XX(2,1) 0 0 XX(3,1) XX(4,1)]';
% Make a shot from s = 0
[s xs0] = ode45('Fun_Der_Catenaria',[0 0.5],x0,options );

% Prepare initial guess
x0     = [rN(1,1) rN(3,1) XX(5,1) XX(6,1) 0 0 XX(7,1) XX(8,1)]';
% Make a shot from s = 1
[s xs1] = ode45('Fun_Der_Catenaria',[1 0.5],x0,options);


F       = xs0(end,:)'-xs1(end,:)';


end