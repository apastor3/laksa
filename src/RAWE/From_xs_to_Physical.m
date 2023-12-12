function [r alpha r_t alpha_t ] = From_xs_to_Physical(tau,xs0,NP,PP)

% Load Numerical Parameter
N_fe   = NP(1);
Nu     = 3*(N_fe - 1) + N_fe +1;
% Compute Boundary conditions
[r0 rN alpha_s_N r0_t rN_t rN_tt m_N_tt] = Compute_BC(tau,PP);

% Initialize the vectors
r        = zeros(N_fe+1,3);
alpha    = zeros(N_fe+1,1);
r_t      = zeros(N_fe+1,3);
alpha_t  = zeros(N_fe+1,1);

% Evaluate r
r(1,:)      = r0';
for j=1:1:3
    r(2:N_fe,j) = xs0((j-1)*(N_fe-1)+1:j*(N_fe-1),1); 
end
r(N_fe+1,:) = rN';
% Evaluate alpha
alpha(:,1)  = xs0(3*(N_fe-1)+1:3*(N_fe-1)+ N_fe + 1 ,1);

% Evaluate r_t
r_t(1,:)      = r0_t';
for j=1:1:3
    r_t(2:N_fe,j) = xs0(Nu + (j-1)*(N_fe-1)+1:Nu + j*(N_fe-1),1); 
end
r_t(N_fe+1,:) = rN_t';

% Evaluate alpha_t
alpha_t(:,1)  = xs0(Nu + 3*(N_fe-1)+1:Nu + 3*(N_fe-1)+ N_fe + 1 ,1);





end