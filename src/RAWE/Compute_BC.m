function [r0 rN alpha_s_N r0_t rN_t rN_tt m_N_tt] = Compute_BC(tau,PP)



% Load Physical Parameters
Gamma0 = PP(5);
eps0   = PP(6);
r1     = PP(14);
omega1 = PP(15);

% Position: Boundary conditions at s = 0
r0        =  zeros(3,1);
r0_t      =  zeros(3,1);

% Position: Boundary conditions at s = 1
rN(1,1)   = -(1+eps0)*cos(Gamma0) - r1*sin(omega1*tau)*sin(Gamma0);
rN(2,1)   =  r1*cos(omega1*tau);
rN(3,1)   = -(1+eps0)*sin(Gamma0) + r1*sin(omega1*tau)*cos(Gamma0);


rN_t(1,1)   =  - omega1*r1*cos(omega1*tau)*sin(Gamma0);
rN_t(2,1)   =  - omega1*r1*sin(omega1*tau);
rN_t(3,1)   =  + omega1*r1*cos(omega1*tau)*cos(Gamma0);

rN_tt(1,1)   =  + omega1^2*r1*sin(omega1*tau)*sin(Gamma0);
rN_tt(2,1)   =  - omega1^2*r1*cos(omega1*tau);
rN_tt(3,1)   =  - omega1^2*r1*sin(omega1*tau)*cos(Gamma0);

% Angle: Boundary conditions at s = 1
if PP(17) == 0 % The torque is fixed
    m_N       =  PP(8); 
    m_N_t     =  0;
    m_N_tt    =  0;
else
    m_N       =  PP(8)*( 1 - exp(- tau/PP(16)) ); 
    m_N_t     = (PP(8)/PP(16))*exp(- tau/PP(16) );
    m_N_tt    = (PP(8)/PP(16)^2)*exp(- tau/PP(16) );
end

alpha_s_N = m_N;


end