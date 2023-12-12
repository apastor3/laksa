function [PP PD] = Load_Physical_Parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Introduce the Physical Parameter of the structure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
g         = 9.81;                % m/s^2
rho0      = 1.225;               % kg/m^3, ar density
% Physical Parameters
L         = 16.12;                  % m, length of the structure
R         = 0.19;                 % m, radius of the structure
M         = 0.51;                 % kg, mass of the structure

la        = 0.012;               % m, typical characteristic length
CD        = 1;                   % Drag coefficient
Vw        = 10;                  % m/s, Wind velocity

% Boundary Conditions
M0        = 0;                   % N*m Torque at s = 0;
M1        = 0;                  % N*m Torque at s = 1;
R1        = 0.;                 % m, circular motion of point s = 1
T1        = 0.0;                   % Period of the circularmotion
Omega1    = 0;%2*acos(-1)/T1 %40*sqrt(g/L);% 28.8;%2*acos(-1)           % rad/s, angular velocity of point 1


% Reference Values
T_ref     = 200;                 % N, Tension force measured in the experiments
Omega_ref = 0*acos(-1)/30;     % Rad/s, Angular velocity measure in the experiment

% Stiffness
EI        = 1.28e4;              % N*m^2, bending stiffness 
EA        = 3.64e4;              % N, axial stiffness 
GJ        = 140.0;               % N*m^2, torsional stiffness


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dimensionless parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu        = EI/(M*g*L^2);
beta      = GJ/(M*g*L^2);
sigma     = EA/(M*g);
delta     = beta*sigma/(2*mu);
nu        = rho0*la*L^2*CD/(2*M);
vw        = Vw/sqrt(g*L);
r1        = R1/L;
omega1    = Omega1*sqrt(L/g);

% Boundary conditions
Type      = 0;                   % 0-> M0 and M1 are fixed, 1-> M0 is given and M1 follows control law
Gamma0    = 45*acos(-1)/180;     % Angle of elevation
eps0      = T_ref/(EA);          % Distance between 0 and 1 is L(1+eps0);
m_0       = M0*L/(GJ);           % External Torque at s = 0
m_N       = M1*L/(GJ);           % External Torque at s = 1
tau0      = 0.0;                 % m_0 = m_N*(1-exp(-tau/tau0)) + m00*exp(-tau/tau_0)  
% Control Variables 
omega_tar = Omega_ref/sqrt(g/L); % Target dimensionless angular velocity
kp        = 0.;
kd        = 0.;


% Construct the vector
PP        = [mu beta sigma delta Gamma0 eps0 m_0 m_N omega_tar kp kd nu vw r1 omega1 tau0 Type];
PD        = [L  sqrt(L/g) M  Omega_ref GJ M1 EA EI R1 Omega1];




end