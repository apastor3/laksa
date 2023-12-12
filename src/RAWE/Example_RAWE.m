function Example_RAWE

clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   This subroutine is the main programme to explore the               %%
%%%          Dynamics of the structure of a rotary AWE machine           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Common/')

format long

global PP NP PND s_fe s_nc
global B_alpha_inv Shape Shape_s Shape_ss Shape_sss
global phi phi_s phi_ss phi_sss varphi varphi_s varphi_ss
global P0_Inv Q0_Inv p0 P4 Q1 P_iN

% Load the Physical Parameters
Flag_Dim    = 1;                                % 1-> Plot results with units
Flag_Movie  = 1;
               
NP             = Load_Numerical_Parameters;
PND.Num.dh     = NP(12);
PND.Num.NewTol = NP(13);
PND.Num.MaxIt  = NP(14);
N_fe           = NP(1);

[PP PD]     = Load_Physical_Parameters;
% Initial Guess
X           = [-0.716095963570029   -0.705956187879401    0.066890506422876   -0.062630792228560   -0.705864774420679   -0.716014907261332   -0.067959164769445    0.062046448577242]';

% Parameter to control the integration
options  = odeset('RelTol',NP(10),'AbsTol',NP(11));
N_movie  = 10000;
T0       = 0;   % Dimensionless Time
TF       = 10;  % Dimensionless time

m_0_tt   = 0;
m_N_tt   = 0;

display('---------------------------------')
display(['mu    = ' num2str(PP(1))])
display(['beta  = ' num2str(PP(2))])
display(['sigma = ' num2str(PP(3))])
display(['delta = ' num2str(PP(4))])
display(['Gamma = ' num2str(PP(5))])
display(['epsilon = ' num2str(PP(6))])
display(['m_0   = ' num2str(PP(7))])
display(['m_N   = ' num2str(PP(8))])
display(['nu   = ' num2str(PP(12))])


% Load the Numerical Parameters 

display('---------------------------------')
display(['N_{fe} = ' num2str(NP(1))])
display(['N_Int  = ' num2str(NP(2))])
display(['Test   = ' num2str(NP(3))])

options  = odeset('RelTol',NP(10),'AbsTol',NP(11));

% Compute the Global function for r and alpha
[s_fe s_nc phi phi_s phi_ss  phi_sss varphi varphi_s varphi_ss Shape Shape_s Shape_ss Shape_sss B_alpha_inv] = Fun_Global(NP(1),NP(2));

% Compute matrices
[P0_Inv Q0_Inv p0 P4 Q1 P00 Q00 P_iN] = Compute_Matrix(NP(1),100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State vector        xs      = [u du/dt m_0];
%                     Dim(xs) = 2*Nu + 1
%                     u       = [x_1... x_{N_fe-1} y_1... y_{N_fe-1} z_1... z_{N_fe-1} alpha_1... alpha_{N_fe+1} ]
%                     Dim(u)  = Nu = 3*(N_fe -1 ) + N_fe + 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the initial guess

% Compute the Global function for r and alpha
[s_fe s_nc phi phi_s phi_ss  phi_sss varphi varphi_s varphi_ss Shape Shape_s Shape_ss Shape_sss B_alpha_inv] = Fun_Global(NP(1),NP(2));
% Compute vector V  for alpha, alpha = sum (alpha_j*varphi_j) + V_alpha
[V_alpha V_alpha_s V_alpha_ss V_alpha_tt] = Fun_Global_Valpha(NP,PP(7),PP(8),B_alpha_inv,Shape,Shape_s,Shape_ss,m_0_tt,m_N_tt);
% Compute matrices
[P0_Inv Q0_Inv p0 P4 Q1 P00 Q00 P_iN] = Compute_Matrix(NP(1),100);
%MM = Mass_RT(0);
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the stationary solution  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X0 Error EXITO] = my_fzero('Fun_Zero_Catenaria',X,PND);
X = X0;
x0               = [0 0  X(1,1)  X(2,1) 0 0  X(3,1)  X(4,1)]';     
[s xs]           = ode45('Fun_Der_Catenaria',s_fe,x0,options);
              

Nu   = 3*(N_fe-1) + N_fe + 1;
% Initialize the vector  
xs0  = zeros(2*Nu+1,1);
for i=1:1:N_fe-1
      xs0(0*(N_fe-1)+i,1) =  xs(i+1,1);
      xs0(1*(N_fe-1)+i,1) =  0;
      xs0(2*(N_fe-1)+i,1) =  xs(i+1,2);
end
X0 = [xs0(1) xs0(3)]';
[x Error EXITO]= my_fzero('Function_Equilibrium',X0,PND);
xs0(0*(N_fe-1)+1:1:1*(N_fe-1),1)    = x(0*(N_fe-1)+1:1*(N_fe-1),1);
xs0(2*(N_fe-1)+1:1:3*(N_fe-1),1)    = x(1*(N_fe-1)+1:2*(N_fe-1),1);

% Find the full vectors
[r alpha r_t alpha_t ] = From_xs_to_Physical(0,xs0,NP,PP);


 
% Make a call to the RHS for checking that ixs0 is a stationary solution    
[DF DG r_N alpha_N r_t alpha_t m_0 kappa tau Ten F_NL RHS_r RHS_alpha] = Fun_ODE_Full_RT(0,xs0);
display('---  Stationary Solution --- ')
display(['Error = ' num2str(max(abs(DF)))])
display(['Error = ' num2str(max(abs(DG)))])

% Compute the Jacobian
Jnum = Jacobian('Fun_ODE_RT_NC',0,xs0(1:end-1,1),PND);
% Compute Eigenvalues and Eigenvectors
[Vec Val ] = eig(Jnum);
display('---  Eigenvalues --- ')
Val        = diag(Val);
for i=1:1:length(Val)
   display(['Eigenvalues ' num2str(i) ' = ' num2str(Val(i))])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the initial condition
if Flag_Movie==1
    Fig = figure(1);
    Ax  = gca;
    Plot_Rotary(Fig,Ax,NP,r,alpha,phi,phi_s,phi_ss,phi_sss,varphi,varphi_s,varphi_ss,V_alpha,V_alpha_s,V_alpha_ss);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Perturba the initial condition
xs0      = xs0 + 1e-3*rand(length(xs0),1);

display('Integration is in progress')
[tau xs] = ode15s('Fun_ODE_RT',T0 + [0:1:N_movie]*(TF-T0)/N_movie,xs0,options);

display('Integration is Finishing')
% Plot the Results
Plot_Results(Fig,Ax,Flag_Dim,Flag_Movie,NP,PP,PD,tau,xs)


end