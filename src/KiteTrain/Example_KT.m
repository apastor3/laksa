

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : Example                                                        %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %%
% Inputs: No inputs                                             %%
%                                                               %%
% Outputs: Integration Results are placed in the workspace      %% 
%                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Define the dimensionless parameters as global variables
global PND

% Add the path of the common folder
addpath('../Common/')

% Load the physical parameters of KiteAcrobat
PD  = Fun_PD_KT;

% Find the dimensionless parameters of KiteAcrobat
PND = Fun_PND_KT(PD);

% Output Dimensions and time
Flag_Dim = 1;    % 0 - > Dimensionless Outputs, 1 -> Outputs with dimensions 
TF       = 10.;  % Final dimensionless integration time
options  = odeset('RelTol',PND.Num.RelTol,'AbsTol',PND.Num.AbsTol);


NK  = PND.Kite.N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Equilibrium                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('---------------------------------------------------------')
display('Computing Equilibrium State')
[u0  Error Flag] = Equilibrium_KT(0,PND);

display('---------------------------------------------------------')
display('Stability Analysis of the Equilibrium State')
Jnum             = Jacobian('Fun_ODE_KT',0,u0,PND);
    
[VecT ValT]   = eig(Jnum);
ValT         = diag(ValT);
figure(24)
hold on
plot(real(ValT),imag(ValT),'go')
xlabel('$Real(\lambda)$','interpreter','latex','fontsize',14)
ylabel('$Imag(\lambda)$','interpreter','latex','fontsize',14)
set(gca,'box','on','fontsize',12)

% Separate Longitudinal and Lateral-Direcional Modes
for i=1:1:4*NK
    for j=1:1:4*NK
       JLat(i,j) = Jnum(2*(i-1)+1,2*(j-1)+1);
       JLon(i,j) = Jnum(2*(i-1)+2,2*(j-1)+2);
    end
end
figure(5)
display(' *****          Longitudinal Modes   ********')
[Vec_Long Val_Long] = eig(JLon);
diag(Val_Long)    
display(' *****          Lateral Modes   ********')
[Vec_Lat Val_Lat] = eig(JLat);
 diag(Val_Lat)
    

%% Analysis of the tensions  
[xc T RBE rk vk ak euler omega omega_p rA rC Lambda TA TC FA MA W alfa beta Error ] = Fun_Post_KT(PD,PND,0,u0,Flag_Dim);
% Compute tension modulus
for i=1:1:PND.Kite.N
    Tension(i) = sqrt(squeeze(TA(1,:,i))*squeeze(TA(1,:,i))');
    Gamma(i)   = u0(4*(i-1)+2,1)*180/pi;
end
figure(5)
subplot(3,1,1)
hold on
plot(-rk(1,:),-rk(3,:),'-b')
plot(-rk(1,:),-rk(3,:),'b+')
xlabel('$-x_E\ (m)$','fontsize',12,'interpreter','latex')
ylabel('$-z_E\ (m)$','fontsize',12,'rotation',0,'interpreter','latex','fontsize',14)
set(gca,'fontsize',12,'box','on')

subplot(3,1,2)
hold on
plot([1:1:PND.Kite.N],alfa,'-b')
plot([1:1:PND.Kite.N],alfa,'b+')
xlabel('Kite number','fontsize',12)
ylabel('$\alpha\ (^o)$','fontsize',12,'rotation',0,'interpreter','latex','fontsize',14)
set(gca,'fontsize',12,'box','on')

subplot(3,1,3)
hold on
plot([1:1:PND.Kite.N],Tension,'-b')
plot([1:1:PND.Kite.N],Tension,'b+')
xlabel('Kite number','fontsize',12)
ylabel('$|T_A^+|\ (N)$','fontsize',12,'rotation',0,'interpreter','latex','fontsize',14)
set(gca,'fontsize',12,'box','on')


display('---------------------------------------------------------')
display('Take-off Simulation')
display('Be patient. The simulation can take several minutes')
xs_amp0 = zeros(8*PD.Kite.N,1);
for i=1:1:PD.Kite.N
    xs_amp0(4*(i-1)+1,1) = 80*pi/180;
    xs_amp0(4*(i-1)+2,1) = 88*pi/180;
    xs_amp0(4*(i-1)+3,1) = 65*pi/180;
    xs_amp0(4*(i-1)+4,1) = -70*pi/180;
end

[T u]   = ode45('Fun_ODE_KT',[0:.01:TF],xs_amp0,options);  
display('Postprocessing the results')

NK     = PND.Kite.N;
NT     = length(T);

T_out   = zeros(NT);
RBE     = zeros(3,3,NK,NT);
rk      = zeros(3,NK,NT);
vk      = zeros(3,NK,NT);
ak      = zeros(3,NK,NT);
euler   = zeros(3,NK,NT);
omega   = zeros(3,NK,NT);
omega_p = zeros(3,NK,NT);
rA      = zeros(2,3,NK,NT);
rC      = zeros(2,3,NK,NT);
Lambda  = zeros(2,NK,NT);
TA      = zeros(2,3,NK,NT);
TC      = zeros(2,3,NK,NT);
FA      = zeros(3,NK,NT);
MA      = zeros(3,NK,NT);
W       = zeros(3,NK,NT);
alfa    = zeros(NK,NT);
beta    = zeros(NK,NT);
Error   = zeros(NT);

for i=1:1:length(T)
     [xc T_out(i) RBE(:,:,:,i) rk(:,:,i) vk(:,:,i) ak(:,:,i) euler(:,:,i) omega(:,:,i) omega_p(:,:,i) rA(:,:,:,i) rC(:,:,:,i) Lambda(:,:,i)...
     TA(:,:,:,i) TC(:,:,:,i) FA(:,:,i) MA(:,:,i) W(:,:,i) alfa(:,i) beta(:,i) Error(i) ] = Fun_Post_KT(PD,PND,T(i),u(i,:)',Flag_Dim);
end
    
display('Click intro to see movie')
pause
hFig     = figure(23);
set(hFig, 'Position', [100 100 1000 600]);
for i=1:100:length(T)
     ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
     cla(ax1)
     hold on
     Plot_KT(T_out(i),rk(:,:,i),RBE(:,:,:,i),PND,PD,Flag_Dim,ax1);
     pause(0.001)
     set(get(ax1, 'title'), 'string','                  ')  
end
 

Flag_Plot = [1 1 1 1 1 ...   %[Kite Position, Velocity, Euler, alfa&beta, Tension at inelastic tethers,...
                 1 1 1 1 1 ...   % Control, Position of elastic tether, Tension of elastic tethers, Elongation, Moments  ]
                 1 1];             % Error Angular Velocity   
Plot_Results_KT(T_out,RBE,rk,vk,ak,euler,omega,omega_p,rA,rC,Lambda,TA,TC,FA,MA,W,alfa,beta,Error,Flag_Dim,Flag_Plot,PND)
     

