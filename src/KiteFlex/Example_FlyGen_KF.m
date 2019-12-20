%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Example for a Fly-Generation system                            %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %%
% Inputs: No inputs                                                       %%
%                                                                         %%
% Outputs: Integration Results are placed in the workspace                %% 
%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

global PND

% Add the path of the common folder
addpath('../common/')


Validation_Type = 3;    % 1 -> Gradients and conservation of the energy
                        % 2 -> Check Hamilton equations
                        % 3 -> Equilibrium
                        
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if Validation_Type == 1
    % Load Kite parameters with dimensions
    PD = Fun_PD_FlyGen_KF;
    % Construct the dimensionless parameters
    PND             = Fun_PND_KF(PD);
    % Recover some parameters
    NT               = 2*PND.Num.N+3;                   % Rods and euler angles
    NT_p             = 2*PND.Num.N+3+PND.Gen.Num;       % Rods, euler angles, and generator
    NC               = 4+PND.Gen.Num+3;                 % control variables
    Ns               = NT+NT_p;                         % Dimension of the full state vector
    % Create artifical state vector
    xs_amp          = rand(Ns,1);
    % Create artificial control vector
    xc_amp          = rand(3*NC,1);
    % Create artificial time
    t               = rand(1);
    % Step foor finite differences
    dh              = 1e-6;
    % State and control vectors 
    [xs xs_p]       = From_xs2Var_KF(xs_amp,PND);
    [xc xc_p xc_pp] = From_xc2Var_KF(xc_amp,PND);
    % Options for the integrator
    %options = odeset('RelTol',PND.Num.IntTol,'AbsTol',PND.Num.IntTol);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check Velocity and angular velocity   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    display('*************************************** ')
    display('--------------------------------------- ')
    display(' Velocities and angular Velocities      ')
    display('--------------------------------------- ')
    display('-')
    [R_KE Grad_R_KE] = Matrix_R_KE_KF(xs);

    [SR CR OmR SK CK OmK SG CG OmG] = Matrix_V_KF(xs,xc,R_KE,PND);

    
    [rQ rR0 rR_Edge0 rK0 rG0 ]     = Compute_Positions_KF(xs,xc,R_KE,PND);
   
    [rR vR omega_R rK vK omega_K rG vG omegaG]  = Compute_Kinematics_KF(xs,xc,xs_p,xc_p,R_KE,PND);
   
    for i=1:1:PND.Num.N
        Error_rR(i) =  max(abs( rR0(:,i)- rR(:,i)));
        vR0         = squeeze(SR(:,:,i))*xs_amp(NT+1:Ns,1) +  CR(:,:,i)*xc_amp(NC+1:NC+4,1);
        Error_VR(i) =  max(abs( vR0- vR(:,i))); 
        omegaR0     = squeeze(OmR(:,:,i))*xs_amp(NT+1:Ns,1);
        Error_OR(i) =  max(abs( omegaR0- omega_R(:,i)));
    end
    
    for i=1:1:PND.Gen.Num
        Error_rG(i) = max(abs( rG0(:,i)- rG(:,i)));
        vG0         = squeeze(SG(:,:,i))*xs_amp(NT+1:Ns,1) +  CG(:,:,i)*xc_amp(NC+1:NC+4,1);
        Error_VG(i) = max(abs( vG0- vG(:,i)));
        omegaG0     = squeeze(OmG(:,:,i))*xs_amp(NT+1:Ns,1);      
        Error_OG(i) = max(abs( omegaG0- omegaG(:,i)));
    end
  
    
    display(['Error in rR = ' num2str(max(Error_rR)) ])
    display(['Error in vR = ' num2str(max(Error_VR)) ])
    display(['Error in omegaR = ' num2str(max(Error_OR)) ])
    display('-')
    display(['Error in rk = ' num2str( max(abs( rK0- rK)) ) ])    
    vK0 = SK(:,:)*xs_amp(NT+1:Ns,1) +  CK(:,:)*xc_amp(NC+1:NC+4,1);
    display(['Error en Vk = ' num2str( max(abs( vK0- vK)) ) ])  
    omega_K0      = OmK*xs_amp(NT+1:Ns,1);
    display(['Error en Omegak = ' num2str( max(abs( omega_K0- omega_K)) ) ])
    display('-')
    if PND.Gen.Num>0
        display(['Error in rG = ' num2str(max(Error_rG)) ])
        display(['Error in vG = ' num2str(max(Error_VG)) ])
        display(['Error in OmegaG = ' num2str(max(Error_OG)) ])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check Bar Gradients using finite-differences  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display('--------------------------------------- ')
    display('Bar Gradients using finite-differences ')
    
    [SR_xs SR_xc CR_xs OmR_xs] = Grad_Rod_KF(xs,xc,PND);

    for j=1:1:NT % Do Jacobian with respect to xs
          xs_amp_Mas                    = xs_amp;
          xs_amp_Mas(j,1)               = xs_amp_Mas(j,1)+dh;
          [xs_Mas xs_p_Mas]             = From_xs2Var_KF(xs_amp_Mas,PND);
          
          [R_KE_Mas Grad_R_KE_Mas]      = Matrix_R_KE_KF(xs_Mas); 
          [SR0 CR0 OmR0 SK0 CK0 OmK0 SG0 CG0 OmG0] = Matrix_V_KF(xs_Mas,xc,R_KE_Mas,PND);  
      
          SR_xs0(:,:,:,j)               = (SR0-SR)/dh; 
          CR_xs0(:,:,:,j)               = (CR0-CR)/dh; 
          OmR_xs0(:,:,:,j)              = (OmR0-OmR)/dh; 
    end   
    Error_SR_xs  = max(max(max(max(abs(SR_xs-SR_xs0)))));
    Error_CR_xs  = max(max(max(max(abs(CR_xs-CR_xs0)))));
    Error_OmR_xs = max(max(max(max(abs(OmR_xs-OmR_xs0)))));

    for j=1:1:4 % Do Jacobian with respect to xc
          xc_amp_Mas                    = xc_amp;
          xc_amp_Mas(j,1)               = xc_amp_Mas(j,1)+dh;
         [xc_Mas xc_p_Mas]             = From_xc2Var_KF(xc_amp_Mas,PND);
         [SR0 CR0 OmR0 SK0 CK0 OmK0 SG0 CG0 OmG0] = Matrix_V_KF(xs,xc_Mas,R_KE,PND); 
         SR_xc0(:,:,:,j)               = (SR0-SR)/dh; 
    end
    Error_SR_xc = max(max(max(max(abs(SR_xc-SR_xc0)))));

    display(['Error en SR_xs  = ' num2str( max(Error_SR_xs )) ])
    display(['Error en SR_xc  = ' num2str( max(Error_SR_xc )) ])
    display(['Error en CR_xs  = ' num2str( max(Error_CR_xs )) ])
    display(['Error en OmR_xs  = ' num2str( max(Error_OmR_xs )) ])

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check Kite Gradients using finite-differences  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    display('--------------------------------------- ')
    display('Kite Gradients using finite-differences ')
    [SK_xs SK_xc CK_xs CK_xc  OmK_xs] = Grad_Kite_KF(xs,xc,R_KE,Grad_R_KE,PND);
    for j=1:1:NT % Do Jacobian with respect to xs
          xs_amp_Mas                    = xs_amp;
          xs_amp_Mas(j,1)               = xs_amp_Mas(j,1)+dh;
          [xs_Mas xs_p_Mas]             = From_xs2Var_KF(xs_amp_Mas,PND);        
          [R_KE_Mas Grad_R_KE_Mas ]     = Matrix_R_KE_KF(xs_Mas);
          [SR0 CR0 OmR0 SK0 CK0 OmK0 SG0 CG0 OmG0] = Matrix_V_KF(xs_Mas,xc,R_KE_Mas,PND);
          SK_xs0(:,:,j)                 = (SK0-SK)/dh; 
          CK_xs0(:,:,j)                 = (CK0-CK)/dh; 
          OmK_xs0(:,:,j)                = (OmK0-OmK)/dh; 
    end

    Error_SK_xs = max(max(max(max(abs(SK_xs-SK_xs0)))));
    Error_CK_xs = max(max(max(max(abs(CK_xs-CK_xs0)))));
    Error_OmK_xs = max(max(max(max(abs(OmK_xs-OmK_xs0)))));


    for j=1:1:4 % Do Jacobian with respect to xc
          xc_amp_Mas                    = xc_amp;
          xc_amp_Mas(j,1)               = xc_amp_Mas(j,1)+dh;
          [xc_Mas xc_p_Mas]             = From_xc2Var_KF(xc_amp_Mas,PND);
          [SR0 CR0 OmR0 SK0 CK0 OmK0 SG0 CG0 OmG0] = Matrix_V_KF(xs,xc_Mas,R_KE,PND);    
          SK_xc0(:,:,j)                 = (SK0-SK)/dh; 
          CK_xc0(:,:,j)                 = (CK0-CK)/dh; 
    end
    Error_SK_xc = max(max(max(max(abs(SK_xc-SK_xc0)))));
    Error_CK_xc = max(max(max(max(abs(CK_xc-CK_xc0)))));

    display(['Error en Sk_xs  = ' num2str( max(Error_SK_xs )) ])
    display(['Error en Sk_xc  = ' num2str( max(Error_SK_xc )) ])
    display(['Error en Ck_xs  = ' num2str( max(Error_CK_xs )) ])
    display(['Error en Ck_xc  = ' num2str( max(Error_CK_xc )) ])
    display(['Error en Omk_xs  = ' num2str( max(Error_OmK_xs )) ])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check Generator Gradients using finite-differences  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    display('--------------------------------------- ')
    display('Generators Gradients using finite-differences ')
    [SG_xs SG_xc CG_xs CG_xc ] = Grad_Rotors_KF(xs,xc,R_KE,Grad_R_KE,SK_xs,SK_xc,CK_xs,CK_xc,PND);
   
    for j=1:1:NT % Do Jacobian with respect to xs
          xs_amp_Mas                    = xs_amp;
          xs_amp_Mas(j,1)               = xs_amp_Mas(j,1)+dh;
          [xs_Mas xs_p_Mas]             = From_xs2Var_KF(xs_amp_Mas,PND);  
          
          [R_KE_Mas Grad_R_KE_Mas ]                = Matrix_R_KE_KF(xs_Mas);
          [SR0 CR0 OmR0 SK0 CK0 OmK0 SG0 CG0 OmG0] = Matrix_V_KF(xs_Mas,xc,R_KE_Mas,PND);   
            
          SG_xs0(:,:,:,j)                 = (SG0-SG)/dh; 
          CG_xs0(:,:,:,j)                 = (CG0-CG)/dh; 
    end

    Error_SG_xs = max(max(max(max(abs(SG_xs-SG_xs0)))));
    Error_CG_xs = max(max(max(max(abs(CG_xs-CG_xs0)))));
   

    for j=1:1:4 % Do Jacobian with respect to xc
          xc_amp_Mas                    = xc_amp;
          xc_amp_Mas(j,1)               = xc_amp_Mas(j,1)+dh;
          [xc_Mas xc_p_Mas]             = From_xc2Var_KF(xc_amp_Mas,PND);
          [SR0 CR0 OmR0 SK0 CK0 OmK0 SG0 CG0 OmG0] = Matrix_V_KF(xs,xc_Mas,R_KE,PND);    
          SG_xc0(:,:,:,j)               = (SG0-SG)/dh; 
          CG_xc0(:,:,:,j)               = (CG0-CG)/dh; 
    end
    Error_SG_xc = max(max(max(max(abs(SG_xc-SG_xc0)))));
    Error_CG_xc = max(max(max(max(abs(CG_xc-CG_xc0)))));

    display(['Error en SG_xs  = ' num2str( max(Error_SG_xs )) ])
    display(['Error en SG_xc  = ' num2str( max(Error_SG_xc )) ])
    display(['Error en CG_xs  = ' num2str( max(Error_CG_xs )) ])
    display(['Error en CG_xc  = ' num2str( max(Error_CG_xc )) ])
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Check Matrices  and their gradients       %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display('--------------------------------------- ')
    display('Matrices Gradients using finite-differences ')
    [Ms Msc Mc] = Matrix_M_KF(SK,CK,OmK,SR,CR,OmR,SG,CG,xs,xc,PND);
    [Ms_xs Ms_xc Msc_xs Msc_xc Mc_xs Mc_xc] = Grad_Matrix_M_KF(R_KE,Grad_R_KE,SR,CR,OmR,SK,CK,OmK,SG,CG,xs,xc,1,PND);
         
    for j=1:1: NT % Do Jacobian with respect to xs
          xs_amp_Mas                          = xs_amp;
          xs_amp_Mas(j,1)                     = xs_amp_Mas(j,1)+dh;
          [xs_Mas xs_p_Mas]                   = From_xs2Var_KF(xs_amp_Mas,PND);

          [R_KE_M Grad_R_KE_M ] = Matrix_R_KE_KF(xs_Mas);
          [SR_M CR_M OmR_M SK_M CK_M OmK_M SG_M CG_M OmG_M]= Matrix_V_KF(xs_Mas,xc,R_KE_M,PND); 
        
          [Ms_Mas Msc_Mas Mc_Mas]             = Matrix_M_KF(SK_M,CK_M,OmK_M,SR_M,CR_M,OmR_M,SG_M,CG_M,xs_Mas,xc,PND);
         
          Ms_xs0(:,:,j)                       = (Ms_Mas-Ms)/dh; 
          Msc_xs0(:,:,j)                      = (Msc_Mas-Msc)/dh; 
          Mc_xs0(:,:,j)                       = (Mc_Mas-Mc)/dh; 
    end
 
    Error_Ms_xs = max(max(max(max(abs(Ms_xs-Ms_xs0)))));
    Error_Msc_xs = max(max(max(max(abs(Msc_xs-Msc_xs0)))));
    Error_Mc_xs = max(max(max(max(abs(Mc_xs-Mc_xs0)))));

    display(['Error en Ms_xs  = ' num2str( max(Error_Ms_xs )) ])
    display(['Error en Msc_xs  = ' num2str( max(Error_Msc_xs )) ])
    display(['Error en Mc_xs  = ' num2str( max(Error_Mc_xs )) ])

    for j=1:1:4 % Do Jacobian with respect to xc
          xc_amp_Mas                          = xc_amp;
          xc_amp_Mas(j,1)                     = xc_amp_Mas(j,1)+dh;
          [xc_Mas xc_p_Mas]                   = From_xc2Var_KF(xc_amp_Mas,PND);
      
          [R_KE Grad_R_KE ]                   = Matrix_R_KE_KF(xs);
      
          [SR_M CR_M OmR_M SK_M CK_M OmK_M SG_M CG_M OmG_M] = Matrix_V_KF(xs,xc_Mas,R_KE,PND);   
          [Ms_Mas Msc_Mas Mc_Mas]                      = Matrix_M_KF(SK_M,CK_M,OmK_M,SR_M,CR_M,OmR_M,SG_M,CG_M,xs,xc_Mas,PND);
          
          Ms_xc0(:,:,j)                       = (Ms_Mas-Ms)/dh;  
          Msc_xc0(:,:,j)                      = (Msc_Mas-Msc)/dh;
          Mc_xc0(:,:,j)                       = (Mc_Mas-Mc)/dh;
    end
    Error_Ms_xc = max(max(max(max(abs(Ms_xc-Ms_xc0)))));
    Error_Msc_xc = max(max(max(max(abs(Msc_xc-Msc_xc0)))));
    Error_Mc_xc  = max(max(max(max(abs(Mc_xc-Mc_xc0)))));

    display(['Error en Ms_xc  = ' num2str( max(Error_Ms_xc )) ])
    display(['Error en Msc_xc  = ' num2str( max(Error_Msc_xc )) ])
    display(['Error en Mc_xc  = ' num2str( max(Error_Mc_xc )) ])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Energy                                    %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display('--------------------------------------- ')
    display('Energy Balance ')
    [Em Ec Ep H] = Compute_Energy_KF(SR,CR,OmR,SK,CK,OmK,SG,CG,xs_amp,xc_amp,PND);
    U_xs = Grad_U_KF(xs,xc,PND);
   
    for j=1:1:NT % Do Jacobian with respect to xs
          xs_amp_Mas                          = xs_amp;
          xs_amp_Mas(j,1)                     = xs_amp_Mas(j,1)+dh;
         [xs_Mas xs_p_Mas]                    = From_xs2Var_KF(xs_amp_Mas,PND);
         
         [R_KE_M Grad_R_KE_M ]                = Matrix_R_KE_KF(xs_Mas);
         [SR_M CR_M OmR_M SK_M CK_M OmK_M SG_M CG_M OmG_M] = Matrix_V_KF(xs_Mas,xc,R_KE_M,PND); 
         
          [Em_Mas Ec_Mas Ep_Mas H_Mas]  = Compute_Energy_KF(SR_M,CR_M,OmR_M,SK_M,CK_M,OmK_M,SG_M,CG_M,xs_amp_Mas,xc_amp,PND);
          U_xs0(j,1)                    = (Ep_Mas-Ep)/dh; 
    end
    Error_U_xs = max(max(max(max(abs(U_xs-U_xs0)))));
    display(['Error en U_xs  = ' num2str( max(Error_U_xs )) ])

end



if Validation_Type == 2
          
    % Compute equilibrium with the number of bars set in Fun_PD_KF_FlyGen
    PD = Fun_PD_FlyGen_KF;
    PND = Fun_PND_KF(PD);
    
    [X  Error Flag PND delta_a_Eq xi_Rotor_Eq]=Equilibrium_FlyGen_KF(0,PND);
       
    
    %% Load the parameters
    
    options = odeset('RelTol',PND.Num.IntTol,'AbsTol',PND.Num.IntTol','Stats','on');
    Nvar    = 2*PND.Num.N+3;               % Number of variables 
    Nvar_p  = 2*PND.Num.N+3+PND.Gen.Num;   % Number of variables (dot) 
    Nc      = 4+PND.Gen.Num ;              % Number of control parameters
    xc_amp                       = Fun_Control_KF(0,X,PND);
            
    [xs xs_p]                    = From_xs2Var_KF(X,PND);
    [xc xc_p xc_pp]              = From_xc2Var_KF(xc_amp,PND);
   
    xs_amp = X;
  
   [R_KE Grad_R_KE]             = Matrix_R_KE_KF(xs);
   [SR CR OmR SK CK OmK SG CG ] = Matrix_V_KF(xs,xc,R_KE,PND);
  
   [Em Ec Ep H] = Compute_Energy_KF(SR,CR,OmR,SK,CK,OmK,SG,CG,xs_amp,xc_amp,PND);
 
    X0_Lagrange = [xs_amp;H];
    TF = 1.0;
   
    %------------------------- Integration  ---------------------------------
    display('Lagrangian Integration')
    % Lagrangian integraiton 
    tic
    [T_Lan Y_Lan]   = ode45('Fun_ODE_Lag_KF',[0  TF],X0_Lagrange,options);
     Time_Lag =  toc
  
    % Hamiltonian integration
    %Hamiltonian initial condition
    X0_Hamiltonian   = Lag2Ham_KF(0,X0_Lagrange,PND);
    X0_Hamiltonian  =  [X0_Hamiltonian;H]; 
    display('Hamiltonian Integration')
    tic
    [T_Ham Y_Ham]                = ode45('Fun_ODE_Ham_KF',T_Lan,X0_Hamiltonian,options);
    Time_Ham =  toc
  
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Recover Lagrangan variables and plot
    for i=1:1:length(T_Lan)
        Y_Lan2(i,:) = Ham2Lag_KF(T_Lan(i),Y_Ham(i,:)',PND);
        xs_amp                      = Y_Lan2(i,1:Nvar+Nvar_p)';
        xc_amp                      = Fun_Control_KF(T_Lan(i),xs_amp,PND);
       [xs xs_p]                    = From_xs2Var_KF(xs_amp,PND);
       [xc xc_p xc_pp]              = From_xc2Var_KF(xc_amp,PND);
       [R_KE Grad_R_KE]             = Matrix_R_KE_KF(xs);
       [SR CR OmR SK CK OmK SG CG OmG]  = Matrix_V_KF(xs,xc,R_KE,PND);
       [Em Ec Ep H(i)]              = Compute_Energy_KF(SR,CR,OmR,SK,CK,OmK,SG,CG,xs_amp,xc_amp,PND);
        
    end
    figure(1)
    subplot(2,1,1)
    hold on
    plot(T_Lan,Y_Lan(:,1)*180/pi,'+')
    plot(T_Ham,Y_Ham(:,1)*180/pi,'or')    
    ylabel('\Gamma')
    xlabel('\tau')
   
    subplot(2,1,2)
    hold on
    plot(T_Lan,Y_Lan(:,Nvar+1),'+')
    plot(T_Ham,Y_Lan2(:,Nvar+1),'or')
    ylabel('d\Gamma/d\tau')
    xlabel('\tau')
    
   figure(3)
   semilogy(T_Lan,abs(Y_Lan(:,end)-H'),'b+')
   hold on
   semilogy(T_Ham,abs(Y_Ham(:,end)-H'),'ro')
   legend('Lagrangian', 'Hamiltonian')
end


if Validation_Type == 3
    Flag_Dim = 1;
    % We try to find an equilibrium state with the following characteristics 
    
    % a) Symmetrically flight: varphi = 0, psi=0, phi = 0. 
    % b) Rotors spinning at the target angular velocity lambda_p = PND.Control.Ome_Rotor
    % c) The reaction torques of the motor controllers are compensated by the ailerons.
    
    % This scheme only works if nu = 0 (otherwise a rudder deflection would be also needed)
    
    % Load Kite parameters with dimensions
    PD = Fun_PD_FlyGen_KF;
    % Construct the dimensionless parameters
    PND             = Fun_PND_KF(PD);
    
    
    [u0  Error Flag PND delta_a_Eq xi_Rotor_Eq]=Equilibrium_FlyGen_KF(0,PND);  
  
    PND.Target.delta_a  = delta_a_Eq;
    PND.Target.xi_Rotor = xi_Rotor_Eq;        
    
    
    [T rR_Edge rR vR aR omegaR gR FA_R Tension rQ R_KE ...
     rK vK aK euler omegaK gK FA_K FR_K FG_K  MA_K MR_K MG_K MMC_K alfa_K beta_K... 
     rG vG aG omegaG gG FA_G  FK_G MA_G MK_G  MMC_G xc_out xs_target_out Error_C] = Fun_Post_KF(PD,PND,0,u0,Flag_Dim);
           
    display('Integrating Equations of motion')
    display('Initial condition = Equilibrium state + perturbationt')
    [T X] = ode45('Fun_ODE_Lag_KF',[0:0.01:1],u0+1e-3*rand(length(u0),1));
    
    
    
    
    for i=1:1:length(T)
          [Tout(i) rR_Edge(:,:,i) rR(:,:,i) vR(:,:,i) aR(:,:,i) omegaR(:,:,i) gR(:,:,i) FA_R(:,:,i) Tension(:,:,i) rQ(:,i) R_KE(:,:,i) ...
          rK(:,i) vK(:,i) aK(:,i) euler(:,i) omegaK(:,i) gK(:,i) FA_K(:,i) FR_K(:,i) FG_K(:,i)  MA_K(:,i) MR_K(:,i) MG_K(:,i) MMC_K(:,i) alfa_K(i) beta_K(i)... 
          rG(:,:,i) vG(:,:,i) aG(:,:,i) omegaG(:,:,i) gG(:,:,i) FA_G(:,:,i)  FK_G(:,:,i) MA_G(:,:,i) MK_G(:,:,i)  MMC_G(:,:,i) xc_out(:,i) xs_target_out(:,i) Error_C(i)] = Fun_Post_KF(PD,PND,T(i),X(i,:)',Flag_Dim);
            
         Plot_FlyGen_KF(X(i,:)',Tout(i),rQ(:,i),rR(:,:,i),rK(:,i),rG(:,:,i),R_KE(:,:,i),rR_Edge(:,:,i),PND,Flag_Dim,PD);
         pause(0.01)
         title('')
    end
   
                %[Kite Position, Velocity, Euler, alfa&beta, Tension, Control, Rotor angular velocity, Error_C
    Flag_Plot = [1             ,     1   ,   1  ,     1    ,    1   ,    1   ,       1               ,   1  ];
              
    
    Plot_Results_KF(Tout,rR_Edge, rR, vR, aR, omegaR, gR, FA_R, Tension, rQ, R_KE, ...
          rK, vK, aK, euler, omegaK, gK, FA_K, FR_K, FG_K,  MA_K, MR_K, MG_K, MMC_K, alfa_K, beta_K,... 
          rG, vG, aG, omegaG, gG, FA_G,  FK_G, MA_G, MK_G,  MMC_G, xc_out,xs_target_out,Error_C, Flag_Dim, Flag_Plot,PD)
   

   
end