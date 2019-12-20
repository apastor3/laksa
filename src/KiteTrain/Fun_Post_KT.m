function [xc T RBE rK vK aK euler omegaK omegaK_p rA rC Lambda TA TC  FA MA W alfa beta Error ] = Fun_Post_KT(PD,PND,t,xs_amp,Flag_Dim)
        

%-------------------------------------------------------------------------------
% Project   : LAKSA                                                            %
% Authors   : Gonzalo Sanchez-Arriaga and  Jose A. Serrano-Iglesia             %
% Language  : Matlab                                                           %
% Synopsis  : Postprocess the state vector to find all the physical quantities %
% Copyright : Universidad Carlos III de Madrid, 2017. All rights reserved      %
%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Inputs: PD        -> Physical parameters of the system            %%
%         PND       -> Dimensionless parameters                     %%
%         t         -> Dimensionless time                           %% 
%         xs_amp    -> Extended state vector                        %% 
%         Flag_Dim  -> set 0 for dimensionless outputs              %%
%                     set 1 for outputs with dimensions             %%
%                                                                   %%
% Outputs: t        -> Time                                         %%
%          Kinematics quantities of the kite                        %%
%                      RBE     -> Body-Earth rotation matrix        %%
%                      rk      -> S_E-comp of the position,         %%
%                      vk       velocity and acceleration vectors   %%
%                      ak       of the center of mass of the kite   %% 
%                      euler   -> Euler angles                      %%
%                      omegaK  -> S_B-comp of the kite              %%
%                      omegaK_p  angular velocity and acceleration  %%
%          Kite forces and moments                                  %%
%                      Lambda_A -> Tether Tensions along            %%
%                                tether directions                  %%
%                      TA       -> SE-comp Forces exerted by the    %%
%                      TC          tethers linked at A+- and C+-    %%
%                      FA       -> SB-comp of the Aerodynamic force %% 
%                      MA          and torque                       %% 
%                      W        -> SE-comp kite weights             %%
%           Others                                                  %%
%                      alfa -> Angle of attack                      %%  
%                      beta -> Sideslip angle                       %%
%                      Error-> Comparison with classical mechanics  %%
%                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NK     = PND.Kite.N;
Nvar   = 4*NK; 
% Recover the state vector and its derivative
xs     = xs_amp(1:Nvar,1);
xs_p   = xs_amp(Nvar+1:2*Nvar,1);
% Evaluate the RHS
[xc RBE rK vK omegaK FA MA alfa beta Ups Ups_xs Phi Phi_xs Q DF RHS]=  Fun_ODE_Full_Output_KT(t,xs_amp,PND);

%% Time
T = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Kinematics calculations                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration and angular acceleration
xs_pp   = DF(Nvar+1:2*Nvar,1);
aK      = zeros(3,NK);
omegaK_p = zeros(3,NK);
for i=1:1:NK
    aK(:,i)       = Ups(:,:,i)*xs_pp;
    omegaK_p(:,i) = Phi(:,:,i)*xs_pp;
    for k=1:1:Nvar
        aK(:,i)       = aK(:,i)      + xs_p(k)*Ups_xs(:,:,k,i)*xs_p;
        omegaK_p(:,i) = omegaK_p(:,i) + xs_p(k)*Phi_xs(:,:,k,i)*xs_p;
    end
end
% Attachment points of the tethers 
rA  = zeros(2,3,NK);  
rC  = zeros(2,3,NK);  
for i=1:1:NK
    
    ATT_C(:,1,i)   = [PND.Tether.XC(i)   PND.Tether.YC(i)   PND.Tether.ZC(i)]';
    ATT_C(:,2,i)   = [PND.Tether.XC(i)  -PND.Tether.YC(i)   PND.Tether.ZC(i)]';
  
    ATT_A(:,1,i)   = [PND.Tether.XA(i)   PND.Tether.YA(i)   PND.Tether.ZA(i)]'; 
    ATT_A(:,2,i)   = [PND.Tether.XA(i)  -PND.Tether.YA(i)   PND.Tether.ZA(i)]'; 
     
    for k=1:1:2
       rA(k,:,i) = rK(:,i) + squeeze(RBE(:,:,i))'*squeeze(ATT_A(:,k,i));
       rC(k,:,i) = rK(:,i) + squeeze(RBE(:,:,i))'*squeeze(ATT_C(:,k,i));
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit vectors along tether tensions  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uA  = zeros(2,3,NK);
% First kite
for i=1:1:2 %Plus and Minus
  uA(i,:,1)  = -squeeze(rA(i,:,1));  % Point A 
  uA(i,:,1)  = uA(i,:,1)/sqrt(uA(i,:,1)*uA(i,:,1)');
end
% From Kite 2 to NK-1 
for i=2:1:NK
    for k=1:1:2
       uA(k,:,i)  = squeeze(rC(k,:,i-1))-squeeze(rA(k,:,i));  % Point A 
       uA(k,:,i)  = uA(k,:,i)/sqrt(uA(k,:,i)*uA(k,:,i)');
    end
end
  
%%%%%%%%%%%%%%%%%  
% Euler angles %%
%%%%%%%%%%%%%%%%%
euler            =  zeros(3,NK);
for i=1:1:NK
   euler(2,i)  = -asin(RBE(1,3,i));                                              % Theta (cabeceo)
   euler(1,i)  =  atan2(RBE(1,2,i)/cos(euler(2,i)),RBE(1,1,i)/cos(euler(2,i)));  % Psi (guiñada)
   euler(3,i)  =  atan2(RBE(2,3,i)/cos(euler(2,i)),RBE(3,3,i)/cos(euler(2,i)));  % Phi (balance)  
end

%%%%%%%%%%%%%%%%%
% Kite Weight  %%
%%%%%%%%%%%%%%%%%
for i=1:1:NK % S_E components of the Weight 
    W(:,i)  = PND.Kite.sigma(i)*[0 0 1]';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S_E components of the Aerodynamic force %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:NK
    FA_Earth(:,i) = squeeze(RBE(:,:,i))'*FA(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Find  the tension       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TC      = zeros(2,3,NK);
TC_Body = zeros(2,3,NK);
for i=NK:-1:1
  
    % Linear Momentum
    Mat    = [1 uA(1,:,i)*uA(2,:,i)';uA(1,:,i)*uA(2,:,i)' 1 ];
       
    for k=1:1:2
        b(k,1) = uA(k,:,i)*( PND.Kite.sigma(i)*aK(:,i)  -  FA_Earth(:,i)  -  W(:,i) - squeeze(TC(1,:,i))' - squeeze(TC(2,:,i))' );
    end
    % Tension along the tether's directions
    Sol            = Mat\b;
    Lambda(1,i)    = Sol(1);
    Lambda(2,i)    = Sol(2);
    % Compute tension components in the Earth frame
    for k=1:1:2
        TA(k,:,i)      = Sol(k)*uA(k,:,i)';
        TA_Body(k,:,i) = squeeze(RBE(:,:,i))*squeeze(TA(k,:,i))';
        if i>1
           TC(k,:,i-1)      = -TA(k,:,i);
           TC_Body(k,:,i-1) = squeeze(RBE(:,:,i-1))*squeeze(TC(k,:,i-1))';
        end
    end
 
    Error_Clas_Linear(i) = max(abs(PND.Kite.sigma(i)*aK(:,i) -  FA_Earth(:,i)  -  W(:,i)...
                             - TA(1,:,i)' - TA(2,:,i)' - TC(1,:,i)' - TC(2,:,i)'  ));
    % Angular Momentum
    for k=1:1:2
        Torque_A(k,:,i) = cross(squeeze(ATT_A(:,k,i))',squeeze(TA_Body(k,:,i)));
        Torque_C(k,:,i) = cross(squeeze(ATT_C(:,k,i))',squeeze(TC_Body(k,:,i)));
    end
    Hg = squeeze(PND.Kite.iota(:,:,i))*omegaK(:,i);

    Error_Clas_Angul(i) = max(abs( squeeze(PND.Kite.iota(:,:,i))*omegaK_p(:,i) + cross(omegaK(:,i),Hg)...
                                     - MA(:,i) - squeeze(Torque_A(1,:,i))' - squeeze(Torque_A(2,:,i))'...
                                                 - squeeze(Torque_C(1,:,i))' - squeeze(Torque_C(2,:,i))' ));
end

Error = max(max([Error_Clas_Angul Error_Clas_Linear]));


if Flag_Dim==1  % Put the outputs with dimensions
   
    % Characteristic Values
    M_ast    = PD.Kite.m(1);        % Mass
    L_ast    = PD.Tether.L(1);                     % Distance
    V_ast    = sqrt(L_ast*PD.Env.g);      % Velocity
    A_ast    = PD.Env.g;                         % Acceleration
    T_ast    = sqrt(L_ast/PD.Env.g);      % Time 
    F_ast    = M_ast*PD.Env.g;               % Force
    M_ast    = M_ast*PD.Env.g*L_ast;  % Torque
      
    T         = T_ast*t;
    rK        = L_ast*rK; 
    vK        = V_ast*vK;
    aK        = A_ast*aK;
    
    euler     = (180/pi)*euler;
    omegaK    = omegaK/T_ast;
    omegaK_p  = omegaK_p/T_ast^2;
    
    rA        = L_ast*rA;
    rC        = L_ast*rC;
    
    Lambda    = F_ast*Lambda;
    TA        = F_ast*TA;
    TC        = F_ast*TC;
    FA        = F_ast*FA;
    MA        = M_ast*MA;
    W         = F_ast*W;
  
    alfa      =  alfa*180/pi;
    beta      =  beta*180/pi;
  
    xc        =  xc*180/pi;
   
     
end





end