function  [rR vR omegaR rK vK omegaK rG vG omegaG]  = Compute_Kinematics_KF(xs,xc,xs_p,xc_p,R_KE,PND)
      
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute kinematic quantities                                   %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%----------------------------------------------------------------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                  %%
%      xs     -> state vector                                             %%
%      xc     -> control vector                                           %%
%      xs_p   -> dx_s/d tau                                               %%
%      xc_p   -> dx_c/d tau                                               %%
%      R_KE   -> SE-SK rotation matrix                                    %%
%      PND    -> Dimensionless parameters                                 %%
% Outputs                                                                 %%
%      rR     -> Rod Position vector (SE components)                      %%
%      vR     -> Rod velocity vector (SE components)                      %%
%      omegaR -> Rod angular velocity (SR components)                     %%
%      rK     -> Kite Position vector (SE components)                     %%
%      vK     -> Kite velocity vector (SE components)                     %%
%      omegaK -> Kite angular velocity (SK components)                    %%
%      rG     -> Rotor Position vector (SE components)                    %%
%      vG     -> Rotor velocity vector (SE components)                    %%
%      omegaG -> Rotor angular velocity (SK components)                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%        Tether            %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vR     = zeros(3,NR);
rR     = zeros(3,NR);
omegaR = zeros(3,NR);
for i=1:1:NR
    for j=1:1:i
        if j==i
            g = 1/2;
        else
            g = 1;
        end
        vR(1,i) = vR(1,i) - xc_p.lt*g*cos(xs.gamma(j))*cos(xs.varphi(j));
        vR(2,i) = vR(2,i) - xc_p.lt*g*cos(xs.gamma(j))*sin(xs.varphi(j));
        vR(3,i) = vR(3,i) - xc_p.lt*g*sin(xs.gamma(j));
        
        vR(1,i) = vR(1,i) + xc.lt*g*xs_p.gamma(j)*sin(xs.gamma(j))*cos(xs.varphi(j));
        vR(2,i) = vR(2,i) + xc.lt*g*xs_p.gamma(j)*sin(xs.gamma(j))*sin(xs.varphi(j));
        vR(3,i) = vR(3,i) - xc.lt*g*xs_p.gamma(j)*cos(xs.gamma(j));
        
        vR(1,i) = vR(1,i) + xc.lt*g*xs_p.varphi(j)*cos(xs.gamma(j))*sin(xs.varphi(j));
        vR(2,i) = vR(2,i) - xc.lt*g*xs_p.varphi(j)*cos(xs.gamma(j))*cos(xs.varphi(j));
        
        rR(1,i) = rR(1,i) - xc.lt*g*cos(xs.gamma(j))*cos(xs.varphi(j));
        rR(2,i) = rR(2,i) - xc.lt*g*cos(xs.gamma(j))*sin(xs.varphi(j));
        rR(3,i) = rR(3,i) - xc.lt*g*sin(xs.gamma(j));
        
    end
    omegaR(1,i) =  xs_p.varphi(i)*sin(xs.gamma(i));
    omegaR(2,i) = -xs_p.gamma(i);
    omegaR(3,i) =  xs_p.varphi(i)*cos(xs.gamma(i));
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%        KITE              %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

%% Kite Calculations
rK     = zeros(3,1);
vK     = zeros(3,1);
omegaK = zeros(3,1);
for j=1:1:NR
  rK(1,1) = rK(1,1)-xc.lt*cos(xs.gamma(j))*cos(xs.varphi(j));
  rK(2,1) = rK(2,1)-xc.lt*cos(xs.gamma(j))*sin(xs.varphi(j));
  rK(3,1) = rK(3,1)-xc.lt*sin(xs.gamma(j));

  % gamma_dot
  vK(1,1) = vK(1,1)+xc.lt*xs_p.gamma(j)*sin(xs.gamma(j))*cos(xs.varphi(j));
  vK(2,1) = vK(2,1)+xc.lt*xs_p.gamma(j)*sin(xs.gamma(j))*sin(xs.varphi(j));
  vK(3,1) = vK(3,1)-xc.lt*xs_p.gamma(j)*cos(xs.gamma(j));
 
  % varphi_dot 
  vK(1,1) = vK(1,1)+xc.lt*xs_p.varphi(j)*cos(xs.gamma(j))*sin(xs.varphi(j));
  vK(2,1) = vK(2,1)-xc.lt*xs_p.varphi(j)*cos(xs.gamma(j))*cos(xs.varphi(j));
 
  % lt_t
  vK(1,1) = vK(1,1)-xc_p.lt*cos(xs.gamma(j))*cos(xs.varphi(j));
  vK(2,1) = vK(2,1)-xc_p.lt*cos(xs.gamma(j))*sin(xs.varphi(j));
  vK(3,1) = vK(3,1)-xc_p.lt*sin(xs.gamma(j));
  
end

% Complete the position vector with the QG contribution
QG = -xc.lb*[cos(xc.delta)*cos(xc.eta) cos(xc.delta)*sin(xc.eta) sin(xc.delta)]';
rK = rK + R_KE'*QG; 

% Velocity vector contribution due to QG 
% Components in the body frame
v_QG(1,1)  = -xc_p.lb*cos(xc.delta)*cos(xc.eta) + xc.lb*xc_p.delta*sin(xc.delta)*cos(xc.eta) + xc.lb*xc_p.eta*cos(xc.delta)*sin(xc.eta);
v_QG(2,1)  = -xc_p.lb*cos(xc.delta)*sin(xc.eta) + xc.lb*xc_p.delta*sin(xc.delta)*sin(xc.eta) - xc.lb*xc_p.eta*cos(xc.delta)*cos(xc.eta);
v_QG(3,1)  = -xc_p.lb*sin(xc.delta)             - xc.lb*xc_p.delta*cos(xc.delta);
% Components in the Earth Frame
v_QG       = R_KE'*v_QG;


% QG components in the Earth frame

% Angular velocity componentsin the Kite frame
omegaK(1,1) = xs_p.phi               - xs_p.psi*sin(xs.theta);
omegaK(2,1) = xs_p.theta*cos(xs.phi) + xs_p.psi*cos(xs.theta)*sin(xs.phi);
omegaK(3,1) = xs_p.psi*cos(xs.theta)*cos(xs.phi)-xs_p.theta*sin(xs.phi);

% Kite velocity Components in the Earth frame
vK         = vK+v_QG+R_KE'*cross(omegaK,QG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%        Generators        %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
rG     = zeros(3,NG);
vG     = zeros(3,NG);
omegaG = zeros(3,NG);
if NG>0
    nu     = PND.Gen.nu;
    for i=1:1:NG
       OKOG        = [PND.Gen.x(i)  PND.Gen.y(i)  PND.Gen.z(i) ]';

       rG(:,i)     = rK + R_KE'*OKOG;                      % Components in the Earth frame
       vG(:,i)     = vK + R_KE'*cross(omegaK,OKOG);        % Components in the Earth frame 
       omegaG(:,i) = omegaK + xs_p.lambda(i)*[cos(nu(i)), 0, -sin(nu(i))]'; % Component in the kite frame 
    end
end

end
