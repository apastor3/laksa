function [SR CR OmR SK CK OmK SG CG OmG] = Matrix_V_KF(xs,xc,R_KE,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Find matrices to compute the velocity                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                                    %%
%       State vectors   xs -> [xs.gamma xs.varxs.phi thxc.eta xs.psi xs.phi]%% 
%       Control vector  xc -> [xc.lt xc.lb dexc.lta xc.eta]                 %%
%       Rotation matrix R_KE                                                %%
%       Dimensionless parameters -> PND                                     %%
% Output:                                                                   %%  
%                                                                           %%
%       Rod  velocity in SE                    -> v_R     =SR*xs_p+CR*xc_p  %%
%       Rod  angular velocity in the Rod frame -> omega_R = OmR*xs_p        %%
%       Kite velocity in the SE                -> v_k     = SK*xs_p+CK*xc_p %%
%       Kite angular velocity in SK            -> omega_K = OmK*xs_p        %%
%       Generator velocity in SE               -> v_G     = SG*xs_p+CG*xc_p %%
%       Generator angular velocity in SK       -> omega_R = OmR*xs_p        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Number of bars
NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators

Nc0 = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%              Tether      %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SR  = zeros(3,2*NR+3+NG,NR); % Initialize the matrix
CR  = zeros(3,Nc0,NR);       % Initialize the matrix
OmR = zeros(3,2*NR+3+NG,NR); % Initialize the matrix
for i=1:1:NR
    for j=1:1:i
        if j==i
            g = 0.5;
        else
            g = 1;
        end
        SR(1,j,i)    =  xc.lt*g*sin(xs.gamma(j))*cos(xs.varphi(j));
        SR(2,j,i)    =  xc.lt*g*sin(xs.gamma(j))*sin(xs.varphi(j));
        SR(3,j,i)    = -xc.lt*g*cos(xs.gamma(j));
      
        SR(1,NR+j,i) =  xc.lt*g*cos(xs.gamma(j))*sin(xs.varphi(j));
        SR(2,NR+j,i) = -xc.lt*g*cos(xs.gamma(j))*cos(xs.varphi(j));
      
        CR(1,1,i)    =  CR(1,1,i)-g*cos(xs.gamma(j))*cos(xs.varphi(j));  
        CR(2,1,i)    =  CR(2,1,i)-g*cos(xs.gamma(j))*sin(xs.varphi(j));  
        CR(3,1,i)    =  CR(3,1,i)-g*sin(xs.gamma(j));  
    end
    OmR(1,NR+i,i) = sin(xs.gamma(i));
    OmR(2,i,i)    = -1;
    OmR(3,NR+i,i) = cos(xs.gamma(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%         Kite      %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SK = zeros(3,2*NR+3+NG); % Initialize the matrix
CK = zeros(3,Nc0);      % Initialize the matrix
for j=1:1:NR
      
        SK(1,j)    =  xc.lt*sin(xs.gamma(j))*cos(xs.varphi(j));
        SK(2,j)    =  xc.lt*sin(xs.gamma(j))*sin(xs.varphi(j));
        SK(3,j)    = -xc.lt*cos(xs.gamma(j));
    
        SK(1,NR+j) =  xc.lt*cos(xs.gamma(j))*sin(xs.varphi(j));
        SK(2,NR+j) = -xc.lt*cos(xs.gamma(j))*cos(xs.varphi(j));
        
       
        CK(1,1)   =  CK(1,1)-cos(xs.gamma(j))*cos(xs.varphi(j));  
        CK(2,1)   =  CK(2,1)-cos(xs.gamma(j))*sin(xs.varphi(j));  
        CK(3,1)   =  CK(3,1)-sin(xs.gamma(j));  
end

% d QG/ dtau
AUX(1,:)    = [-cos(xc.delta)*cos(xc.eta) xc.lb*sin(xc.delta)*cos(xc.eta)  xc.lb*cos(xc.delta)*sin(xc.eta)]; 
AUX(2,:)    = [-cos(xc.delta)*sin(xc.eta) xc.lb*sin(xc.delta)*sin(xc.eta) -xc.lb*cos(xc.delta)*cos(xc.eta)]; 
AUX(3,:)    = [-sin(xc.delta)            -xc.lb*cos(xc.delta)                0]; 

CK(:,2:4)   = R_KE'*AUX;

% Term omega cross QG
AUX0      = zeros(3,3);

AUX0(1,1) = sin(xc.delta)*cos(xs.phi)+cos(xc.delta)*sin(xs.phi)*sin(xc.eta);
AUX0(1,2) = sin(xc.delta)*cos(xs.theta)*sin(xs.phi)-cos(xc.delta)*sin(xc.eta)*cos(xs.phi)*cos(xs.theta);

AUX0(2,1) = -cos(xc.delta)*sin(xs.phi)*cos(xc.eta);
AUX0(2,2) = cos(xc.delta)*cos(xc.eta)*cos(xs.phi)*cos(xs.theta)+sin(xc.delta)*sin(xs.theta);
AUX0(2,3) = -sin(xc.delta);

AUX0(3,1)  = -cos(xc.delta)*cos(xc.eta)*cos(xs.phi);
AUX0(3,2)  = -cos(xc.delta)*sin(xc.eta)*sin(xs.theta)-cos(xc.delta)*cos(xc.eta)*cos(xs.theta)*sin(xs.phi);
AUX0(3,3)  =  cos(xc.delta)*sin(xc.eta);

SK(:,2*NR+1:2*NR+3) = -xc.lb*R_KE'*AUX0;

%% Kite Angular velocity Matrix (Kite frame)
OmK            = zeros(3,2*NR+3+NG);

OmK(1,2*NR+2)   = -sin(xs.theta);
OmK(1,2*NR+3)   = 1;

OmK(2,2*NR+1)   = cos(xs.phi);
OmK(2,2*NR+2)   = cos(xs.theta)*sin(xs.phi);

OmK(3,2*NR+1)   = -sin(xs.phi);
OmK(3,2*NR+2)   = cos(xs.theta)*cos(xs.phi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%         Generators            %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SG   = zeros(3,2*NR+3+NG,NG);  % Initialize the matrix
CG   = zeros(3,Nc0,NG);        % Initialize the matrix
OmG  = zeros(3,2*NR+3+NG,NR);  % Initialize the matrix
if NG>0
   nu   = PND.Gen.nu;
    for i=1:1:NG

        SG(:,:,i) = SK;  
        SG0       = zeros(3,3);

        SG0(1,1)  = PND.Gen.z(i)*cos(xs.phi)+PND.Gen.y(i)*sin(xs.phi);
        SG0(1,2)  = PND.Gen.z(i)*cos(xs.theta)*sin(xs.phi)-PND.Gen.y(i)*cos(xs.theta)*cos(xs.phi);

        SG0(2,1)  = -PND.Gen.x(i)*sin(xs.phi);
        SG0(2,2)  =  PND.Gen.x(i)*cos(xs.theta)*cos(xs.phi)+PND.Gen.z(i)*sin(xs.theta);
        SG0(2,3)  = -PND.Gen.z(i);

        SG0(3,1)  = -PND.Gen.x(i)*cos(xs.phi);
        SG0(3,2)  = -PND.Gen.y(i)*sin(xs.theta)-PND.Gen.x(i)*cos(xs.theta)*sin(xs.phi);
        SG0(3,3)  =  PND.Gen.y(i);

        SG(:,2*NR+1:2*NR+3,i) = SG(:,2*NR+1:2*NR+3,i) + R_KE'*SG0; 

        CG(:,:,i)   = CK;    

        % Matrix OmG (Kite frame)
        OmG(:,:,i)        =  OmK;
        OmG(1,2*NR+3+i,i) =  cos(nu(i));
        OmG(3,2*NR+3+i,i) = -sin(nu(i));
    end
end


