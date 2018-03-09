function [SK_xs SK_xc CK_xs CK_xc OmK_xs] = Grad_Kite_KF(xs,xc,R_KE,Grad_R_KE,PND)
       
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Gradient of Kite matrices                                      %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%----------------------------------------------------------------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                                    %%
%    State vectors   xs -> [xs.gamma xs.varxs.phi thxc.eta xs.psi xs.phi ]  %%
%    Control vector  xc -> [xc.lt xc.lb dexc.lta xc.eta]                    %%
%    Rotation matrix R_KE and its gradient Grad_R_KE                        %%
%    PND                -> Dimensionless parameters                         %%
% Output:                                                                   %%
%    Gradients                                                              %%
%       Sk_xs  = partial Sk/partial x_s                                     %%
%       Sk_xc  = partial Sk/partial x_c                                     %%
%       Ck_xs  = partial Ck/partial x_s                                     %%
%       Ck_xc  = partial Ck/partial x_c                                     %%
%       Omk_xs = partial Omk/partial x_s                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Number of bars and generators

NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators

% State and control vector dimensions
Nvar   = 2*NR+3;
Nvar_p = 2*NR+3+NG;
Nc0    = 4;

SK_xs  = zeros(3,Nvar_p,Nvar);  % Initialize the matrix
SK_xc  = zeros(3,Nvar_p,Nc0);    % Initialize the matrix
CK_xs  = zeros(3,Nc0,Nvar);      % Initialize the matrix
CK_xc  = zeros(3,Nc0,Nc0);        % Initialize the matrix
OmK_xs = zeros(3,Nvar_p,Nvar);  % Initialize the matrix


%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matrix Sk
%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:1:NR
    %%%%%%%%%%%%%%%%%
    %%% Sk matrix  %%
    %%%%%%%%%%%%%%%%%
    %Derivative with respect to gamma  
    SK_xs(1,j,j)      =  xc.lt*cos(xs.gamma(j))*cos(xs.varphi(j));
    SK_xs(1,NR+j,j)   = -xc.lt*sin(xs.gamma(j))*sin(xs.varphi(j));
           
    SK_xs(2,j,j)      =  xc.lt*cos(xs.gamma(j))*sin(xs.varphi(j));
    SK_xs(2,NR+j,j)   =  xc.lt*sin(xs.gamma(j))*cos(xs.varphi(j));
        
    SK_xs(3,j,j)      =  xc.lt*sin(xs.gamma(j));
    %Derivative with respect to varphi  
    SK_xs(1,j,NR+j)   = -xc.lt*sin(xs.gamma(j))*sin(xs.varphi(j));
    SK_xs(1,NR+j,NR+j) =  xc.lt*cos(xs.gamma(j))*cos(xs.varphi(j));
           
    SK_xs(2,j,NR+j)    =  xc.lt*sin(xs.gamma(j))*cos(xs.varphi(j));
    SK_xs(2,NR+j,NR+j) =  xc.lt*cos(xs.gamma(j))*sin(xs.varphi(j));     
    
    % Derivatives with respect to the control 
    SK_xc(1,j,1)      =  sin(xs.gamma(j))*cos(xs.varphi(j));
    SK_xc(1,NR+j,1)   =  cos(xs.gamma(j))*sin(xs.varphi(j));
           
    SK_xc(2,j,1)      =  sin(xs.gamma(j))*sin(xs.varphi(j));
    SK_xc(2,NR+j,1)   = -cos(xs.gamma(j))*cos(xs.varphi(j));
      
    SK_xc(3,j,1)      = -cos(xs.gamma(j));
    
    %%%%%%%%%%%%%%%%%
    %% Ck matrix  %%%
    %%%%%%%%%%%%%%%%%
    % Derivative with respect to gamma
    CK_xs(1,1,j)     =   sin(xs.gamma(j))*cos(xs.varphi(j));  
    CK_xs(2,1,j)     =   sin(xs.gamma(j))*sin(xs.varphi(j));  
    CK_xs(3,1,j)     =  -cos(xs.gamma(j));  
    
     % Derivative with respect to varphi
    CK_xs(1,1,NR+j)   =   cos(xs.gamma(j))*sin(xs.varphi(j));  
    CK_xs(2,1,NR+j)   =  -cos(xs.gamma(j))*cos(xs.varphi(j));  
    
         
end

% Vector QG in Body Frame
QG(1,1) = -xc.lb*cos(xc.delta)*cos(xc.eta);
QG(2,1) = -xc.lb*cos(xc.delta)*sin(xc.eta);
QG(3,1) = -xc.lb*sin(xc.delta);

% Vector QG derivatiives with respect to the control vector
QG_xc   = zeros(3,4);
   % Derivative with respect to lb
    QG_xc(1,2)   = -cos(xc.delta)*cos(xc.eta);
    QG_xc(2,2)   = -cos(xc.delta)*sin(xc.eta);
    QG_xc(3,2)   = -sin(xc.delta);
   % Derivative with respect to delta
    QG_xc(1,3)   =  xc.lb*sin(xc.delta)*cos(xc.eta);
    QG_xc(2,3)   =  xc.lb*sin(xc.delta)*sin(xc.eta);
    QG_xc(3,3)   = -xc.lb*cos(xc.delta);
   % Derivative with respect to eta
    QG_xc(1,4) =  xc.lb*cos(xc.delta)*sin(xc.eta);
    QG_xc(2,4) = -xc.lb*cos(xc.delta)*cos(xc.eta);
    
% Gradient of QG in the Earth Frame
QG_xc   = R_KE'*QG_xc; %(3,4) Matrix
QG_xs   = zeros(3,3);
for i=1:1:3 % Derivative with respect to theta,psi and phi
   QG_xs(:,i) = (squeeze(Grad_R_KE(:,:,i)))'*QG;
end
% QG in the Earth Frame
QG = R_KE'*QG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Complete Sk_xs and Ck_xs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Derivatives with respect to theta, psi and phi of QG
   for i=1:1:3
        SK_xs(1,2*NR+1,2*NR+i) =   QG_xs(3,i)*cos(xs.psi);
        SK_xs(1,2*NR+2,2*NR+i) =  -QG_xs(2,i);
        SK_xs(1,2*NR+3,2*NR+i) =   QG_xs(3,i)*cos(xs.theta)*sin(xs.psi)+QG_xs(2,i)*sin(xs.theta);

        SK_xs(2,2*NR+1,2*NR+i) =   QG_xs(3,i)*sin(xs.psi);
        SK_xs(2,2*NR+2,2*NR+i) =   QG_xs(1,i);
        SK_xs(2,2*NR+3,2*NR+i) =  -QG_xs(1,i)*sin(xs.theta)-QG_xs(3,i)*cos(xs.theta)*cos(xs.psi);

        SK_xs(3,2*NR+1,2*NR+i) =  -QG_xs(2,i)*sin(xs.psi)-QG_xs(1,i)*cos(xs.psi);
        SK_xs(3,2*NR+3,2*NR+i) =   cos(xs.theta)*(QG_xs(2,i)*cos(xs.psi)-QG_xs(1,i)*sin(xs.psi));
        
        CK_xs(1,2,2*NR+i)     = -cos(xc.delta)*(Grad_R_KE(1,1,i)*cos(xc.eta)+Grad_R_KE(2,1,i)*sin(xc.eta))-Grad_R_KE(3,1,i)*sin(xc.delta);     
        CK_xs(2,2,2*NR+i)     = -cos(xc.delta)*(Grad_R_KE(1,2,i)*cos(xc.eta)+Grad_R_KE(2,2,i)*sin(xc.eta))-Grad_R_KE(3,2,i)*sin(xc.delta);     
        CK_xs(3,2,2*NR+i)     = -cos(xc.delta)*(Grad_R_KE(1,3,i)*cos(xc.eta)+Grad_R_KE(2,3,i)*sin(xc.eta))-Grad_R_KE(3,3,i)*sin(xc.delta);     

        CK_xs(1,3,2*NR+i)     = xc.lb*sin(xc.delta)*(Grad_R_KE(1,1,i)*cos(xc.eta)+Grad_R_KE(2,1,i)*sin(xc.eta))-xc.lb*Grad_R_KE(3,1,i)*cos(xc.delta);     
        CK_xs(2,3,2*NR+i)     = xc.lb*sin(xc.delta)*(Grad_R_KE(1,2,i)*cos(xc.eta)+Grad_R_KE(2,2,i)*sin(xc.eta))-xc.lb*Grad_R_KE(3,2,i)*cos(xc.delta);    
        CK_xs(3,3,2*NR+i)     = xc.lb*sin(xc.delta)*(Grad_R_KE(1,3,i)*cos(xc.eta)+Grad_R_KE(2,3,i)*sin(xc.eta))-xc.lb*Grad_R_KE(3,3,i)*cos(xc.delta);    

        CK_xs(1,4,2*NR+i)     = xc.lb*cos(xc.delta)*(Grad_R_KE(1,1,i)*sin(xc.eta)-Grad_R_KE(2,1,i)*cos(xc.eta));     
        CK_xs(2,4,2*NR+i)     = xc.lb*cos(xc.delta)*(Grad_R_KE(1,2,i)*sin(xc.eta)-Grad_R_KE(2,2,i)*cos(xc.eta));    
        CK_xs(3,4,2*NR+i)     = xc.lb*cos(xc.delta)*(Grad_R_KE(1,3,i)*sin(xc.eta)-Grad_R_KE(2,3,i)*cos(xc.eta));          
   end
   % Derivatives with respec to theta
   SK_xs(1,2*NR+3,2*NR+1) =  SK_xs(1,2*NR+3,2*NR+1) - QG(3)*sin(xs.theta)*sin(xs.psi)+QG(2)*cos(xs.theta);
   SK_xs(2,2*NR+3,2*NR+1) =  SK_xs(2,2*NR+3,2*NR+1) - QG(1)*cos(xs.theta)+QG(3)*sin(xs.theta)*cos(xs.psi);
   SK_xs(3,2*NR+3,2*NR+1) =  SK_xs(3,2*NR+3,2*NR+1) - sin(xs.theta)*(QG(2)*cos(xs.psi)-QG(1)*sin(xs.psi));
    % Derivatives with respec to psi
   SK_xs(1,2*NR+1,2*NR+2) =  SK_xs(1,2*NR+1,2*NR+2) - QG(3)*sin(xs.psi);
   SK_xs(1,2*NR+3,2*NR+2) =  SK_xs(1,2*NR+3,2*NR+2) + QG(3)*cos(xs.theta)*cos(xs.psi);
  
   SK_xs(2,2*NR+1,2*NR+2) =  SK_xs(2,2*NR+1,2*NR+2) + QG(3)*cos(xs.psi);
   SK_xs(2,2*NR+3,2*NR+2) =  SK_xs(2,2*NR+3,2*NR+2) + QG(3)*cos(xs.theta)*sin(xs.psi);

   SK_xs(3,2*NR+1,2*NR+2) =  SK_xs(3,2*NR+1,2*NR+2) - QG(2)*cos(xs.psi)+QG(1)*sin(xs.psi);
   SK_xs(3,2*NR+3,2*NR+2) =  SK_xs(3,2*NR+3,2*NR+2) - cos(xs.theta)*(QG(2)*sin(xs.psi)+QG(1)*cos(xs.psi));
  
% Complete Sk_xc
 % Derivatives with respec to theta, psi and phi of QG
   for i=2:1:4
        SK_xc(1,2*NR+1,i) =   QG_xc(3,i)*cos(xs.psi);
        SK_xc(1,2*NR+2,i) =  -QG_xc(2,i);
        SK_xc(1,2*NR+3,i) =   QG_xc(3,i)*cos(xs.theta)*sin(xs.psi)+QG_xc(2,i)*sin(xs.theta);

        SK_xc(2,2*NR+1,i) =   QG_xc(3,i)*sin(xs.psi);
        SK_xc(2,2*NR+2,i) =   QG_xc(1,i);
        SK_xc(2,2*NR+3,i) =  -QG_xc(1,i)*sin(xs.theta)-QG_xc(3,i)*cos(xs.theta)*cos(xs.psi);

        SK_xc(3,2*NR+1,i) =  -QG_xc(2,i)*sin(xs.psi)-QG_xc(1,i)*cos(xs.psi);
        SK_xc(3,2*NR+3,i) =   cos(xs.theta)*(QG_xc(2,i)*cos(xs.psi)-QG_xc(1,i)*sin(xs.psi));
   end

%%%%%%%%%%%%%%%%%%%%%   
%%%% Do Ck_xc   %%%%%
%%%%%%%%%%%%%%%%%%%%%
% Derivative with respect to lb
CK_xc(1,3,2)     = sin(xc.delta)*(R_KE(1,1)*cos(xc.eta)+R_KE(2,1)*sin(xc.eta))-R_KE(3,1)*cos(xc.delta);     
CK_xc(2,3,2)     = sin(xc.delta)*(R_KE(1,2)*cos(xc.eta)+R_KE(2,2)*sin(xc.eta))-R_KE(3,2)*cos(xc.delta);    
CK_xc(3,3,2)     = sin(xc.delta)*(R_KE(1,3)*cos(xc.eta)+R_KE(2,3)*sin(xc.eta))-R_KE(3,3)*cos(xc.delta);    

CK_xc(1,4,2)     = cos(xc.delta)*(R_KE(1,1)*sin(xc.eta)-R_KE(2,1)*cos(xc.eta));     
CK_xc(2,4,2)     = cos(xc.delta)*(R_KE(1,2)*sin(xc.eta)-R_KE(2,2)*cos(xc.eta));    
CK_xc(3,4,2)     = cos(xc.delta)*(R_KE(1,3)*sin(xc.eta)-R_KE(2,3)*cos(xc.eta));  
% Derivative with respect to delta
CK_xc(1,2,3)     = sin(xc.delta)*(R_KE(1,1)*cos(xc.eta)+R_KE(2,1)*sin(xc.eta))-R_KE(3,1)*cos(xc.delta);     
CK_xc(2,2,3)     = sin(xc.delta)*(R_KE(1,2)*cos(xc.eta)+R_KE(2,2)*sin(xc.eta))-R_KE(3,2)*cos(xc.delta);     
CK_xc(3,2,3)     = sin(xc.delta)*(R_KE(1,3)*cos(xc.eta)+R_KE(2,3)*sin(xc.eta))-R_KE(3,3)*cos(xc.delta);     

CK_xc(1,3,3)     = xc.lb*cos(xc.delta)*(R_KE(1,1)*cos(xc.eta)+R_KE(2,1)*sin(xc.eta))+xc.lb*R_KE(3,1)*sin(xc.delta);     
CK_xc(2,3,3)     = xc.lb*cos(xc.delta)*(R_KE(1,2)*cos(xc.eta)+R_KE(2,2)*sin(xc.eta))+xc.lb*R_KE(3,2)*sin(xc.delta);    
CK_xc(3,3,3)     = xc.lb*cos(xc.delta)*(R_KE(1,3)*cos(xc.eta)+R_KE(2,3)*sin(xc.eta))+xc.lb*R_KE(3,3)*sin(xc.delta);    

CK_xc(1,4,3)     = -xc.lb*sin(xc.delta)*(R_KE(1,1)*sin(xc.eta)-R_KE(2,1)*cos(xc.eta));     
CK_xc(2,4,3)     = -xc.lb*sin(xc.delta)*(R_KE(1,2)*sin(xc.eta)-R_KE(2,2)*cos(xc.eta));    
CK_xc(3,4,3)     = -xc.lb*sin(xc.delta)*(R_KE(1,3)*sin(xc.eta)-R_KE(2,3)*cos(xc.eta)); 
% Derivative with respect to eta
CK_xc(1,2,4)     = cos(xc.delta)*(R_KE(1,1)*sin(xc.eta)-R_KE(2,1)*cos(xc.eta));     
CK_xc(2,2,4)     = cos(xc.delta)*(R_KE(1,2)*sin(xc.eta)-R_KE(2,2)*cos(xc.eta));     
CK_xc(3,2,4)     = cos(xc.delta)*(R_KE(1,3)*sin(xc.eta)-R_KE(2,3)*cos(xc.eta));     

CK_xc(1,3,4)     = -xc.lb*sin(xc.delta)*(R_KE(1,1)*sin(xc.eta)-R_KE(2,1)*cos(xc.eta));     
CK_xc(2,3,4)     = -xc.lb*sin(xc.delta)*(R_KE(1,2)*sin(xc.eta)-R_KE(2,2)*cos(xc.eta));    
CK_xc(3,3,4)     = -xc.lb*sin(xc.delta)*(R_KE(1,3)*sin(xc.eta)-R_KE(2,3)*cos(xc.eta));    

CK_xc(1,4,4)     = xc.lb*cos(xc.delta)*(R_KE(1,1)*cos(xc.eta)+R_KE(2,1)*sin(xc.eta));     
CK_xc(2,4,4)     = xc.lb*cos(xc.delta)*(R_KE(1,2)*cos(xc.eta)+R_KE(2,2)*sin(xc.eta));    
CK_xc(3,4,4)     = xc.lb*cos(xc.delta)*(R_KE(1,3)*cos(xc.eta)+R_KE(2,3)*sin(xc.eta)); 
  
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Angular velocity   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative with respect to Theta
OmK_xs(1,2*NR+2,2*NR+1)   = -cos(xs.theta);
OmK_xs(2,2*NR+2,2*NR+1)   = -sin(xs.theta)*sin(xs.phi);
OmK_xs(3,2*NR+2,2*NR+1)   = -sin(xs.theta)*cos(xs.phi);
% Derivative with respect to Phi
OmK_xs(2,2*NR+1,2*NR+3)   = -sin(xs.phi);
OmK_xs(2,2*NR+2,2*NR+3)   = cos(xs.theta)*cos(xs.phi);

OmK_xs(3,2*NR+1,2*NR+3)   = -cos(xs.phi);
OmK_xs(3,2*NR+2,2*NR+3)   = -cos(xs.theta)*sin(xs.phi);

end

