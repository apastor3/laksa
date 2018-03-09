function [Ms Msc Mc] = Matrix_M_KF(SK,CK,OmK,SR,CR,OmR,SG,CG,xs,xc,PND)
           
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : MAtrix M to compute the kinetic energy                         %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: SK,CK,OmK,SR,CR,OmR,SG,CG - > Kinematic Matrices                   %% 
%        xs, xc                    - > state and control vectors            %%
%        Dimensionless Parameters  - > PND                                  %%
% Output:  Matrices Ms Msc and Mc  - > T = 0.5*( xs_p'*Ms*xs_p  +           %%
%                                               2xs_p'*Msc*xc_p +           %%
%                                                xc_p'*Mc*xc_p )            %%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
NR   = PND.Num.N ;   % Number of Rods
NG   = PND.Gen.Num;  % Number of Generators
Ns   = 2*NR+3;       % Rods and Euler
Ns_p = 2*NR+3+NG;    % Rods, Eulers, and generators

% Kite Parameters
iK        = PND.Kite.ik;
% Tether Parameter
sigmaR    = PND.Tether.Sigma;
iR        = PND.Tether.Ups;

if NG>0
    % Generator Parameters
    iG        = PND.Gen.iota; % Generator tensor of inertia (divided by Sigma*lG^2)
    sigmaG    = PND.Gen.Sigma;
    lG        = PND.Gen.l;
    nu        = PND.Gen.nu;
end
% Kite Contribution
Ms  = SK'*SK+OmK'*iK*OmK;
Msc = SK'*CK;
Mc  = CK'*CK;

% Tether Contribution
for j=1:1:NR
  Aux1 = (squeeze(SR(:,:,j)))'*squeeze(SR(:,:,j));
  Aux2 = (squeeze(OmR(:,:,j)))'*iR*squeeze(OmR(:,:,j));
  Ms   = Ms + sigmaR*xc.lt*(Aux1+xc.lt^2*Aux2);
  
  Aux1 = (squeeze(SR(:,:,j)))'*(squeeze(CR(:,:,j)));
  Msc  = Msc + sigmaR*xc.lt*Aux1;
  
  Aux1 = (squeeze(CR(:,:,j)))'*(squeeze(CR(:,:,j)));
  Mc   = Mc + sigmaR*xc.lt*Aux1;  
end


% Generator Contribution
theta = xs.theta;
psi   = xs.psi;
phi   = xs.phi;


for j=1:1:NG
    
  Aux1 = (squeeze(SG(:,:,j)))'*squeeze(SG(:,:,j));
 
  % term omega_G * iG * omega_G 
  Aux2 = zeros(Ns_p,Ns_p);
  % Terms lambda_dot^2
  Aux2(2*NR+3+j,2*NR+3+j)      =  iG(1,1);
  % Terms lambda_dot by (theta_dot, psi_dot, and phi_dot )
  Aux2(2*NR+3+j,2*NR+1:2*NR+3) =  iG(1,1)*[sin(nu(j))*sin(phi) -cos(nu(j))*sin(theta)-sin(nu(j))*cos(theta)*cos(phi) cos(nu(j))];
  Aux2(2*NR+1:2*NR+3,2*NR+3+j) =  iG(1,1)*[sin(nu(j))*sin(phi) -cos(nu(j))*sin(theta)-sin(nu(j))*cos(theta)*cos(phi) cos(nu(j))]';
 
  % Terms (theta_dot, psi_dot, and phi_dot )*(theta_dot, psi_dot, and phi_dot )
  B(1,1:3)  = [ iG(1,1)*(cos(nu(j)))^2+iG(2,2)*(sin(nu(j)))^2                  0      -iG(1,1)*sin(nu(j))*cos(nu(j))+iG(2,2)*sin(nu(j))*cos(nu(j)) ];
  B(2,1:3)  = [0                                                            iG(2,2)               0 ];
  B(3,1:3)  = [-iG(1,1)*sin(nu(j))*cos(nu(j))+iG(2,2)*sin(nu(j))*cos(nu(j))    0       iG(1,1)*(sin(nu(j)))^2+iG(2,2)*(cos(nu(j)))^2];
  Aux2(2*NR+1:2*NR+3,2*NR+1:2*NR+3) = OmK(:,2*NR+1:2*NR+3)'*B*OmK(:,2*NR+1:2*NR+3);
    
  Ms   = Ms + sigmaG*(Aux1+lG^2*Aux2);
  
  Aux1 = (squeeze(SG(:,:,j)))'*(squeeze(CG(:,:,j)));
  Msc  = Msc + sigmaG*Aux1;
  
  Aux1 = (squeeze(CG(:,:,j)))'*(squeeze(CG(:,:,j)));
  Mc   = Mc + sigmaG*Aux1;  
end

end