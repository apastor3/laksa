function    [Em Ec Ep H] = Compute_Energy_KF(SR,CR,OmR,SK,CK,OmK,SG,CG,xs_amp,xc_amp,PND)
 
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute the energy of the system                               %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %%
%     SR,CR,OmR,SK,CK,OmK,SG,CG - > Kinematic Matrices                    %%
%     xs_amp,xc_amp             - > State and control vectors             %%
%     PND                       - > dimensionless Parameters              %%
% Output                                                                  %%
%     Em = Ec+Ep -> Mechanical energy                                     %% 
%     Ec         -> Kinetic Energy                                        %%
%     Ep         -> Potential Energy                                      %%
%     H          -> Hamiltonian                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recover variables
[xs xs_p] = From_xs2Var_KF(xs_amp,PND);
[xc xc_p] = From_xc2Var_KF(xc_amp,PND);

NR = PND.Num.N ;   % Number of Rods
NG = PND.Gen.Num;  % Number of Generators

%
Nvar      = 2*NR+3;
Nvar_p    = 2*NR+3+NG;
Nc        = 4+NG+3;
Nc0       = 4;

% Kite Parameters
iK        = PND.Kite.ik;
% Tether Parameter
sigmaR    = PND.Tether.Sigma;
iR        = PND.Tether.Ups;
% Generator Parameters
if NG>0
    iG        = PND.Gen.iota; % Generator tensor of inertia (divided by Sigma*lG^2)
    sigmaG    = PND.Gen.Sigma;
    lG        = PND.Gen.l;
else
    sigmaG    = 0;
end

[Ms Msc Mc] = Matrix_M_KF(SK,CK,OmK,SR,CR,OmR,SG,CG,xs,xc,PND);

us        = xs_amp(1:Nvar,1);
us_p      = xs_amp(Nvar+1:Nvar+Nvar_p,1);
uc        = xc_amp(1:Nc,1);
uc_p      = xc_amp(Nc+1:Nc+4,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Kinetic energy   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
Ec  = (1/2)*(us_p'*Ms*us_p + 2*us_p'*Msc*uc_p + uc_p'*Mc*uc_p );

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Potential energy  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kite and generators contribution from QG
Ep  = (1+NG*sigmaG)*xc.lb*(cos(xs.theta)*(cos(xc.delta)*sin(xc.eta)*sin(xs.phi)+sin(xc.delta)*cos(xs.phi))-sin(xs.theta)*cos(xc.delta)*cos(xc.eta));
% Tether and Kite and generators contributions from rQ
for j=1:1:NR
    
    Ep = Ep + (1+NG*sigmaG)*xc.lt*sin(xs.gamma(j));
    
    for k=1:1:j
        if k==j
            g = 1/2;
        else
            g = 1;
        end
        Ep = Ep + sigmaR*xc.lt^2*g*sin(xs.gamma(k));
    end
end
% Generators
for i=1:1:NG
     Ep = Ep + sigmaG*(PND.Gen.x(i)*sin(xs.theta)-PND.Gen.y(i)*cos(xs.theta)*sin(xs.phi)-PND.Gen.z(i)*cos(xs.theta)*cos(xs.phi));
end


% Uncomment this to check the potential energy
%[R_KE Grad_R_KE R_GK Grad_R_GK] = Matrix_R(xs.theta,xs.psi,xs.phi,xs.lambda,PND);
%[rQ rR0 rK0 rG0 ]               = Compute_Positions(xs,xc,R_KE,PND);
    
%Ep0 = -rK0(3);
%for i=1:1:NR
%    Ep0 = Ep0 - sigmaR*xc.lt*rR0(3,i);
%end
%for i=1:1:NG
%    Ep0 = Ep0 - sigmaG*rG0(3,i);
%end
%Ep-Ep0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Total Mechanical Energy %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Em = Ec+Ep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Hamiltonian    %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H  = 1/2*(us_p'*Ms*us_p - uc_p'*Mc*uc_p)+Ep; 

end
