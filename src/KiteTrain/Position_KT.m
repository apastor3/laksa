function rk   = Position_KT(R1E,R21,RB2,ZX,PND)
%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : SB-S2 rotation matrix                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  xs                                                      %%
%          R1E   -> S1-SE rotation matrix                          %% 
%          R21   -> S2-S1 rotation matrix                          %%
%          RB2   -> SB-S2 rotation matrix                          %%
%          ZX    -> [zeta xi]                                      %%
%          PND   -> Dimensionless parameters                       %%
% Outputs: rk    -> SE-components of kite's position vectors       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
rk         = zeros(3,PND.Kite.N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Compute center of mass positions  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Kite
zeta     = ZX(1,1);
xi       = ZX(2,1);

R2E      = R21(:,:,1)*R1E(:,:,1);
RBE      = RB2(:,:,1)*R2E;
rk(:,1)  = -R2E'*[0 zeta xi]'-RBE'*[PND.Tether.XA(1) 0 PND.Tether.ZA(1) ]';
% Other Kites
for i=2:1:PND.Kite.N
    zeta     = ZX(1,i);
    xi       = ZX(2,i);
    
    R2E     = squeeze(R21(:,:,i))*squeeze(R1E(:,:,i));
    RBE     = squeeze(RB2(:,:,i))*R2E;
    rk(:,i) = rk(:,i-1) -R2E'*[0 zeta xi]'-RBE'*[PND.Tether.XA(i) 0 PND.Tether.ZA(i) ]'; 
end




end