function vw  = Fun_Wind(t,rg,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute Wind velocity                                          %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %%
% Inputs:  t                  -> dimensionless time                %%
%          rg                 -> dimensionless vector position of  %%
%                                the center of mass of the kite    %% 
%                                (components in SE)                %%
%          PND                -> dimensionless parameters          %%
%                                                                  %%
% Outputs: vw                 -> components in SE of the           %%
%          dimensionless  wind velocity vector                     %%
%                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch PND.Env.Type
    case 0 % Constant wind speed
      
        vw = -PND.Env.vw*[1 0 0]';
      
    case 1 % Wind Speed = Vw*(h/H0)^alfa*(1+eps*sin(Omega*t))
      
        H  = -rg(3);  % Normalized kite altitude
        vw = -PND.Env.vw*(H/PND.Env.H0)^PND.Env.alfa*(1+PND.Env.eps*sin(PND.Env.Omega*t))*[1 0 0]';
end


end