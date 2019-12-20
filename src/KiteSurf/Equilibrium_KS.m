function [u0  Error Flag]=Equilibrium_KS(u0,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Equilibrium state                                              %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%                                                                    %%
% Inputs:  u0       -> Initial Guess                                 %%
%          PND      -> Dimensionless parameters                      %%
%                                                                    %%
% Outputs: u0       -> Vector that makes | function(0,u0) <Error     %%
%          Error    -> Error of the solution                         %%
%          Flag     -> 1-> Success, 0 -> The method did not converge %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save original control and wind types
Control_Type = PND.Ctr.Type;
Wind_Type    = PND.Env.Type;
 
% Set steady wind and no control
PND.Ctr.Type = 0;
PND.Env.Type = 0;

if u0 == 0 % The user did not insert any initial guess 
    xs0              = [pi/4 -pi/6]'; 
end
[xu0 Error Flag] = my_fzero(@Fun_ODE_KS_EQ_SP,xs0,PND);
u0               = [0 xu0(1) 0 xu0(2) 0    0 0 0 0 0 ]';   

% Set original variables
PND.Ctr.Type = Control_Type;
PND.Env.Type = Wind_Type;

    function DF = Fun_ODE_KS_EQ_SP(t,X)
        w00     = [0; X(1); 0; X(2);0; zeros(5,1)];
        DF0     = Fun_ODE_KS(t,w00);
        DF(1,1) = DF0(7,1);
        DF(2,1) = DF0(9,1);
    end



end