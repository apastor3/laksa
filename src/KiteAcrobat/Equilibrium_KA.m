function [u0  Error Flag]=Equilibrium_KA(u0,PND)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute equilibrium state                                      %
% Copyright :  Universidad Carlos III de Madrid, 2017. All rights reserved   %
%----------------------------------------------------------------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %%
% Inputs:  u0       -> Initial Guess                                 %%
%          PND      -> Dimensionless parameters                      %%
%                                                                    %%
% Outputs: u0       -> Vector that makes | function(0,X) <Error      %%
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
    X  = [75 -65 ]'*pi/180;
else
    X(1,1) = u0(1,1);
    X(2,1) = u0(3,1);
end

[X Error Flag]=my_fzero(@Equilibrium_KA,X,PND);

u0 = [0; X(1); 0; X(2);zeros(4,1)];

 
%% Check Equilibrium
if Flag == 0
    display('Newton method did not converge')
    display(['Error = ' num2str(Error)])
else
    [Tout RBE rk  vk  ak euler omega omega_p Lambda FAP FAM MAP MAM FA MA W alfa  beta LP LM] = Fun_Post_KA(0,PND,0,u0,0); 
                                                                                            
    Ok = 1;
    if Lambda(1)<0
         display('Tether tensions are negative')
         Ok = 0;
    end
    if alfa > PND.Aero.alfa_s
         display('Attack Angle is above the stall ')
         Ok = 0;
    end
    if alfa <0
         display('Attack Angle is negative ')
         Ok = 0
    end 
end

if Flag==0 || Ok==0
     display('Press intro to continue')
     pause
end

% Set original control and wind types
PND.Ctr.Type = Control_Type;
PND.Env.Type = Wind_Type;


    function DF = Equilibrium_KA(t,X)
              
        w00     = [0; X(1); 0; X(2); zeros(4,1)];
        
        DF0     = Fun_ODE_Ham_KA(t,w00);
        
        DF(1,1) = DF0(6,1);
        DF(2,1) = DF0(8,1);
    end

end