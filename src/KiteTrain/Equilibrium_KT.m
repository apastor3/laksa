function [u0  Error Flag]=Equilibrium_KT(u0,PND)


%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : Matrix M                                                       %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
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


if u0 == 0 % The user did not insert any initial guess
    % Do it for one kite
    NK             = PND.Kite.N;
    PND.Kite.N     = 1;
    gamma          = 45*pi/180;
    alfa           = 10*pi/180;
    X              = [gamma alfa-gamma ]';
    [X Error Flag] = my_fzero(@Equilibrium_KA,X,PND);
    PND.Kite.N     = NK;
    for i=1:1:PND.Kite.N
        X(2*(i-1)+1:2*i,1)  = [X(1) X(2) ]';
    end
    
   % X = [0.398660178875186 -0.263975143869782 0.398660178875186 -0.263975143869782]';
   % DF = Equilibrium_KA(0,X)
   % pause
else
    for i=1:1:PND.Kite.N
      X(2*(i-1)+1:2*i,1) = [u0(4*(i-1)+1,1) u0(4*(i-1)+4,1)]';
    end
end

[X Error Flag]=my_fzero(@Equilibrium_KA,X,PND);



for i=1:1:PND.Kite.N
  u0(4*(i-1)+1:4*i,1) = [0 X(2*i-1) 0 X(2*i)]';
end
u0 = [u0;zeros(4*PND.Kite.N ,1)];
 
%% Check Equilibrium
%if Flag == 0
%    display('Newton method did not converge')
%    display(['Error = ' num2str(Error)])
%else
%    [Tout RBE rk  vk  ak euler omega omega_p Lambda FAP FAM MAP MAM FA MA W alfa  beta LP LM] = Fun_Post_KA(0,PND,0,u0,0); 
%                                                                                            
%    Ok = 1;
%    if Lambda(1)<0
%         display('Tether tensions are negative')
%         Ok = 0;
%    end
%    if alfa > PND.Aero.alfa_s
%         display('Attack Angle is above the stall ')
%         Ok = 0;
%    end
%    if alfa <0
%         display('Attack Angle is negative ')
%         Ok = 0
%    end 
%end

%if Flag==0 || Ok==0
%     display('Press intro to continue')
%     pause
%end

    function DF = Equilibrium_KA(t,X)
        for i=1:1:PND.Kite.N      
            w00(4*(i-1)+1:4*i,1) = [0 X(2*i-1) 0 X(2*i)]';
        end
        w00 = [w00;zeros(4*PND.Kite.N ,1)];
        
        [xc RBE rK vK omegaK FA_K MA_K alfa beta Ups Ups_xs Phi Phi_xs Q DF0 RHS]=  Fun_ODE_Full_Output_KT(0,w00,PND);
        
        for i=1:1:PND.Kite.N  
            DF(2*i-1,1) = RHS( 4*(i-1)+2,1);
            DF(2*i,1)   = RHS( 4*(i-1)+4,1);
        end
    
    end

end