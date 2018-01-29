function [Ypo,Error,Floquet,flag] = Corrector_Non_Autonomous(Fun_DF,Y0,T0,Itmax,Tol,TolRel,TolAbs,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Periodic Orbit Corrector for Non-autonomous system             %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
% Inputs                                                                     %
% Fun_DF                    - > Right-hand side                              % 
% Y0                        - > Initial Condition guess                      %
% T0                        - > Period of the periodic orbit                 %
% Itmax, Tol, TolRel,TolAbs - > Numerical Parameters                         %
% PND                       - > Dimensionless Parameters                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                                                                     %
% Yp0                      - > Initial condition                             % 
% Error                    - > Periodic Orbit Error                          %
% Floquet                  - > Floquet Multipliers                           %
% Flag                     - > Exit Flag                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize some variables and options
N        = length(Y0);
Iden     = eye(N);
%options  = odeset('RelTol',TolRel,'AbsTol',ones(1,N*(N+1))*TolAbs);
%options2 = odeset('RelTol',TolRel,'AbsTol',ones(1,N)*TolAbs);
Floquet  = zeros(1,N);
M        = zeros(N,N);
    
for cont1 = 1:1:Itmax
        Yvar0 = [];
      %  [~, Y]      = ode45(Fun_DF,[0 T0],Y0,options2); 
        [~, Y]      = ode45(Fun_DF,[0 T0],Y0); 
     
        Error       = max(abs(Y(1,1:N) - Y(end,1:N)));
        display(['Corrector: iteration no. ' num2str(cont1) ',  Error = ' num2str(Error)])
        Yvar0(1,1:N) = Y0';
        
        for cont2 = 1:1:N
            Yvar0(1,N + (cont2-1)*N + 1:N + cont2*N) = Iden(cont2,:);
        end
        Yvar0 = Yvar0';

        if Error < Tol      

            if cont1 == 1       
               %[~, Yvar]    = ode45(@FUNvariationalDF,[0 T0],Yvar0,options);
               [~, Yvar]    = ode45(@FUNvariationalDF,[0 T0],Yvar0);
               
               for cont2 = 1:1:N
                   for cont3 = 1:1:N
                       M(cont2,cont3)   = Yvar(end,N + (cont2-1)*N + cont3);
                   end
                end
                M = M';
                [~, Val] = eig(M);
                for cont2 = 1:1:N
                    Floquet(1,cont2)  = Val(cont2,cont2);
                end
            end
            Ypo = Y0;
            flag = 1;
            break
        end % If Err<Tol

     %   [~, Yvar]    = ode45(@FUNvariationalDF,[0 T0],Yvar0,options);
        [~, Yvar]    = ode45(@FUNvariationalDF,[0 T0],Yvar0);
    
        YF        = Yvar(end,1:N*(N+1))';
        for cont2 = 1:1:N
            for cont3 = 1:1:N
                M(cont2,cont3) = Yvar(end,N+(cont2-1)*N + cont3);
            end
        end
        M         = M';
        A         = M - eye(N);
        B         = -(YF(1:N) - Y0);
        Correc    = A\B;
        for cont2 = 1:1:N
            Y0(cont2,1)     = Y0(cont2,1) + Correc(cont2);
        end
        [~, Val] = eig(M);
        for cont2 = 1:1:N
            Floquet(1,cont2)   = Val(cont2,cont2);
        end
        if cont1 == Itmax
            Ypo  = Y0;
            flag = -1;    
            return
        end
end

    function  DF = FUNvariationalDF(t,YV)
        
        Jac       = Jacobian(Fun_DF,t,YV(1:N,1),PND);      
        DF(1:N,1) = feval(Fun_DF,t,YV(1:N,1));   

          
        for cont1 = 1:1:N                                                            % Loop for Xi variational
            for cont2 = 1:1:N                                                        % Index of the Xi
               DF(N + (cont1-1)*N + cont2) = 0;
               for cont3 = 1:1:N                                                    % Loop for the sum
                 DF(N + (cont1-1)*N + cont2) = DF(N + (cont1-1)*N + cont2) + Jac(cont2,cont3)*YV(N + N*(cont1-1) + cont3);
               end
            end
        end
    end

    

end  %% End Corrector
