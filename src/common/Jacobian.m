function Jnum=Jacobian(funcion,t,X,PND)

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga, Alejandro Pastor-Rodriguez,           %
% Language  : Matlab                                                         %
% Synopsis  : Compute Jacobian                                               %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %%
% Inputs:  fun   -> function(t,x)                                   %%
%          t     -> time                                            %%
%          X     -> State vector                                    %%
%          PND   -> Dimensionless parameters                        %%
%                                                                   %%
% Outputs: Jnum  -> Jacobian of the function                        %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:length(X)
    
    XM    = X;
    XM(i) = XM(i)+PND.Num.dh/2;
    fM    = feval(funcion,t,XM);
    
    Xm    = X;
    Xm(i) = Xm(i)-PND.Num.dh/2;
    fm    = feval(funcion,t,Xm);
       
    Jnum(:,i)=(fM-fm)/PND.Num.dh;
end