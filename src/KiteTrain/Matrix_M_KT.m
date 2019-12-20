function M  = Matrix_M_KT(Ups,Phi,PND)
    

%-----------------------------------------------------------------------------
% Project   : LAKSA                                                          %
% Authors   : Gonzalo Sanchez-Arriaga and Jose A. Serrano-Iglesia            %
% Language  : Matlab                                                         %
% Synopsis  : Matrix M                                                       %
% Copyright:  Universidad Carlos III de Madrid, 2017. All rights reserved    %
%-----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: Ups, Phi                  - > Kinematic Matrices                   %% 
%        xs                        - > state and control vectors            %%
%        Dimensionless Parameters  - > PND                                  %%
% Output:  Matrix M                - > T = 0.5*xs_p'*M*xs_p                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
M = zeros(4*PND.Kite.N,4*PND.Kite.N);

for i=1:1:PND.Kite.N
   sigma   = PND.Kite.sigma(i);
   M       = M +  sigma*squeeze(Ups(:,:,i))'*squeeze(Ups(:,:,i)) + ... 
                  squeeze(Phi(:,:,i))'*squeeze(PND.Kite.iota(:,:,i))*squeeze(Phi(:,:,i));
end



end