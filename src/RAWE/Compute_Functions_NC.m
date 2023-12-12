function [r_NC r_NC_s r_NC_ss r_NC_sss alpha_NC alpha_NC_s alpha_NC_ss ] = Compute_Functions_NC(NP,r_N,alpha_N,phi,phi_s,phi_ss,phi_sss,varphi,varphi_s,varphi_ss,V_alpha,V_alpha_s,V_alpha_ss); 

N_nc           = NP(4);

r_NC        = zeros(N_nc,3);  
r_NC_s      = zeros(N_nc,3);  
r_NC_ss     = zeros(N_nc,3);
r_NC_sss    = zeros(N_nc,3);
alpha_NC    = zeros(N_nc,1);
alpha_NC_s  = zeros(N_nc,1);
alpha_NC_ss = zeros(N_nc,1);


for i=1:1:NP(1)+1
         % Function
         r_NC(:,1) = r_NC(:,1) + r_N(i,1)*phi(i,:)';
         r_NC(:,2) = r_NC(:,2) + r_N(i,2)*phi(i,:)';
         r_NC(:,3) = r_NC(:,3) + r_N(i,3)*phi(i,:)';
       
         % First Derivative of the Function
         r_NC_s(:,1) = r_NC_s(:,1) + r_N(i,1)*phi_s(i,:)';
         r_NC_s(:,2) = r_NC_s(:,2) + r_N(i,2)*phi_s(i,:)';
         r_NC_s(:,3) = r_NC_s(:,3) + r_N(i,3)*phi_s(i,:)';
    
         % Second Derivative of the Function
         r_NC_ss(:,1) = r_NC_ss(:,1) + r_N(i,1)*phi_ss(i,:)';
         r_NC_ss(:,2) = r_NC_ss(:,2) + r_N(i,2)*phi_ss(i,:)';
         r_NC_ss(:,3) = r_NC_ss(:,3) + r_N(i,3)*phi_ss(i,:)';        
   
         % Second Derivative of the Function
         r_NC_sss(:,1) = r_NC_sss(:,1) + r_N(i,1)*phi_sss(i,:)';
         r_NC_sss(:,2) = r_NC_sss(:,2) + r_N(i,2)*phi_sss(i,:)';
         r_NC_sss(:,3) = r_NC_sss(:,3) + r_N(i,3)*phi_sss(i,:)';        
           
         % Alpha
         alpha_NC    = alpha_NC   + alpha_N(i)*varphi(i,:)'    + V_alpha(i,:)';
         
         % Alpha_s
         alpha_NC_s  = alpha_NC_s + alpha_N(i)*varphi_s(i,:)'  + V_alpha_s(i,:)';
         
         % Alpha_ss
         alpha_NC_ss = alpha_NC_ss + alpha_N(i)*varphi_ss(i,:)'+ V_alpha_ss(i,:)';
end



end