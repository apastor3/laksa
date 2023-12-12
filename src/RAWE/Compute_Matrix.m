function [P0_Inv Q0_Inv p0 P4 Q1 P00 Q0 P_iN] = Compute_Matrix(N_fe,N_Int)


[s_fe s_nc phi phi_s phi_ss  phi_sss varphi varphi_s varphi_ss Shape Shape_s Shape_ss Shape_sss B_alpha_inv] = Fun_Global(N_fe,N_Int);

P0  = zeros(N_fe+1,N_fe+1);
P4  = zeros(N_fe+1,N_fe+1);
p0  = zeros(N_fe+1,1);

Q0  = zeros(N_fe+1,N_fe+1);
Q1  = zeros(N_fe+1,N_fe+1);

for i=1:1:N_fe+1
    for j=1:1:N_fe+1
       
        P0(i,j)  = trapz(s_nc,phi(j,:).*phi(i,:));
        P4(i,j)  = trapz(s_nc,phi_ss(j,:).*phi_ss(i,:));
        
        Q0(i,j)  = trapz(s_nc,varphi(j,:).*varphi(i,:));
        Q1(i,j)  = trapz(s_nc,varphi_s(j,:).*varphi_s(i,:));
    end
    
    p0(i,1)      = trapz(s_nc,phi(i,:));        
end
P00    = P0(2:N_fe,2:N_fe);
P0_Inv = inv(P00);
Q0_Inv = inv(Q0);

P_iN   = squeeze(P0(2:N_fe,N_fe+1));

end