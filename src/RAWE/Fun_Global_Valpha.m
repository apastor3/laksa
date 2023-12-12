function [V_alpha V_alpha_s V_alpha_ss V_alpha_tt] = Fun_Global_Valpha(NP,m_0,m_N,B_alpha_inv,Shape,Shape_s,Shape_ss,m_0_tt,m_N_tt)

% Recover Parameters
N_fe  = NP(1);
N_Int = NP(2);
   
% Find the vector with the BC
Vec_BC           = zeros(N_fe+1,1);
Vec_BC(1,1)      = - m_0/N_fe;
Vec_BC(end,1)    =   m_N/N_fe;

% Find the vector with the BC
Vec_BC_tt        = zeros(N_fe+1,1);
Vec_BC_tt(1,1)   = - m_0_tt/N_fe;
Vec_BC_tt(end,1) =   m_N_tt/N_fe;

% Multiply by the inverse (N_fe+1,1)
BB               = (B_alpha_inv*Vec_BC);
BB_tt            = (B_alpha_inv*Vec_BC_tt);
%%%% Compute Global Functions
V_alpha    = zeros(N_fe+1,N_fe*N_Int+1);
V_alpha_s  = zeros(N_fe+1,N_fe*N_Int+1);
V_alpha_ss = zeros(N_fe+1,N_fe*N_Int+1);
V_alpha_tt = zeros(N_fe+1,N_fe*N_Int+1);
%% Do i=1
i   = 1;
j   = i;
mm0 = 1;
for m = mm0:1:N_Int+1  % node     
    Ind  = (j-1)*N_Int+m;    
    V_alpha(i,Ind)    = BB(i)*Shape(2,m);
    V_alpha_s(i,Ind)  = BB(i)*Shape_s(2,m);
    V_alpha_ss(i,Ind) = BB(i)*Shape_ss(2,m);
    V_alpha_tt(i,Ind) = BB_tt(i)*Shape(2,m);
end

% Do i = 2 to N_fe
for i=2:1:N_fe % Function
    
    j= i-1;
    if j==1
       mm0 = 1;
    else
       mm0 = 2;
    end
    for m = mm0:1:N_Int+1  % node
        Ind  = (j-1)*N_Int+m;                                                                                                  
        V_alpha(i,Ind)      = BB(i)*Shape(4,m);
        V_alpha_s(i,Ind)    = BB(i)*Shape_s(4,m);
        V_alpha_ss(i,Ind)   = BB(i)*Shape_ss(4,m);
        V_alpha_tt(i,Ind)   = BB_tt(i)*Shape(4,m);
    end
    
    j = i;
    mm0 = 2;
    for m = mm0:1:N_Int+1  % node
       Ind  = (j-1)*N_Int+m;    
       V_alpha(i,Ind)    = BB(i)*Shape(2,m);
       V_alpha_s(i,Ind)  = BB(i)*Shape_s(2,m);
       V_alpha_ss(i,Ind) = BB(i)*Shape_ss(2,m);
       V_alpha_tt(i,Ind) = BB_tt(i)*Shape(2,m);
    end
end

% Do i = N_fe+1
i   = N_fe+1;
j   = i-1;
mm0 = 2;
for m = mm0:1:N_Int+1  % node
        Ind  = (j-1)*N_Int+m;                                                                                                  
        V_alpha(i,Ind)      = BB(i)*Shape(4,m);
        V_alpha_s(i,Ind)    = BB(i)*Shape_s(4,m);
        V_alpha_ss(i,Ind)   = BB(i)*Shape_ss(4,m);
        V_alpha_tt(i,Ind)   = BB_tt(i)*Shape(4,m);
end
    




end