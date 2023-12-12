function [s_fe s_nc phi phi_s phi_ss  phi_sss varphi varphi_s varphi_ss Shape Shape_s Shape_ss Shape_sss B_alpha_inv] = Fun_Global(N_fe,N_Int)
        
% Recover Parameters
L_fe         = 1/N_fe;                     % Length of finite element

% Compute Grids
s_fe    = [0:1:N_fe]'/N_fe;                % Grid of finite elements
s_nc    = [0:1:N_Int*N_fe]'/(N_Int*N_fe);  % Grid of Newton-Cotes integration

% Compute Splines
chi                                = [0:1:N_Int]/N_Int;              % Newton-Cotes local grid values
[Shape Shape_s Shape_ss Shape_sss] = Spline_Poly(chi,L_fe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     Compute Global Function for r                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute matrix to write x'' = A_M1*x

B_r      = zeros(N_fe-1,N_fe-1);
C_r      = zeros(N_fe-1,N_fe+1);

for j=1:1:N_fe-1
    if j>1
       B_r(j,j-1) = 1/6;
    end
    B_r(j,j)   = 2/3;
    if j<N_fe-1
       B_r(j,j+1) = 1/6;
    end
    C_r(j,j)      =  1;
    C_r(j,j+1)    = -2;
    C_r(j,j+2)    =  1;
end

A_M1 = inv(B_r)*C_r;
A_M1 = [zeros(1,N_fe+1);A_M1;zeros(1,N_fe+1)];

%%%% Compute Global Functions

phi     = zeros(N_fe+1,N_fe*N_Int+1);
phi_s   = zeros(N_fe+1,N_fe*N_Int+1);
phi_ss  = zeros(N_fe+1,N_fe*N_Int+1);
phi_sss = zeros(N_fe+1,N_fe*N_Int+1);



%% Function 2-N_fe

for i=1:1:N_fe+1 % Function
    for j=1:1:N_fe % Element
        if j==1
            mm0 = 1;
        else
            mm0 = 2;
        end
        for m = mm0:1:N_Int+1  % node
            Ind  = (j-1)*N_Int+m;
  
            phi(i,Ind)         = phi(i,Ind)     +  A_M1(j,i)*Shape(2,m)     + A_M1(j+1,i)*Shape(4,m);                                                 
            phi_s(i,Ind)       = phi_s(i,Ind)   +  A_M1(j,i)*Shape_s(2,m)   + A_M1(j+1,i)*Shape_s(4,m);
            phi_ss(i,Ind)      = phi_ss(i,Ind)  +  A_M1(j,i)*Shape_ss(2,m)  + A_M1(j+1,i)*Shape_ss(4,m);
            phi_sss(i,Ind)     = phi_sss(i,Ind) +  A_M1(j,i)*Shape_sss(2,m) + A_M1(j+1,i)*Shape_sss(4,m);
                                                                                                            
            if i-1==j
                phi(i,Ind)     = phi(i,Ind)     + Shape(3,m);
                phi_s(i,Ind)   = phi_s(i,Ind)   + Shape_s(3,m);
                phi_ss(i,Ind)  = phi_ss(i,Ind)  + Shape_ss(3,m);
                phi_sss(i,Ind) = phi_sss(i,Ind) + Shape_sss(3,m);
            end
            if i==j
                phi(i,Ind)     = phi(i,Ind)     + Shape(1,m);
                phi_s(i,Ind)   = phi_s(i,Ind)   + Shape_s(1,m);
                phi_ss(i,Ind)  = phi_ss(i,Ind)  + Shape_ss(1,m);
                phi_sss(i,Ind) = phi_sss(i,Ind) + Shape_sss(1,m);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     Compute Global Function for alpha               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute matrix to write x'' = A_M1*x

B_alpha      = zeros(N_fe+1,N_fe+1);
C_alpha      = zeros(N_fe+1,N_fe+1);

% First Row
B_alpha(1,1) = 1/3;
B_alpha(1,2) = 1/6;

C_alpha(1,1) = -1;
C_alpha(1,2) =  1;

for j=2:1:N_fe
    B_alpha(j,j-1) = 1/6;
    B_alpha(j,j)   = 2/3;
    B_alpha(j,j+1) = 1/6;
    
    C_alpha(j,j-1) =  1;
    C_alpha(j,j)   = -2;
    C_alpha(j,j+1) =  1;
end

% Last Row
B_alpha(N_fe+1,N_fe)   = 1/6;
B_alpha(N_fe+1,N_fe+1) = 1/3;

C_alpha(N_fe+1,N_fe)   =  1;
C_alpha(N_fe+1,N_fe+1) = -1;

% Do inverse
B_alpha_inv = inv(B_alpha);
% Dot product
BC_alpha    = B_alpha_inv*C_alpha;


%%%% Compute Global Functions

varphi     = zeros(N_fe+1,N_fe*N_Int+1);
varphi_s   = zeros(N_fe+1,N_fe*N_Int+1);
varphi_ss  = zeros(N_fe+1,N_fe*N_Int+1);
 
%% Function 2-N_fe
for i=1:1:N_fe+1 % Function
    for j=1:1:N_fe % Element
        if j==1
            mm0 = 1;
        else
            mm0 = 2;
        end
        for m = mm0:1:N_Int+1  % node
            Ind  = (j-1)*N_Int+m;
  
            varphi(i,Ind)         = varphi(i,Ind)     +  BC_alpha(j,i)*Shape(2,m)     + BC_alpha(j+1,i)*Shape(4,m);                                                 
            varphi_s(i,Ind)       = varphi_s(i,Ind)   +  BC_alpha(j,i)*Shape_s(2,m)   + BC_alpha(j+1,i)*Shape_s(4,m);
            varphi_ss(i,Ind)      = varphi_ss(i,Ind)  +  BC_alpha(j,i)*Shape_ss(2,m)  + BC_alpha(j+1,i)*Shape_ss(4,m);
                                                                                                            
            if i-1==j
                varphi(i,Ind)     = varphi(i,Ind)     + Shape(3,m);
                varphi_s(i,Ind)   = varphi_s(i,Ind)   + Shape_s(3,m);
                varphi_ss(i,Ind)  = varphi_ss(i,Ind)  + Shape_ss(3,m);
            end
            if i==j
                varphi(i,Ind)     = varphi(i,Ind)     + Shape(1,m);
                varphi_s(i,Ind)   = varphi_s(i,Ind)   + Shape_s(1,m);
                varphi_ss(i,Ind)  = varphi_ss(i,Ind)  + Shape_ss(1,m);
            end
        end
    end
end





end