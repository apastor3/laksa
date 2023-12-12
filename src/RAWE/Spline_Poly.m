function [Shape Shape_s Shape_ss Shape_sss] = Spline_Poly(chi,L_fe)


Shape      = zeros(4,length(chi));
Shape_s    = zeros(4,length(chi)); 
Shape_ss   = zeros(4,length(chi)); 
Shape_sss  = zeros(4,length(chi)); 

Shape(1,:) = -chi+1; 
Shape(2,:) = 1/6*(-chi.^3+3*chi.^2-2*chi); 
Shape(3,:) = chi; 
Shape(4,:) = 1/6*(chi.^3-chi); 

Shape_s(1,:) = -1; 
Shape_s(2,:) = 1/6*(-3*chi.^2+6*chi-2); 
Shape_s(3,:) =  1; 
Shape_s(4,:) = 1/6*(3*chi.^2-1); 

Shape_ss(1,:) =  0; 
Shape_ss(2,:) =  -chi+1; 
Shape_ss(3,:) =  0; 
Shape_ss(4,:) =  chi; 

Shape_sss(1,:) =  0; 
Shape_sss(2,:) = -1; 
Shape_sss(3,:) =  0; 
Shape_sss(4,:) =  1; 


Shape_s   = Shape_s/L_fe;
Shape_ss  = Shape_ss/L_fe^2;
Shape_sss = Shape_sss/L_fe^3;


end