function [ZX R2 Aux_S1 Aux_S2 Aux_Sq] = Circle(x1,y1,r1,x2,y2,r2)

R2     = (x2-x1)^2+(y2-y1)^2;
Aux_S1  = (r1+r2)^2-R2;
Aux_S2 =  R2-(r1-r2)^2;
Aux_Sq = sqrt(Aux_S1*Aux_S2);

x      = 0.5*(x1+x2)+(r1^2-r2^2)/(2*R2)*(x2-x1)+0.5*Aux_Sq*(y2-y1)/R2; 
y      = 0.5*(y1+y2)+(r1^2-r2^2)/(2*R2)*(y2-y1)+0.5*Aux_Sq*(x1-x2)/R2;

ZX     = [x y]';


E(1) = r1^2-(ZX(1)-x1)^2-(ZX(2)-y1)^2;
E(2) = r2^2-(ZX(1)-x2)^2-(ZX(2)-y2)^2;
if max(abs(E))>1e-7 || abs(R2)<1e-10
    display(['Error in function Circle = ' num2str(max(abs(E)))]) 
    R2
end

end