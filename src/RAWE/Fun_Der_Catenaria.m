function DF = Fun_Der_Catenaria(s,xs)

global PP

mu    = PP(1);
sigma = PP(3);
nu    = PP(12);
vw    = PP(13);

x0      = xs(1,1);
z0      = xs(2,1);
x1      = xs(3,1);
z1      = xs(4,1);
x2      = xs(5,1);
z2      = xs(6,1);
x3      = xs(7,1);
z3      = xs(8,1);

nu0     = sqrt(x1^2 + z1^2);
kappa   = sqrt(x2^2 + z2^2);
Ten     = sigma*(nu0-1);


DF(1,1) = x1;
DF(2,1) = z1;
DF(3,1) = x2;
DF(4,1) = z2;
DF(5,1) = x3;
DF(6,1) = z3;
DF(7,1) = ( -nu*vw^2  + ( (sigma/nu0)*( x1*x2 + z1*z2) - 2*mu*( x2*x3 + z2*z3) )*x1 + (Ten - mu*kappa^2)*x2)/mu;
DF(8,1) = ( 1         + ( (sigma/nu0)*( x1*x2 + z1*z2) - 2*mu*( x2*x3 + z2*z3) )*z1 + (Ten - mu*kappa^2)*z2)/mu;




end