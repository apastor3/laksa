function DF = Fun_ODE_RT_NC(tau,xs)

xs_amp = [xs;0];

[DF_amp DG r_N alpha_N r_t alpha_t m0 kappa tau Ten F_NL RHS_r RHS_alpha] = Fun_ODE_Full_RT(tau,xs_amp);

DF = DF_amp(1:end-1,1);

end