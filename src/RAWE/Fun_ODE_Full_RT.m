function [DF DG r_N alpha_N r_t alpha_t m_0 kappa tau Ten F_NL RHS_r RHS_alpha] = Fun_ODE_Full_RT(tau0,xs)


E_m_0_tt = 0;
[DF DG r_N alpha_N r_t alpha_t m_0 kappa tau Ten F_NL RHS_r RHS_alpha] = Fun_ODE_RAWE(tau0,xs,E_m_0_tt);

end