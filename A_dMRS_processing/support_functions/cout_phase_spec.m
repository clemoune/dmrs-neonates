function c=cout_phase_spec(param)

global spec_ref spec_p_lb factor_for_phi
global c

phi=param;

c=[real(factor_for_phi*exp(1i*phi)*spec_p_lb-spec_ref) imag(factor_for_phi*exp(1i*phi)*spec_p_lb-spec_ref)];
