function c=cout_nu(param)

global spec_ref fid_p_lb NAA_inf NAA_sup time
global c

delta_nu=param(1);
% factor=param(2);

% correction=factor*exp(2*pi*1i*delta_nu*time);
correction=exp(2*pi*1i*delta_nu*time);

fid_shift_lb=fid_p_lb.*correction;

spec_shift_lb=fftshift(fft(fid_shift_lb));
spec_shift_lb=spec_shift_lb(NAA_inf:NAA_sup);

c=[abs(spec_shift_lb)-abs(spec_ref)];
