function n = indexAlXGaZAs(lambda,x)
%% Refractive index of Al(x)Ga(1-x)As
% From Kokubo and Ohta 1997
% lambda is the wavelength in units of meters

h = 6.62606957*1e-34; % [Js]
c0 = 299792458; %[m/s]
e = 1.60217657*1e-19;

E = h*c0./lambda;
E = E/e;

n = 3.3 + 0.09*x - (0.08 + 0.7*x)*E + (0.19 + 0.16*x)*E.^2 + 0.00023./((E - (1.42 + 1.25*x)).^2 + 0.003);

end