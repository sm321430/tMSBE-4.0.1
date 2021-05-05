global um;
global ps;
global hbar;
global e;
global c0;

fs = 1.0e-15;
ps = 1.0e-12;
um = 1.0e-6;
ns = 1.0e-9;
cm = 1.0e-2;
nm = 1.0e-9;

hbar = 1.054589e-34;
e = 1.602189e-19;
c0   = 2.99792458E+08;
mu0  = (4.0e-7)*pi;
eps0 = 1.0/(mu0*c0*c0);
n=3.4453; %Background refractive index (change as needed)

m0   = 9.109389754e-31;
me   = 0.06085*m0;
mh   = 0.234*m0;
a0   = 1.062146e-08;
mr = me*mh/(me+mh);
Eb = (hbar^2 / (2*mr*a0*a0));