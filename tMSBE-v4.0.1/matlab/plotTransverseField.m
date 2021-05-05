function plotTransverseField()


clear all
close all

global um;
global ps;

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

set(0,'defaulttextinterpreter','tex') %Default
%set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontName', 'AvantGarde')
set(0,'defaultAxesFontsize', 18)
set(0,'defaultTextFontsize', 18)
set(0,'defaultlinelinewidth',2) %Thin lines


outKey = '../run/refSpecQW__';

w0 = loadD([outKey,'w0.dat']);
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])


round_trip_time = loadD([outKey,'round_trip_time.dat']);
transverse_grid_y = loadD([outKey,'transverse_grid_y.dat']);

NUM_TRANSVERSE = length(dir([outKey,num2str(0),'_E_re_CAVP_T*.dat']))
plot_num = 0;

t = loadD([outKey,num2str(plot_num),'_t.dat']);
Eout = zeros(NUM_TRANSVERSE,length(t));
for i = 0:(NUM_TRANSVERSE-1)
    pulse_re = loadD([outKey,num2str(plot_num),'_E_re_CAVP_T',num2str(i),'.dat']);
    pulse_im = loadD([outKey,num2str(plot_num),'_E_im_CAVP_T',num2str(i),'.dat']);
    Eout(1+i,:) = (pulse_re + 1i*pulse_im).*exp(-1i*t*w0);
end

figure
hold on
surf(t/ps,transverse_grid_y/um,(0.5*eps0*c0*abs(Eout).^2)*cm*cm/1e6,'edgecolor','none')
hold off
xlabel('t [ps]')
ylabel('y [um]')
zlabel('I(t) [MW/cm^2]')
grid on
view(43,12)
%xlim([0,5])

asd

figure
hold on
surf(t/ps,transverse_grid_y/um,(0.5*eps0*c0*abs(Eout).^2)*cm*cm/1e6,'edgecolor','none')
hold off
xlabel('t [ps]')
ylabel('y [um]')
zlabel('I(t) [MW/cm^2]')
xlim([0,5])
xl = xlim;
yl = ylim;
hold on
for i = 1:length(transverse_grid_y)
   % plot3(xl,transverse_grid_y(i).*[1,1]/um,(0.5*eps0*c0*max(abs(Eout(i,:))).^2)*[1,1]*cm*cm/1e6,'k-');
end
hold off
grid on
view(90,90)








end