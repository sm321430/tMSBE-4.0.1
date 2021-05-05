function plotGainSpectrums()

close all
clear all

fs = 1.0e-15;
ps = 1.0e-12;
um = 1.0e-6;
nm = 1.0e-9;

global hbar;
global e;
global c0;

hbar = 1.054589e-34;
e = 1.602189e-19;
c0   = 2.99792458E+08;
eps0 = 8.854187817620E-12;

m0   = 9.109389754e-31;
me   = 0.06085*m0;
mh   = 0.234*m0;
a0   = 1.062146e-08;
mr = me*mh/(me+mh);
Eb = (hbar^2 / (2*mr*a0*a0));


set(0,'defaulttextinterpreter','tex') %Default
%set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontsize', 50)
set(0,'defaultAxesLinewidth', 5)
set(0,'defaultTextFontsize', 50)
set(0,'defaultlineMarkerSize',26)
set(0,'defaultlinelinewidth',5) %Thin lines
maps %Load maps file for color scheme
POS=[1,1,1200,800];
POS2=[1,1,500,900];
%POS_plot=[0.14,0.1100,0.53125,0.85];
%POS_plot=[0.1330,0.1100,0.5930,0.8150];
my_lineStyle={'-kd','-bo','-r+', '-gx', '-c*','-mo','-k+','-kx','-k*','-b+','-bx','-b*','-g+'};

w0 = loadD('refSpecQW__w0.dat');
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])
disp(['lambda   = ',num2str(2*pi*c0/w0/nm,'%.2f'),' [nm]'])

% Set Frequency limits
W_LIM = w0+e*0.1/hbar*[-1.0,1.0];
W_NUM = 600;

%Load plot limits (GD, GDD, GDDD)
n_s = loadStruct('refSpecQW__system_structure.dat');
DBRind = n_s(1:2);
DBRind = [2.9459, 3.4345];
% Approx width of DBR
DF = w0*(4/pi)*asin((min(DBRind)-max(DBRind))/(sum(DBRind)));
fa = w0 + DF/2;
fb = w0 - DF/2;
DBR_LIM = [fa,fb];
%DBR_LIM = W_LIM;

%% Plot all electric fields

%% Reflection spectrum from QW
outKey = 'refSpecQW__';
n = 0;
t = loadD([outKey,num2str(n),'_t.dat']);
e_cav_re = loadD([outKey,num2str(n),'_E_re_CAVP_T0.dat']);
e_cav_im = loadD([outKey,num2str(n),'_E_im_CAVP_T0.dat']);

E_cav = (e_cav_re + 1i*e_cav_im).*exp(-1i*w0*t);
E_cav_abs = abs(E_cav);

%{
figure
plot(t/ps,E_cav_abs*um,'k-')
xlabel('t [ps]')
ylabel('E [V/um]')
grid on
legend('CAV')
title('QW: Input fields')
xlim([min(t)/ps, max(t)/ps])
%}

%e_qw_re = loadD([outKey,num2str(n),'_E_re_QW1.dat']);
%e_qw_im = loadD([outKey,num2str(n),'_E_im_QW1.dat']);
%E_qw_abs = abs(e_qw_re + 1i*e_qw_im);

% Find peaks
[pk,pk_i] = findpeaks(E_cav_abs,'sortstr','descend');
pk(3:end) = [];
pk_i(3:end) = [];
pk_i = sort(pk_i);

sizeIl = 200*fs;
sizeIr = 200*fs;

% Focus in on important part only
% Input to QW
x0 = t-t(pk_i(1));
indp = x0>=-sizeIl;
indm = x0<=sizeIr;
ind3 = indp==indm;
t0_qw = t(ind3);
y0_qw = E_cav(ind3);

sizeIl = 200*fs;
sizeIr = 6000*fs;

% Output from QW
x0 = t-t(pk_i(2));
indp = x0>=-sizeIl;
indm = x0<=sizeIr;
ind3 = indp==indm;
t1_qw = t(ind3);
y1_qw = E_cav(ind3);

%{
figure
plot(t/ps,E_cav_abs*um,'k-',...
     t0_qw/ps,abs(y0_qw)*um,'rs',...
     t1_qw/ps,abs(y1_qw)*um,'gx');
xlabel('t [ps]')
ylabel('E [V/um]')
grid on
legend('CAV')
title('QW: Input fields')
xlim([min(t)/ps, max(t)/ps])
%}

%{
ph = phase(E_cav.*exp(1i*w0*t))*180/pi;
figure
plot(t/ps,ph,'b-')
grid on
xlabel('t [ps]')
max(ph)-min(ph)
%}


% QW
%{
[F1_qw,Y1_qw] = getFT(t0_qw,y0_qw,W_LIM,W_NUM); %Input field
[F2_qw,Y2_qw] = getFT(t1_qw,y1_qw,W_LIM,W_NUM); %Result from QW

rat1 = Y2_qw./Y1_qw;


figure
plot(hbar*F2_qw/e,(abs(rat1)-1),'b-')
xlim(hbar*W_LIM/e)
yl =ylim;
ylim([0,yl(2)])
grid on



l0 = length(t0_qw);
l1 = length(t1_qw);

ld = floor((l1-l0)/2);

dt = t0_qw(2) - t0_qw(1);
[w,Y1_qw] = getFFT(y0_qw,dt,ld+0*floor(length(t1_qw)/2));
[w,Y2_qw] = getFFT(y1_qw,dt,0 +0*floor(length(t1_qw)/2));

figure
rat1 = Y2_qw./Y1_qw;
plot(hbar*w/e,(abs(rat1)-1),'b-')
xlim(hbar*W_LIM/e)
yl =ylim;
ylim([0,yl(2)])
grid on


asd
%}




%% Reflection spectrum from ABS
outKey = 'refSpecABS__';
n = 0;
t = loadD([outKey,num2str(n),'_t.dat']);
e_cav_re = loadD([outKey,num2str(n),'_E_re_CAVP_T0.dat']);
e_cav_im = loadD([outKey,num2str(n),'_E_im_CAVP_T0.dat']);

E_cav = (e_cav_re + 1i*e_cav_im).*exp(-1i*w0*t);
E_cav_abs = abs(E_cav);

%{
e_abs_re = loadD([outKey,num2str(n),'_E_re_ABS1.dat']);
e_abs_im = loadD([outKey,num2str(n),'_E_im_ABS1.dat']);
E_abs_abs = abs(e_abs_re + 1i*e_abs_im);

figure(10)
plot(t/ps,E_cav_abs*um,'k-',...
     t/ps,E_abs_abs*um,'r-')

 asd
%}

% Find peaks
[pk,pk_i] = findpeaks(E_cav_abs,'sortstr','descend');
pk(3:end) = [];
pk_i(3:end) = [];

sizeIl = 200*fs;
sizeIr = 200*fs;

% Focus in on important part only
% Input to QW
x0 = t-t(pk_i(1));
indp = x0>=-sizeIl;
indm = x0<=sizeIr;
ind3 = indp==indm;
t0_abs = t(ind3);
y0_abs = E_cav(ind3);

sizeIl = 200*fs;
sizeIr = 6000*fs;

% Output from QW
x0 = t-t(pk_i(2));
indp = x0>=-sizeIl;
indm = x0<=sizeIr;
ind3 = indp==indm;
t1_abs = t(ind3);
y1_abs = E_cav(ind3);

%{
figure
plot(t/ps,E_cav_abs*um,'k-',...
     t0_abs/ps,abs(y0_abs)*um,'rs',...
     t1_abs/ps,abs(y1_abs)*um,'gx');
xlabel('t [ps]')
ylabel('E [V/um]')
grid on
legend('CAV')
title('ABS: Input fields')
xlim([min(t)/ps, max(t)/ps])
%}


%[F1_abs,Y1_abs] = getFT(t0_abs,y0_abs,W_LIM,W_NUM); %Input field
%[F2_abs,Y2_abs] = getFT(t1_abs,y1_abs,W_LIM,W_NUM); %Result from QW
%{
[w,Y1_abs,Y2_abs] = getSpectrums(y0_abs,y1_abs,t0_abs(2)-t0_abs(1),t1_abs(2)-t1_abs(1));

rat2 = Y2_abs./Y1_abs;

figure
plot(hbar*w/e,(1-abs(rat2)),'b-')
xlim(hbar*W_LIM/e)
yl =ylim;
%ylim([0,yl(2)])
grid on

%cd('../pics')
%saveas(gcf,'absorbtion_abs_left','eps2c')

asd
%}


%% Process data

% QW
%[F1_qw,Y1_qw] = getFT(t0_qw,y0_qw,W_LIM,W_NUM); %Input field
%[F2_qw,Y2_qw] = getFT(t1_qw,y1_qw,W_LIM,W_NUM); %Result from QW
[w_qw,Y1_qw,Y2_qw] = getSpectrums(y0_qw,y1_qw,t0_qw(2)-t0_qw(1),t1_qw(2)-t1_qw(1));

rat1 = Y2_qw./Y1_qw;


%asd
%}

% ABSORBER
%[F1_abs,Y1_abs] = getFT(t0_abs,y0_abs,W_LIM,W_NUM); %Input field
%[F2_abs,Y2_abs] = getFT(t1_abs,y1_abs,W_LIM,W_NUM); %Result from QW
[w_abs,Y1_abs,Y2_abs] = getSpectrums(y0_abs,y1_abs,t0_abs(2)-t0_abs(1),t1_abs(2)-t1_abs(1));

rat2 = Y2_abs./Y1_abs;

%% CLEAN UP DATA
[rat1,w_qw]  = cleanUp(w_qw,rat1,W_LIM);
[rat2,w_abs] = cleanUp(w_abs,rat2,W_LIM);


%% Reflection spectrums
refSpecQW = (abs(rat1)-1);
refSpecABS = 1-abs(rat2);

%{
figure
plot(hbar*w/e,refSpecQW,'b-')
xlim(hbar*W_LIM/e)
grid on

cd('../matlab');
asd
%}


if length(w_qw)>length(w_abs)
   
    w = w_abs;
    refSpecQW = interp1(w_qw,refSpecQW,w)';
else
    w = w_qw;
    refSpecABS = interp1(w_abs,refSpecABS,w)';
end




[~,indW0] = min(abs(w-w0));
amplRef = (1-sqrt(loadD('refSpecABS__reflection_right.dat')));
disp(['Total absorbtion at w0 = ',num2str(100*refSpecABS(indW0),'%.6f'),' %'])
disp(['Output coupling        = ',num2str(100*amplRef),' %'])
amplAbs = refSpecABS(indW0)-amplRef;
disp(['Absorber               = ',num2str(100*amplAbs,'%.6f'),' %'])
disp(['Saturable ratio        = ',num2str(100*(amplAbs/refSpecABS(indW0))),' %'])

R = 0.5;
target_total = 2.4289; % Requested total
abs_base = 1.228874; % Base absorber absorbtion

abs_ratio     = target_total*R/abs_base;
output_target = (1-R)*target_total;
disp(['To get output = ',num2str(R)])
disp(['Set output    = ',num2str(1-output_target/100)])
disp(['Set abs       = ',num2str(abs_ratio)])

%{
figure
hold on
plot((2*pi*c0./w)/nm ,100*refSpecQW ,'b-')
hold off
%xlim(hbar*W_LIM/e)
%xlim([1.333, 1.533])
grid on

asd
%}

figure
diffGain = refSpecQW-refSpecABS;


diffGain_pks = max(diffGain);
[GAIN_WIDTH,GAIN_WIDTH_x0] = findAllGainWidths(diffGain,hbar*w/e);
hold on
plot(hbar*w/e ,100*refSpecQW ,'b-',...
     hbar*w/e,100*refSpecABS,'k-',...
     hbar*w0/e*[1,1],[0,max(refSpecQW)]*100,'k--',...
     hbar*w/e,diffGain*100,'r--')
hold off
xlim(hbar*W_LIM/e)
%xlim([1.333, 1.533])
xlabel('Energy [eV]')
ylabel('|R(w)| [\%]')
lg = legend('$$R_{qw}-1$$','$$1-R_{abs}$$','$$\hbar \omega_0$$','$$\Delta R$$');
set(lg, 'Interpreter', 'latex')
grid on

yl = get(gca,'ylim');
yl(2) = 10;
ylim([0,yl(2)])


abs1_deph = read_deph_scale('material_ABS1.config');
abs1_deph = abs1_deph/(2.0*Eb/hbar)/fs;

qw1_deph = read_deph_scale('material_QW1.config');
qw1_deph = qw1_deph/(2.0*Eb/hbar)/fs;


xl = xlim;
dy = yl(2);
yspace = linspace(yl(2)-0.55*dy,yl(2)-dy*0.05,7);
qw1_dens  = readDensity('material_QW1.config');
abs1_dens = readDensity('material_ABS1.config');
qw1_bg  = read_BG('material_QW1.config');
abs1_bg = read_BG('material_ABS1.config');
text(xl(1)+0.005,yspace(7),['$$n^{qw}$$ = $$' ,num2str( qw1_dens/1.0e16,'%.3f'),'\cdot 10^{16}$$ [$$m^{-2}$$]'],'BackgroundColor','w')
text(xl(1)+0.005,yspace(6),['$$n^{abs}$$ = $$',num2str(abs1_dens/1.0e14,'%.3f'),'\cdot 10^{14}$$ [$$m^{-2}$$]'],'BackgroundColor','w')
text(xl(1)+0.005,yspace(5),['$$E^{qw}_g$$ = ',num2str(qw1_bg,'%.3f'),' [eV]'],'BackgroundColor','w')
text(xl(1)+0.005,yspace(4),['$$E^{abs}_g$$ = ',num2str(abs1_bg,'%.3f'),' [eV]'],'BackgroundColor','w')
text(xl(1)+0.005,yspace(3),['$$\hbar \omega_0$$ = ',num2str(hbar*w0/e,'%.3f'),' [eV]'],'BackgroundColor','w')
text(xl(1)+0.005,yspace(2),['Gain Width = ',num2str(GAIN_WIDTH*1000,'%.0f '),' [meV]'],'BackgroundColor','w')
text(xl(1)+0.005,yspace(1),['Peak Net Gain  = ',num2str(diffGain_pks*100,'%.3f'), ' [\%]'],'BackgroundColor','w')

text_pos = (xl(1) + 0.65*(xl(2)-xl(1)));
yspace2 = [yspace(3), yspace(4)];
text(text_pos,yspace2(2),['$$\tau^{qw}_{deph}$$ = ',num2str(qw1_deph,'%.2f'), ' [fs]'],'BackgroundColor','w')
text(text_pos,yspace2(1),['$$\tau^{abs}_{deph}$$ = ',num2str(abs1_deph,'%.2f'), ' [fs]'],'BackgroundColor','w')



if ~exist('../pics','dir')
    mkdir('../pics')
end
oldDirPic = cd('../pics');
saveas(gcf,'plotSpectrum.eps','eps2c')


% Output data to files for later plotting requirements
fid = fopen('plot_eq_gain_region.txt','w');
for i = 1:length(w)
    fprintf(fid,'%.4e\t%.4e\n',hbar*w(i)/e, diffGain(i));
end
fclose(fid);

fid = fopen('plot_eq_gain.txt','w');
for i = 1:length(w)
    fprintf(fid,'%.4e\t%.4e\n',hbar*w(i)/e, refSpecQW(i));
end
fclose(fid);


fid = fopen('plot_eq_absorbtion.txt','w');
for i = 1:length(w)
    fprintf(fid,'%.4e\t%.4e\n',hbar*w(i)/e, refSpecABS(i));
end
fclose(fid);

%

%{
dE = Y2_qw./Y1_qw;

nw = dE;
w  = hbar*F1_qw/e;

figure
plot(w,real(nw),'b-')
grid on
title('re')


figure
plot(w,imag(nw),'b-')
grid on
title('im')

figure
plot(w,abs(nw),'b-')
grid on
title('abs')
asd
%}


%{
figure
diffGain = refSpecQW-refSpecABS;
[diffGain_pks,pks_i] = findpeaks(diffGain,'sortstr','descend');
GAIN_WIDTH = findAllGainWidths(diffGain,hbar*F1_qw/e);
hold on
plot((2*pi*c0./F1_qw)/nm ,refSpecQW*100 ,'b-',...
     (2*pi*c0./F1_abs)/nm,refSpecABS*100,'k-',...
     (2*pi*c0./F1_qw)/nm,(refSpecQW-refSpecABS)*100,'r--')
hold off
xlabel('Wavelength [nm]')
ylabel('Spectrum [\%]')
lg = legend('$$R_{qw}-1$$','$$1-R_{abs}$$','$$\hbar \omega_0$$','$$\Delta R$$');
set(lg, 'Interpreter', 'latex')
grid on


asd
%}




%% Dispersion: QW

phi = atan2(imag(rat1),real(rat1)); %Phase of complex signal
%phi = angle(rat1);
phi = unwrap(phi);    %Correct for 2pi phase jump
phi = smooth(phi);

if length(w_qw)>length(w_abs)
    w = w_abs;
    phi = interp1(w_qw,phi,w)';
else
    w = w_qw;
end
dw = w(2)-w(1);

if (numel(GAIN_WIDTH)==1)
   
    if (GAIN_WIDTH > 0)
       
        x_min = GAIN_WIDTH_x0;
        x_max = GAIN_WIDTH_x0 + GAIN_WIDTH;
        
        x_min = x_min - GAIN_WIDTH;
        x_max = x_max + GAIN_WIDTH;
    end
else
    x_min = xl(1);
    x_max = xl(2);
end

%[phi, F1_qw] = cleanUp(w , phi  , e*[x_min*0.5, 2*x_max]/hbar);

gd_qw   = (phi(2:end) - phi(1:end-1))/(dw);
gdd_qw  = (phi(3:end) - 2*phi(2:end-1) + phi(1:end-2))/(dw^2);
gddd_qw = (0.5*phi(5:end) - phi(4:end-1) + phi(2:end-3) - 0.5*phi(1:end-4))/(dw^3);
%const = (2*pi*c0./F1_qw(2:end-1)).^2;
%total_disp = gdd_qw.*const/(-2*pi*c0);

%Clean up picture based on FT
gd_x   = hbar*w(2:end  )/e;
gdd_x  = hbar*w(2:end-1)/e;
gddd_x = hbar*w(3:end-2)/e;

%{
figure(11)
plot(hbar*F1_qw/e,phi,'b-')
xlim([x_min, x_max])
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(phi),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('Phase')
saveas(gca,'dispersion-qw-order0','eps2c')
grid on
%}

%{
figure(12)
plot(gd_x,gd_qw/(fs),'b-')
%xlim([x_min, x_max])
xlim(hbar*DBR_LIM/e)
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(gd_qw/(fs)),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('GD [$$fs$$]')
title('QW: GD')
saveas(gca,'dispersion-qw-order1','eps2c')
grid on

figure(13)
gdd_qw = smooth(gdd_qw);
plot(gdd_x,gdd_qw/(fs*fs),'b-')
%xlim([x_min, x_max])
xlim(hbar*DBR_LIM/e)
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(gdd_qw/(fs*fs)),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('GDD [$$fs^2$$]')
title('QW: GDD')
grid on
saveas(gca,'dispersion-qw-order2','eps2c')
%}


%{
figure(113)
plot((2*pi*c0./(gdd_x*e/hbar))/nm ,gdd_qw/(fs*fs) ,'b-')
xlabel('Wavelength [nm]')
ylabel('GDD [$$fs^2$$]')
grid on
%xlim(2*pi*c0./([x_max, x_min]*e/hbar)/nm)
%yl = ylim;
xlim([950 ,1020])
ylim([-500, 200])
asd
%saveas(gca,'dispersion-qw-order2-Wavelength-WidthOutQW','eps2c')


fid = fopen('plot_Wavelength_GDD-WidthOutQW.txt','w');
x = (2*pi*c0./(gdd_x*e/hbar)); x = x(end:-1:1);
y = gdd_qw; y = y(end:-1:1);
for i = 1:length(x)
    fprintf(fid,'%.4f\t%.4e\n',x(i)/nm, y(i)/(fs*fs));
end
fclose(fid);


asd
%}

%{
figure(14)
gddd_qw = smooth(gddd_qw,10);
plot(gddd_x,gddd_qw/(fs*fs*fs),'b-')
%xlim([x_min, x_max])
xlim(hbar*DBR_LIM/e)
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(gddd_qw/(fs*fs*fs)),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('3rd Order [$$fs^3$$]')
title('QW: GDDD')
saveas(gca,'dispersion-qw-order3','eps2c')
grid on
%}

reflection = abs(rat1);


%% Dispersion: ABS

phi = atan2(imag(rat2),real(rat2)); %Phase of complex signal
phi = unwrap(phi);    %Correct for 2pi phase jump
%phi = smooth(phi);

if length(w_qw)>length(w_abs)
    w = w_abs;
else
    w = w_qw;
    phi = interp1(w_abs,phi,w)';
end
dw = w(2)-w(1);

%[phi, F1_qw] = cleanUp(F1_qw  , phi  , e*[x_min*0.5, 2*x_max]/hbar);

gd_abs   = (phi(2:end) - phi(1:end-1))/(dw);
gdd_abs  = (phi(3:end) - 2*phi(2:end-1) + phi(1:end-2))/(dw^2);
gddd_abs = (0.5*phi(5:end) - phi(4:end-1) + phi(2:end-3) - 0.5*phi(1:end-4))/(dw^3);
%const = (2*pi*c0./F1_qw(2:end-1)).^2;
%total_disp = gdd_abs.*const/(-2*pi*c0);

%Clean up picture based on FT
gd_x   = hbar*w(2:end  )/e;
gdd_x  = hbar*w(2:end-1)/e;
gddd_x = hbar*w(3:end-2)/e;

%{
figure(21)
plot(hbar*F1_qw/e,phi,'b-')
xlim([x_min, x_max])
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(phi),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('Phase')
saveas(gca,'dispersion-abs-order0','eps2c')
grid on
%}

%{
figure(22)
plot(gd_x,gd_abs/(fs),'b-')
%xlim([x_min, x_max])
xlim(hbar*DBR_LIM/e)
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(gd_abs/(fs)),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('GD [$$fs$$]')
title('ABS: GD')
saveas(gca,'dispersion-abs-order1','eps2c')
grid on

figure(23)
gdd_abs = smooth(gdd_abs);
plot(gdd_x,gdd_abs/(fs*fs),'b-')
%xlim([x_min, x_max])
xlim(hbar*DBR_LIM/e)
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(gdd_abs/(fs*fs)),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('GDD [$$fs^2$$]')
title('ABS: GDD')
saveas(gca,'dispersion-abs-order2','eps2c')
grid on

%{
figure(113)
plot((2*pi*c0./(gdd_x*e/hbar))/nm ,gdd_abs/(fs*fs) ,'b-')
xlabel('Wavelength [nm]')
ylabel('GDD [$$fs^2$$]')
grid on
%xlim(2*pi*c0./([x_max, x_min]*e/hbar)/nm)
%yl = ylim;
xlim([950 ,1020])
ylim([-500, 200])


figure(115)
plot(gdd_x ,gdd_abs/(fs*fs) ,'r--',...
     gdd_x ,gdd_qw/(fs*fs) ,'b-')
xlabel('Energy [eV]')
ylabel('GDD [$$fs^2$$]')
grid on
%xlim(2*pi*c0./([x_max, x_min]*e/hbar)/nm)
%yl = ylim;
xlim(hbar*DBR_LIM/e)
ylim([-500, 200])
legend('abs','gain')
hold on
plotLinearGainRegion(yl);
plot(hbar*w0*[1,1]/e,[-1,1]*max(gdd_qw/(fs*fs)),'k--')
hold off
saveas(gca,'dispersion-total-order2-wavelength','eps2c')

asd

figure(114)
plot((2*pi*c0./(gdd_x*e/hbar))/nm ,gdd_abs/(fs*fs) ,'r--',...
     (2*pi*c0./(gdd_x*e/hbar))/nm ,gdd_qw/(fs*fs) ,'b-')
xlabel('Wavelength [nm]')
ylabel('GDD [$$fs^2$$]')
grid on
%xlim(2*pi*c0./([x_max, x_min]*e/hbar)/nm)
%yl = ylim;
xlim([950 ,1020])
ylim([-500, 200])
legend('abs','gain')
hold on
plot((2*pi*c0/w0)*[1,1]/nm,[-1,1]*max(gdd_abs/(fs*fs)),'k--')
hold off


asd
%saveas(gca,'dispersion-qw-order2-Wavelength-WidthOutQW','eps2c')


fid = fopen('plot_Wavelength_GDD-WidthOutQW.txt','w');
x = (2*pi*c0./(gdd_x*e/hbar)); x = x(end:-1:1);
y = gdd_qw; y = y(end:-1:1);
for i = 1:length(x)
    fprintf(fid,'%.4f\t%.4e\n',x(i)/nm, y(i)/(fs*fs));
end
fclose(fid);


asd
%}



figure(24)
gddd_abs = smooth(gddd_abs,10);
plot(gddd_x,gddd_abs/(fs*fs*fs),'b-')
%xlim([x_min, x_max])
xlim(hbar*DBR_LIM/e)
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(gddd_abs/(fs*fs*fs)),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('3rd Order [$$fs^3$$]')
title('ABS: GDDD')
saveas(gca,'dispersion-abs-order3','eps2c')
grid on
%}

%% Dispersion: SUM

gd   = gd_qw + gd_abs;
gdd  = gdd_qw + gdd_abs;
gddd = gddd_qw + gddd_abs;

%{
figure(32)
plot(gd_x,gd/(fs),'b-')
%xlim([x_min, x_max])
xlim(hbar*DBR_LIM/e)
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(gd/(fs)),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('GD [$$fs$$]')
title('Total: GD')
saveas(gca,'dispersion-total-order1','eps2c')
grid on


figure(33)
gdd = smooth(gdd);
plot(gdd_x,gdd/(fs*fs),'b-')
%xlim([x_min, x_max])
xlim(hbar*DBR_LIM/e)
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(gdd/(fs*fs)),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('GDD [$$fs^2$$]')
title('Total: GDD')
saveas(gca,'dispersion-total-order2','eps2c')
grid on

figure(34)
gddd = smooth(gddd,10);
plot(gddd_x,gddd/(fs*fs*fs),'b-')
%xlim([x_min, x_max])
xlim(hbar*DBR_LIM/e)
yl = ylim;
hold on
plotLinearGainRegion(yl);
plot(hbar*w0/e*[1,1],[-1,1]*max(gddd/(fs*fs*fs)),'k--')
hold off
ylim(yl)
xlabel('Energy [eV]')
ylabel('3rd Order [$$fs^3$$]')
title('Total: GDDD')
saveas(gca,'dispersion-total-order3','eps2c')
grid on
%}

save('gdd_stuff.mat','reflection','gd_x','gd_qw','gdd_x','gdd_qw','gddd_x','gddd_qw','gd_abs','gdd_abs','gddd_abs','DBR_LIM','w0')


%close all

gdd_x = gdd_x*e/hbar;

figure
plot(hbar*gdd_x/e ,gdd_abs/(fs*fs) ,'r--',...
     hbar*gdd_x/e ,gdd_qw/(fs*fs) ,'b-',...
     hbar*gdd_x/e ,(gdd_qw+gdd_abs)/(fs*fs) ,'m-')
xlabel('Energy [eV]')
ylabel('GDD [$$fs^2$$]')
grid on
%xlim(2*pi*c0./([x_max, x_min]*e/hbar)/nm)
%yl = ylim;
xlim(hbar*DBR_LIM/e)
ylim([-250, 250])
hold on
%plotLinearGainRegion(yl);
plot(hbar*w0*[1,1]/e,[-1,1]*max(gdd_qw/(fs*fs)),'k--')
hold off
legend('abs','gain','sum','w0')

saveas(gca,'dispersion-total-order2-Zoom-energy','eps2c')

gdd_x_w = 2.0*pi*c0./gdd_x;


figure
plot(gdd_x_w/nm ,gdd_abs/(fs*fs) ,'r--',...
     gdd_x_w/nm ,gdd_qw/(fs*fs) ,'b-',...
     gdd_x_w/nm ,(gdd_qw+gdd_abs)/(fs*fs) ,'m-')
xlabel('Wavelength [nm]')
ylabel('GDD [$$fs^2$$]')
grid on
%xlim(2*pi*c0./([x_max, x_min]*e/hbar)/nm)
%yl = ylim;
ylim([-250, 250])
xlim([990,1060])

yl = ylim;

hold on
plot((2*pi*c0/w0)*[1,1]/nm,yl,'k--')
%plot(mqw2222(:,1),mqw2222(:,2),'ks--')
hold off
legend('abs','gain','sum','\lambda_0')

saveas(gca,'dispersion-total-order2-Zoom-wavelength','eps2c')


w_gain = w_qw;
R_gain = abs(rat1).^2;
R_abs  = abs(rat2).^2;

save('gdd_pulsed_passive.mat','w_gain','w_abs','R_gain','R_abs','refSpecQW','refSpecABS','gdd_x','gdd_x_w','gdd_abs','gdd_qw');




cd(oldDirPic)

cd(oldDir)

end

function [x,y] = d1(X,Y)
%% First difference

y = (Y(2:end) - Y(1:end-1))./(2*(X(2:end) - X(1:end-1)));
x = (X(2:end) + X(1:end-1))/2;


end


function v = loadD(name)
% Read a single double or a list of doubles into v
% v is a row vector

fid = fopen(name,'rb');
v = fread(fid,'double');
fclose(fid);

end

function a = loadDA(name,size)
% Read an 2d array of doubles 
% Name is name of file
% Size is [time-steps, K-steps]

if numel(size)~=2
   disp('loadDA(): Cannot read more or less than 2d arrays')
   name
   size
   asd
end

fid = fopen(name);
a = fread(fid,size(2:-1:1),'double');
fclose(fid);

a = a';

end


function [F,Y] = getFT(t,y,WLIM,N)

N_w   = 6000;
w_min = 6.87638e+14;
w_max = 2.98899e+15;

global hbar;
global e;

%F = linspace(w_min, w_max, N_w);
%F = linspace(1.3, 1.6, 600)*e/hbar;
F = linspace(WLIM(1), WLIM(2), N);
Y = zeros(1,length(F));

dt_max = 1/(2*max(F));
disp(['Nyquist-Shannon limit: dt < ',num2str(dt_max/1.0e-15,'%.3f'),' [fs]'])

dt = t(2)-t(1);

for i = 1:length(F)
    for j = 1:length(y)

        Y(i) = Y(i) + y(j)*exp(1i*F(i)*t(j))*dt;

    end
end


end


function [w,Y1,Y2] = getSpectrums(y1,y2,dt1,dt2)
%% Get spectrums of two signals sampled at two rates
% The joint spectrums will be on the same axis


l0 = length(y1);
l1 = length(y2);
ld = floor(abs(l1-l0)/2);

maxL = floor(max([length(y1),length(y2)])/2);

if (length(y1)>length(y2))
   ld1 = 0;
   ld2 = ld;
else
    ld1 = ld;
    ld2 = 0;
end

[w1,Y1] = getFFT(y1,dt1,ld1 +maxL);
[w,Y2] = getFFT(y2,dt2,ld2 +maxL);

% Interpolate w1,Y1 over to w
Y1 = transpose(interp1(w1,Y1,w));


end


function [w,y] = getFFT(y,dt,PAD)
%% Find fourier transform of y sampled with timestep dt
% PAD offers options to zeropad the data for increased resolution and to
% compare with other signals of different length
% 

Fs = 1/dt;  % Sampling frequency
y = [zeros(PAD,1);y;zeros(PAD,1)];

L = length(y);        
NFFT = 2^nextpow2(L); %Zero padding of signal for speed
Y = fftshift(fft(y,NFFT)/L); % Move data to correct ordering
f = Fs/2*linspace(-1,0,NFFT/2+1); % Definition of DFT in matlab means my freq is on negative side
y = sqrt(2)*Y(1:NFFT/2+1);

w = -2*pi*f; % Go to omega and reverse frequency axis

w = w(end:-1:1);
y = y(end:-1:1);

end


function [F,Y] = getIFT(t,y,WLIM,N)

global hbar;
global e;

%F = linspace(w_min, w_max, N_w);
%F = linspace(1.3, 1.6, 600)*e/hbar;
F = linspace(WLIM(1), WLIM(2), N);
Y = zeros(1,length(F));

dt_max = 1/(2*max(F));
disp(['Nyquist-Shannon limit: dt < ',num2str(dt_max/1.0e-15,'%.3f'),' [fs]'])

dt = t(2)-t(1);

for i = 1:length(F)
    for j = 1:length(y)

        Y(i) = Y(i) + y(j)*exp(-1i*F(i)*t(j))*dt;

    end
end


end


function [width,width_t0] =  findAllGainWidths(diff,x)
% widht is width of interval
% width_t0 is starting point of each interval
% x_lim will limit where the points are checked

n = length(diff);
% Find Number of islands
width = [];
width_t0 = [];
count = 1;
i0 = 0;
for i = 1:n
    
    if ((diff(i)>0)&&(i0==0))
        i0 = i;
    end
    
    if ((diff(i)<=0)&&(i0>0))
        
        width(count) = x(i)-x(i0);
        width_t0(count) = x(i0);
        count = count + 1;
        i0 = 0;
    end
end


end

function [Y,X] = cleanUp(x, y, XL)

x0 = XL(1);
x1 = XL(2);

ind_m = x <= x0;
ind_p = x >= x1;
ind = ind_m==ind_p;

X = x(ind);
Y = y(ind);

end

function width = plotLinearGainRegion(YL)

width = 0;

if (exist('plot_eq_gain_region.txt','file'))
    
   ymin = YL(1);
   ymax = YL(2);
   gain_lin = load('plot_eq_gain_region.txt');
   ind = gain_lin(:,2)<0;
   gain_lin(ind,:) = []; %Remove negative
   if numel(gain_lin(:,2))>0
       gain_lin = [gain_lin(1,1), gain_lin(end,1)];
       
       width = gain_lin(2)-gain_lin(1);
       if (width > 0)
           H = rectangle('Position',[gain_lin(1),ymin,gain_lin(2)-gain_lin(1),ymax-ymin],'FaceColor',0.9*[1,1,1],'EdgeColor','None');
           uistack(H,'bottom');
       end
   end
end

end


function dens = readDensity(file_name)

fid = fopen(file_name,'r');

num = 0; cnt = 1;

l = fgetl(fid);
while (l ~= -1)
    
    % Clean up line
    % Scan for comment and remove comment
    ind = strfind(l,'//');
    l(ind:end) = [];
    
    % Remove all whitespace
    ind = strfind(l,' ');
    l(ind) = [];
    ind = strfind(l,'\t');
    l(ind) = [];
    
    ind = strfind(l,'\n');
    l(ind) = [];
    
    
    num = str2num(l);
    
    %Count of number of lines to the density
    if (cnt == 13)
        break;
    end
    
    cnt = cnt + 1;
    % Get next line
    l = fgetl(fid);
end

dens = num;

fclose(fid);
end


function bg = read_BG(file_name)

fid = fopen(file_name,'r');

num = 0; cnt = 1;

l = fgetl(fid);
while (l ~= -1)
    
    % Clean up line
    % Scan for comment and remove comment
    ind = strfind(l,'//');
    l(ind:end) = [];
    
    % Remove all whitespace
    ind = strfind(l,' ');
    l(ind) = [];
    ind = strfind(l,'\t');
    l(ind) = [];
    
    ind = strfind(l,'\n');
    l(ind) = [];
    
    
    num = str2num(l);
    
    %Count of number of lines to the density
    if (cnt == 12)
        break;
    end
    
    cnt = cnt + 1;
    % Get next line
    l = fgetl(fid);
end

bg = num;

fclose(fid);
end

function bg = read_deph_scale(file_name)

fid = fopen(file_name,'r');

num = 0; cnt = 1;

l = fgetl(fid);
while (l ~= -1)
    
    % Clean up line
    % Scan for comment and remove comment
    ind = strfind(l,'//');
    l(ind:end) = [];
    
    % Remove all whitespace
    ind = strfind(l,' ');
    l(ind) = [];
    ind = strfind(l,'\t');
    l(ind) = [];
    
    ind = strfind(l,'\n');
    l(ind) = [];
    
    
    num = str2num(l);
    
    %Count of number of lines to the deph_scale
    if (cnt == 10)
        break;
    end
    
    cnt = cnt + 1;
    % Get next line
    l = fgetl(fid);
end

bg = num;

fclose(fid);
end



function [n_s, z0_s, z1_s] = loadStruct(fileName)
%% Plot structure in background
% Input:
% runKey - What to load from
% Output:
% n_s  - Refractive index in layer s 
% z0_s - Position of Left edge of layer s
% z1_s - Position of Right edge of layer s
% qw_pos - Position of QW's

str = load(fileName);
[n,m] = size(str);

% Import refractive indices and positions
z0_s = [];
z1_s = [];
n_s = [];
qw_z0_ind = [];
cnt = 1;
cnt2 = 1;
for i = 1:n
    
        
    if (str(i,1)==2)    % QW
        qw_z0_ind(cnt2) = cnt;
        cnt2 = cnt2 + 1;
    end
    
    if (str(i,1) == 1)  % CAVITY
        
        z0_s(cnt) = str(i,2);
        z1_s(cnt) = str(i,3);
        n_s(cnt)  = str(i,4);
        
        cnt = cnt + 1;
    end
end

%Find and remove air layers, replace them with 2um of air
ind = find(n_s==1);
if (numel(ind) <= length(n_s)/2)
    z0_s(ind) = [];
    z1_s(ind) = [];
    n_s(ind) = [];

    for i = 1:length(qw_z0_ind)

        % Find how many layers before qw_z0_ind(i) is deleted
        num = sum(ind<qw_z0_ind(i));
        qw_z0_ind(i) = qw_z0_ind(i) - num;

    end


    AIR_WIDTH = 1e-6;
    if (z0_s(1)>AIR_WIDTH)

        z1_s = [z0_s(1)          , z1_s];
        z0_s = [z0_s(1)-AIR_WIDTH, z0_s];
        n_s  = [1, n_s];

        % All layers are infront of the first one
        qw_z0_ind = qw_z0_ind + 1;

    else
        z0_s = [z0_s, z1_s(end)];
        z1_s = [z1_s, z0_s(end)+AIR_WIDTH];
        n_s  = [n_s,1];
    end
end

width = z1_s - z0_s;




end

function fwhm = findFWHM(DT,pulse)
% Find FWHM of central pulse. (For interfering pulses)
% Input is a pulse with time difference DT. That is... DT is refrenced from
% 0

[maxp, maxi] = max(pulse); % max pulse

% Find left boundary
i0 = maxi;
while ((i0 > 2)&&(pulse(i0-1) > 0.5*maxp))
    
    i0 = i0 -1;
end

% Find right boundary
i1 = maxi;
while ((i1 < length(pulse))&&(pulse(i1+1) > 0.5*maxp))
    
    i1 = i1 +1;
end
%fwhm = DT(i1) - DT(i0);


% Interpolate data
ind0 = i0-1:i0+1;
T0 = invertYData(DT(ind0),pulse(ind0),0.5*maxp);

ind1 = i1-1:i1+1;
T1 = invertYData(DT(ind1),pulse(ind1),0.5*maxp);


fwhm = abs(T1 - T0);
end

function X1 = invertYData(X0,Y0,Y1)

Xscale = max(X0)-min(X0);
X0 = X0/Xscale;

Yscale = max(Y0)-min(Y0);
Y0 = Y0/Yscale;
Y1 = Y1/Yscale;

guess = mean(X0);
f = @(x) abs(interp1(X0,Y0,x)-Y1);
[X1,~,flag] = fminsearch(f,guess); %#ok<NASGU>

%{
if (flag~=1)
    
    
    
    flag
    guess/1.0e-12
    X0/1.0e-12
    Y0
    
    figure
    plot(X0/1.0e-12,Y0,'b-',...
         X1/1.0e-12,Y1,'ro')
    grid on
    asd
end
%}

X1 = X1*Xscale;

end
