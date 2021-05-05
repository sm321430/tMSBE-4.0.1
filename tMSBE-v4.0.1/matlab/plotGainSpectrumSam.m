%%Function for plotting 1D or transverse pump-probe simulation
%S.A.Mclaren (June 2020)

function plotGainSpectrumSam()

close all
clear all

IO_reflectionPlots=1;

setupPlot;
setupConstants;

% fs = 1.0e-15;
% ps = 1.0e-12;
% um = 1.0e-6;
% nm = 1.0e-9;
% 
% global hbar;
% global e;
% global c0;
% 
% hbar = 1.054589e-34;
% e = 1.602189e-19;
% c0   = 2.99792458E+08;
% eps0 = 8.854187817620E-12;
% 
% m0   = 9.109389754e-31;
% me   = 0.06085*m0;
% mh   = 0.234*m0;
% a0   = 1.062146e-08;
% mr = me*mh/(me+mh);
% Eb = (hbar^2 / (2*mr*a0*a0));
% 
% 
% set(0,'defaulttextinterpreter','tex') %Default
% %set(0,'defaulttextinterpreter','latex')
% set(0,'defaultAxesFontName', 'Arial')
% set(0,'defaultAxesFontsize', 40)
% set(0,'defaultAxesLinewidth', 4)
% set(0,'defaultTextFontsize', 40)
% set(0,'defaultlineMarkerSize',40)
% set(0,'defaultlinelinewidth', 4) %Thin lines
% maps %Load maps file for color scheme
% POS=[1,1,1200,800];
% POS2=[1,1,500,900];
% %POS_plot=[0.14,0.1100,0.53125,0.85];
% %POS_plot=[0.1330,0.1100,0.5930,0.8150];
% my_lineStyle={'-kd','-bo','-r+', '-gx', '-c*','-mo','-k+','-kx','-k*','-b+','-bx','-b*','-g+'};

reflectLim=[2,6]; %Limits for reflection plots
outKey ='/Volumes/SAMbackup/tMSBE-ECAV-data-2021/tMSBE-v3.8-ECAV-04-cyclicSquareCav-OGSGS-1D-theta2-spontEmis-wExpSBE-focus10/run/';
date='032221';
test='tMSBE-ECAV-04-cyclicSquareCav-OGSGS-1D-theta2-spontEmis-wExpSBE-focus10';
outKeyQW = [outKey,'refSpecQW__'];
outKeyABS = [outKey,'refSpecABS__'];
%outKey = '../run/refSpecQW__';
location='CAVM'; %Field location for uploading and saving
point='T0';
plot_num =1; %Output number +1
test_folder='test';
saveKey_local='Fall2020-Summer2021/RingCav/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey = [saveKey_fold,'/',test,'-out',num2str(plot_num-1),'-',location,'-'];
user_entry = input(['saveKey=',saveKey,'? (y to continue): '], 's');
if user_entry~='y'
    gafuggle
end
if exist(saveKey_fold,'dir')~=7
    mkdir(saveKey_fold);
end

OC_loss=0.0; %Artificial constant output coupling loss

w0 = loadD([outKeyQW,'w0.dat']);
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])
disp(['lambda   = ',num2str(2*pi*c0/w0/nm,'%.2f'),' [nm]'])

% Set Frequency limits
W_LIM = w0+e*0.1/hbar*[-1.0,1.0];
W_NUM = 600;

%Load plot limits (GD, GDD, GDDD)
n_s = loadStruct([outKeyQW,'system_structure.dat']);
%DBRind = n_s(1:2);
DBRind = [2.9459, 3.4345];
% Approx width of DBR
DF = w0*(4/pi)*asin((min(DBRind)-max(DBRind))/(sum(DBRind)));
fa = w0 + DF/2;
fb = w0 - DF/2;
DBR_LIM = [fa,fb];
%DBR_LIM = W_LIM;

%% Plot all electric fields

%% Reflection spectrum from QW
n = 0;
t = loadD([outKeyQW,num2str(n),'_t.dat']);
e_cav_re = loadD([outKeyQW,num2str(n),'_E_re_',location,'_',point,'.dat']);
e_cav_im = loadD([outKeyQW,num2str(n),'_E_im_',location,'_',point,'.dat']);

E_cav = (e_cav_re + 1i*e_cav_im).*exp(-1i*w0*t);
E_cav_abs = abs(E_cav);

E_cav_qw=E_cav;

% Findpeaks
[pk,pk_i] = findpeaks(E_cav_abs,'sortstr','descend');
pk(3:end) = [];
pk_i(3:end) = [];
pk_i = sort(pk_i);

sizeIl = 400*fs;
sizeIr = 3000*fs;

% Focus in on important part only
% Input to QW

x0 = t-t(pk_i(1));
indp = x0>=-sizeIl;
indm = x0<=sizeIr;
ind3 = indp==indm;
t0_qw = t(ind3);
y0_qw = E_cav(ind3);

sizeIl = 400*fs;
sizeIr = 3000*fs;

% Output from QW
try
x0 = t-t(pk_i(2));
catch
   x0=t-t(pk_i(1));
end
indp = x0>=-sizeIl;
indm = x0<=sizeIr;
ind3 = indp==indm;
t1_qw = t(ind3);
y1_qw = E_cav(ind3);

%% Reflection spectrum from ABS
n = 0;
t = loadD([outKeyABS,num2str(n),'_t.dat']);
e_cav_re = loadD([outKeyABS,num2str(n),'_E_re_CAVM_T0.dat']);
e_cav_im = loadD([outKeyABS,num2str(n),'_E_im_CAVM_T0.dat']);

E_cav = (e_cav_re + 1i*e_cav_im).*exp(-1i*w0*t);
E_cav_abs = abs(E_cav);

% Find peaks
[pk,pk_i] = findpeaks(E_cav_abs,'sortstr','descend');
pk(3:end) = [];
pk_i(3:end) = [];

sizeIl = 400*fs;
sizeIr = 2000*fs;

% Focus in on important part only
% Input to QW

x0 = t-t(pk_i(1));
indp = x0>=-sizeIl;
indm = x0<=sizeIr;
ind3 = indp==indm;
t0_abs = t(ind3);
y0_abs = E_cav(ind3);

sizeIl = 400*fs;
sizeIr = 3000*fs;

% Output from QW
try
x0 = t-t(pk_i(2));
catch
   x0=t-t(pk_i(1));
end
indp = x0>=-sizeIl;
indm = x0<=sizeIr;
ind3 = indp==indm;
t1_abs = t(ind3);
y1_abs = E_cav(ind3);

%% Process data

% QW
[w_qw,Y1_qw,Y2_qw] = getSpectrums(y0_qw,y1_qw,t0_qw(2)-t0_qw(1),t1_qw(2)-t1_qw(1));

rat1 = Y2_qw./Y1_qw;
% figure
% plot(hbar*w_qw/e,abs(Y1_qw))
% xlim([1.16, 1.26])
% 
% figure
% plot(hbar*w_qw/e,abs(Y2_qw))
% xlim([1.16, 1.26])
% 
% 
% figure
% plot(hbar*w_qw/e,100*(abs(rat1)-1))
% xlim([1.16, 1.26])
% asd

% ABSORBER
[w_abs,Y1_abs,Y2_abs] = getSpectrums(y0_abs,y1_abs,t0_abs(2)-t0_abs(1),t1_abs(2)-t1_abs(1));

rat2 = Y2_abs./Y1_abs;

%% CLEAN UP DATA
[rat1,w_qw]  = cleanUp(w_qw,rat1,W_LIM);
[rat2,w_abs] = cleanUp(w_abs,rat2,W_LIM);

%% Reflection spectrums
refSpecQW = (abs(rat1)-1);
refSpecABS = 1-abs(rat2);
%refSpecAbs=abs(rat2);

% figure
% semilogy(t/ps,E_cav_abs)
% 
% figure
% plot(hbar*w_abs/e,100*abs(rat2))
% 
% figure
% plot(hbar*w_abs/e,angle(rat2))
% asd

if length(w_qw)>length(w_abs)
w = w_abs;
    refSpecQW = interp1(w_qw,refSpecQW,w)';
else
    w = w_qw;
    refSpecABS = interp1(w_abs,refSpecABS,w)';
end

refSpecABS = refSpecABS+OC_loss/100; %Artificial constant addditional OC loss

[~,indW0] = min(abs(w-w0));
amplRef = (1-sqrt(loadD([outKeyABS,'reflection_right.dat'])));
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

diffGain = refSpecQW-refSpecABS;
[GAIN_WIDTH,GAIN_WIDTH_x0] = findAllGainWidths(diffGain,hbar*w/e);
diffGain_pks = max(diffGain);

if IO_reflectionPlots==1
    tmp_fig=figure;
    set(tmp_fig,'Position',POS);
    hold on
    plot(hbar*w/e ,100*refSpecQW ,'b-',...
     hbar*w/e,100*refSpecABS,'k-',...
     hbar*w0/e*[1,1],[min(refSpecQW),max(refSpecQW)]*100,'k--',...
     hbar*w/e,diffGain*100,'r--')
    hold off
    xlim(hbar*W_LIM/e)
    xlim([1.14, 1.27])
    ylim(reflectLim)
    xlabel('Energy [eV]')
    ylabel('|R(w)| [%]')
    lg = legend('R_{qw}-1','1-R_{abs}','\omega_0','\Delta R');
    grid on
    saveas(tmp_fig,[saveKey,'ReflectionSpectrum.png']);
    yl = get(gca,'ylim');


abs1_deph = read_deph_scale([outKey,'material_ABS1.config']);
abs1_deph = abs1_deph/(2.0*Eb/hbar)/fs;

qw1_deph = read_deph_scale([outKey,'material_QW1.config']);
qw1_deph = qw1_deph/(2.0*Eb/hbar)/fs;


xl = xlim;
dy = yl(2);
%yspace = linspace(yl(2)-0.55*dy,yl(2)-dy*0.05,7);
qw1_dens  = readDensity([outKey,'material_QW1.config']);
abs1_dens = readDensity([outKey,'material_ABS1.config']);
qw1_bg  = read_BG([outKey,'material_QW1.config']);
abs1_bg = read_BG([outKey,'material_ABS1.config']);
% text(xl(1)+0.005,yspace(7),['$$n^{qw}$$ = $$' ,num2str( qw1_dens/1.0e16,'%.3f'),'\cdot 10^{16}$$ [$$m^{-2}$$]'],'BackgroundColor','w')
% text(xl(1)+0.005,yspace(6),['$$n^{abs}$$ = $$',num2str(abs1_dens/1.0e14,'%.3f'),'\cdot 10^{14}$$ [$$m^{-2}$$]'],'BackgroundColor','w')
% text(xl(1)+0.005,yspace(5),['$$E^{qw}_g$$ = ',num2str(qw1_bg,'%.3f'),' [eV]'],'BackgroundColor','w')
% text(xl(1)+0.005,yspace(4),['$$E^{abs}_g$$ = ',num2str(abs1_bg,'%.3f'),' [eV]'],'BackgroundColor','w')
% text(xl(1)+0.005,yspace(3),['$$\hbar \omega_0$$ = ',num2str(hbar*w0/e,'%.3f'),' [eV]'],'BackgroundColor','w')
% text(xl(1)+0.005,yspace(2),['Gain Width = ',num2str(GAIN_WIDTH*1000,'%.0f '),' [meV]'],'BackgroundColor','w')
% text(xl(1)+0.005,yspace(1),['Peak Net Gain  = ',num2str(diffGain_pks*100,'%.3f'), ' [\%]'],'BackgroundColor','w')

%text_pos = (xl(1) + 0.65*(xl(2)-xl(1)));
%yspace2 = [yspace(3), yspace(4)];
%text(text_pos,yspace2(2),['$$\tau^{qw}_{deph}$$ = ',num2str(qw1_deph,'%.2f'), ' [fs]'],'BackgroundColor','w')
%text(text_pos,yspace2(1),['$$\tau^{abs}_{deph}$$ = ',num2str(abs1_deph,'%.2f'), ' [fs]'],'BackgroundColor','w')
else
    xl=[1.14,1.27];
end

%% Dispersion: QW

phi = atan2(imag(rat1),real(rat1)); %Phase of complex signal
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

%Clean up picture based on FT
gd_x   = hbar*w(2:end  )/e;
gdd_x  = hbar*w(2:end-1)/e;
gddd_x = hbar*w(3:end-2)/e;

reflection = abs(rat1);


%% Dispersion: ABS

phi = atan2(imag(rat2),real(rat2)); %Phase of complex signal
phi = unwrap(phi);    %Correct for 2pi phase jump
phi = smooth(phi);

if length(w_qw)>length(w_abs)
    w = w_abs;
else
    w = w_qw;
    phi = interp1(w_abs,phi,w)';
end
dw = w(2)-w(1);

gd_abs   = (phi(2:end) - phi(1:end-1))/(dw);
gdd_abs  = (phi(3:end) - 2*phi(2:end-1) + phi(1:end-2))/(dw^2);
gddd_abs = (0.5*phi(5:end) - phi(4:end-1) + phi(2:end-3) - 0.5*phi(1:end-4))/(dw^3);

%Clean up picture based on FT
gd_x   = hbar*w(2:end  )/e;
gdd_x  = hbar*w(2:end-1)/e;
gddd_x = hbar*w(3:end-2)/e;

%% Dispersion: SUM

gd   = gd_qw + gd_abs;
gdd  = gdd_qw + gdd_abs;
gddd = gddd_qw + gddd_abs;
gdd_x = gdd_x*e/hbar;

save([saveKey,'gdd_stuff.mat'],'reflection','gd_x','gd_qw','gdd_x','gdd_qw','gddd_x','gddd_qw','gd_abs','gdd_abs','gddd_abs','DBR_LIM','w0')

if IO_reflectionPlots==1
tmp_fig=figure;
  set(tmp_fig,'Position',POS);
plot(hbar*gdd_x/e ,(gdd_abs)/(fs*fs) ,'r--',...
     hbar*gdd_x/e ,(gdd_qw)/(fs*fs) ,'b-',...
     hbar*gdd_x/e ,(gdd_qw+gdd_abs)/(fs*fs) ,'m-')
xlabel('Energy [eV]')
ylabel('GDD [fs^2]')
grid on
%xlim(2*pi*c0./([x_max, x_min]*e/hbar)/nm)
%yl = ylim;
xlim(hbar*DBR_LIM/e)
ylim([-1000, 1000])
hold on
%plotLinearGainRegion(yl);
plot(hbar*w0*[1,1]/e,[-1,1]*max(gdd_qw/(fs*fs)),'k--')
hold off
legend('abs','gain','sum','\omega_0')
saveas(tmp_fig,[saveKey,'GDD.png']);
end

gdd_x_w = 2.0*pi*c0./gdd_x;

if IO_reflectionPlots==1
tmp_fig=figure;
  set(tmp_fig,'Position',POS);
plot(gdd_x_w/nm ,gdd_abs/(fs*fs) ,'r--',...
     gdd_x_w/nm ,gdd_qw/(fs*fs) ,'b-',...
     gdd_x_w/nm ,(gdd_qw+gdd_abs)/(fs*fs) ,'m-')
xlabel('Wavelength [nm]')
ylabel('GDD [fs^2]')
POS=[1,1,1200,800];
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
saveas(tmp_fig,[saveKey,'GDDlambda.png']);
end

w_gain = w_qw;
R_gain = abs(rat1).^2;
R_abs  = abs(rat2).^2;

save([saveKey,'gdd_pulsed_passive.mat'],'E_cav_qw','w_gain','w_abs','R_gain','R_abs','refSpecQW','refSpecABS','gdd_x','gdd_x_w','gdd_abs','gdd_qw');

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
