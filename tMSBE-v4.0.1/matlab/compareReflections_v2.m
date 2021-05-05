%%Compare 1D pump-probe tests from TMSBE code
%%S.A.McLaren August,2020

clear all
close all

%% Simulation flags
IO_reflectionPlot=1; %Plot reflection spectrum comparisons
IO_energyDispersionPlot=1; %Plot energy GDD comparison
IO_wavelengthDispersionPlot=1; %Plot wavelength GDD comparison
IO_save=1; %Save plots

%Setup
setupPlot
setupConstants

location='CAVM'; %Field location for uploading and saving
OC_loss=1.0; %Artificial constant output coupling loss
point='T0';
date='031921';
test='RCAV-n2p5-n2p9';
test_folder='test';
saveKey_local='Fall2020-Summer2021/RingCAV/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey = [saveKey_fold,'/',date,'-',test,'-out','-',location,'-'];
user_entry = input(['saveKey=',saveKey,'? (y to continue): '], 's');
if user_entry~='y'
    gafuggle
end
if exist(saveKey_fold,'dir')~=7
    mkdir(saveKey_fold);
end

%% Plot limits and parameters
yl=[1.2,3.5]; %Y-limits for reflection spectrum plots (Reflection %)
xl=[1.15,1.25]; %X-limits for reflection spectrum plots (eV)
yl_gdd=[-1000, 1000]; %Y-limit for GDD

%% File locations and legend entries
% compReflect(1).outKey='../../tMSBE-v3.7-VCAV72-VCAV64-1dReflection/run/';
% compReflect(2).outKey='../../tMSBE-v3.7-VCAV71-VCAV63-1DReflection/run/';
% compReflect(3).outKey='../../tMSBE-v3.7-VCAV73-VCAV68-1dReflection/run/';
% compReflect(4).outKey='../../tMSBE-v3.7-VCAV77-VCAV75-1dReflection/run/';
% compReflect(5).outKey='../../tMSBE-v3.7-VCAV74-VCAV69-1dReflection/run/';
% compReflect(6).outKey='../../tMSBE-v3.7-VCAV78-VCAV76-1dReflection/run/';
% 
% compReflect(1).legend='\theta=2^\circ';
% compReflect(2).legend='\theta=4^\circ';
% compReflect(3).legend='\theta=8^\circ';
% compReflect(4).legend='\theta=12^\circ';
% compReflect(5).legend='\theta=16^\circ';
% compReflect(6).legend='\theta=20^\circ';

% compReflect(1).outKey='/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV88-UncoatedChip-1dReflection-theta0-Ny288/run/';
% compReflect(2).outKey='/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV87-1dReflection-theta0-Ny288/run/';
% 
% compReflect(1).legend='Uncoated';
% compReflect(2).legend='Ideal coating';

%compReflect(1).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.7-RCAV6-1d-n2p5-CW/run/';
%compReflect(2).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.7-RCAV6-1d-n2p5-CCW-fullOccNum/run/';
% compReflect(1).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021//tMSBE-v3.7-RCAV6-1D-n2p5-colThresh-2em2-theta0-6400lam-spointEmis-wExpSBE-CCWreflect/run/';
% compReflect(2).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.7-LCAV1-1D-n2p5-reflect1/run/'; 
% compReflect(3).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.7-VCAV93-1d-n2p5-reflect1/run/'; 
% 
% compReflect(1).legend='RCAV';
% %compReflect(2).legend='CCW';
% compReflect(2).legend='LCAV';
% compReflect(3).legend='VCAV';
% 
% compReflect(1).pass=1;
% compReflect(2).pass=1;
% %compReflect(3).pass=1;
% compReflect(3).pass=2;

% compReflect(1).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV9-debug-LCAV/run/';
% compReflect(2).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV9-debug-RCAV-broken/run/';
% %compReflect(3).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV9-debug-RCAV-fixed/run/';
% compReflect(1).legend='LCAV';
% compReflect(2).legend='RCAV-old';
% %compReflect(3).legend='RCAV-fixed';
% compReflect(1).pass=1;
% compReflect(2).pass=1;
% %compReflect(3).pass=1;
% compReflect(1).OC_loss=1;
% compReflect(2).OC_loss=0;
%compReflect(3).OC_loss=0;

compReflect(3).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV92-RCAV84reflection/run/';
compReflect(2).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV93-RCAV45reflection/run/';
compReflect(1).outKey='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV94-RCAV27reflection/run/';
compReflect(1).legend='n_0=2.5e16m^{-2}';
compReflect(2).legend='n_0=2.7e16m^{-2}';
compReflect(3).legend='n_0=2.9e16m^{-2}';
compReflect(1).pass=1;
compReflect(2).pass=1;
compReflect(3).pass=1;
compReflect(1).OC_loss=0;
compReflect(2).OC_loss=0;
compReflect(3).OC_loss=0;

num_sims=length(compReflect);

for j=1:length(compReflect)
   outKeyQW=[compReflect(j).outKey,'refSpecQW__'];
   outKeyABS=[compReflect(j).outKey,'refSpecABS__'];
   outKeyMat=[compReflect(j).outKey,'material_'];
   
   w0 = loadD([outKeyQW,'w0.dat']);
   W_LIM = w0+e*0.1/hbar*[-1.0,1.0];
   W_NUM = 600;
   n_s = loadStruct([outKeyQW,'system_structure.dat']);
   DBRind = [2.9459, 3.4345];
   
   % Approx width of DBR
   DF = w0*(4/pi)*asin((min(DBRind)-max(DBRind))/(sum(DBRind)));
   fa = w0 + DF/2;
   fb = w0 - DF/2;
   DBR_LIM = [fa,fb];
   
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
    e_cav_re = loadD([outKeyABS,num2str(n),'_E_re_',location,'_',point,'.dat']);
    e_cav_im = loadD([outKeyABS,num2str(n),'_E_im_',location,'_',point,'.dat']);

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
    
    % ABSORBER
    [w_abs,Y1_abs,Y2_abs] = getSpectrums(y0_abs,y1_abs,t0_abs(2)-t0_abs(1),t1_abs(2)-t1_abs(1));

    rat2 = Y2_abs./Y1_abs;

    %% CLEAN UP DATA
    [rat1,w_qw]  = cleanUp(w_qw,rat1,W_LIM);
    [rat2,w_abs] = cleanUp(w_abs,rat2,W_LIM);

    %% Reflection spectrums
    %WTF
    refSpecQW = (abs(rat1)-1)/compReflect(j).pass;
    refSpecABS = (1-abs(rat2))/compReflect(j).pass;
    
    if length(w_qw)>length(w_abs)
    w = w_abs;
    refSpecQW = interp1(w_qw,refSpecQW,w)';
    else
    w = w_qw;
    refSpecABS = interp1(w_abs,refSpecABS,w)';
    end

    refSpecABS = refSpecABS+compReflect(j).OC_loss/100; %Artificial constant addditional OC loss
    
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
    
    compReflect(j).refSpecQW=refSpecQW;
    compReflect(j).refSpecABS=refSpecABS;
    compReflect(j).w=w;
    compReflect(j).diffGain=diffGain;
    
    %% Dispersion
    abs1_deph = read_deph_scale([outKeyMat,'ABS1.config']);
    abs1_deph = abs1_deph/(2.0*Eb/hbar)/fs;

    qw1_deph = read_deph_scale([outKeyMat,'QW1.config']);
    qw1_deph = qw1_deph/(2.0*Eb/hbar)/fs;
    
    dy = yl(2);
    qw1_dens  = readDensity([outKeyMat,'QW1.config']);
    abs1_dens = readDensity([outKeyMat,'ABS1.config']);
    qw1_bg  = read_BG([outKeyMat,'QW1.config']);
    abs1_bg = read_BG([outKeyMat,'ABS1.config']);

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

    compReflect(j).gd_qw   = (phi(2:end) - phi(1:end-1))/(dw);
    compReflect(j).gdd_qw  = (phi(3:end) - 2*phi(2:end-1) + phi(1:end-2))/(dw^2);
    compReflect(j).gddd_qw = (0.5*phi(5:end) - phi(4:end-1) + phi(2:end-3) - 0.5*phi(1:end-4))/(dw^3);

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

    compReflect(j).gd_abs   = (phi(2:end) - phi(1:end-1))/(dw);
    compReflect(j).gdd_abs  = (phi(3:end) - 2*phi(2:end-1) + phi(1:end-2))/(dw^2);
    compReflect(j).gddd_abs = (0.5*phi(5:end) - phi(4:end-1) + phi(2:end-3) - 0.5*phi(1:end-4))/(dw^3);

    %Clean up picture based on FT
    compReflect(j).gd_x   = hbar*w(2:end  )/e;
    compReflect(j).gdd_x  = hbar*w(2:end-1)/e;
    compReflect(j).gddd_x = hbar*w(3:end-2)/e;

    %% Dispersion: SUM

    compReflect(j).gd   = compReflect(j).gd_qw + compReflect(j).gd_abs;
    %compReflect(j).gdd  = compReflect(j).gdd_qw + compReflect(j).gdd_abs;
    compReflect(j).gdd = compReflect(j).gdd_abs;
    compReflect(j).gddd = compReflect(j).gddd_qw + compReflect(j).gddd_abs;
    compReflect(j).gdd_x = compReflect(j).gdd_x*e/hbar;
    compReflect(j).gdd_x_w = 2.0*pi*c0./compReflect(j).gdd_x;
end

%% Data visualizations

if IO_reflectionPlot==1
    tmp_fig=figure(104);
    set(tmp_fig,'Name','Reflection spectrum comp');
    set(tmp_fig,'Position',POS);
    hold on
    plot(hbar*w0/e*[1,1],10*[min(compReflect(1).refSpecQW),max(compReflect(1).refSpecQW)]*100,'k--','DisplayName','\omega_0');
    %plot(compReflect(1).w*hbar/e,100*compReflect(1).refSpecABS,'DisplayName','Abs.'); %Assume constant absorption
    plot(compReflect(1).w*hbar/e,100*compReflect(j).refSpecABS,'DisplayName','Absorption'); %Assume constant absorption
    for j=1:num_sims
        plot(compReflect(j).w*hbar/e ,100*compReflect(j).refSpecQW, 'DisplayName',[char(compReflect(j).legend),' gain']);
    end
    hold off
    grid on
    lgd=legend('show');
    set(lgd,'Location','best');
    xlabel('Energy [eV]')
    ylabel('Reflection [%]')
    xlim(xl);
    ylim(yl);
    if IO_save==1
        saveas(tmp_fig,[saveKey,'reflectionSpectrumComp.png']);
    end
end

if IO_energyDispersionPlot==1
    tmp_fig=figure(105);
    set(tmp_fig,'Name','Energy dispersion plot');
    set(tmp_fig,'Position',POS);
    hold on
    plot(hbar*w0*[1,1]/e,[-1,1]*max(compReflect(1).gdd_qw/(fs*fs)),'k--', 'DisplayName','\omega_0')
    for j=1:num_sims
        plot(hbar*compReflect(j).gdd_x/e ,compReflect(j).gdd/(fs*fs)/compReflect(j).pass , 'DisplayName', char(compReflect(j).legend));
    end
    grid on
    lgd=legend('show');
    set(lgd,'Location','best');
    xlabel('Energy [eV]')
    ylabel('GDD [fs^2]')
    xlim(xl);%hbar*DBR_LIM/e)
    xticks([xl(1),1.204,xl(2)]);
    ylim(yl_gdd)
    hold off
    if IO_save==1
        saveas(tmp_fig,[saveKey,'energyDispersionComp.png']);
    end  
end

if IO_wavelengthDispersionPlot==1
    tmp_fig=figure(106);
    set(tmp_fig,'Name','Wavelength dispersion plot');
    set(tmp_fig,'Position',POS);
    hold on
    plot((2*pi*c0/w0)*[1,1]/nm,yl_gdd,'k--', 'DisplayName','\omega_0')
    for j=1:num_sims
        plot(compReflect(j).gdd_x_w/nm,compReflect(j).gdd/(fs*fs)/compReflect(j).pass, 'DisplayName', char(compReflect(j).legend));
    end
    grid on
    lgd=legend('show');
    set(lgd,'Location','best');
    xlabel('Wavelength [nm]')
    ylabel('GDD [fs^2]')
    ylim(yl_gdd)
    xlim([990,1060])
    hold off
    if IO_save==1
        saveas(tmp_fig,[saveKey,'wavelengthDispersionComp.png']);
    end  
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