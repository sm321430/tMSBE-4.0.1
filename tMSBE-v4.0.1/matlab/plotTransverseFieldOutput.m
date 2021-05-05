function plotTransverseFieldOutput()


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


outKey = '../run/out__';

w0 = loadD([outKey,'w0.dat']);
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])


round_trip_time = loadD([outKey,'round_trip_time.dat']);
transverse_grid_y = loadD([outKey,'transverse_grid_y.dat']);
transverse_grid_device_y = loadD([outKey,'transverse_grid_device_y.dat']);
[~, ind_device_y] = intersect(transverse_grid_y, transverse_grid_device_y);

NUM_TRANSVERSE = length(dir([outKey,num2str(0),'_E_re_OUTPUT_T*.dat']))
NUM_TRANSVERSE_DEVICE = length(dir([outKey,num2str(0),'_E_re_QW6_T*.dat']))
plot_num = 19;
FOCUS = 10; % Compensate for lens in output: file output is from SESAM side, which has a spot (w/focus) where w is on the GAIN chip spot




%% Plot final output
t = loadD([outKey,num2str(plot_num-1),'_t.dat']);
out_pulse = zeros(NUM_TRANSVERSE,length(t));
for i = 0:(NUM_TRANSVERSE-1)

    pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_OUTPUT_T',num2str(i),'.dat']);
    pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_OUTPUT_T',num2str(i),'.dat']);
    pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*w0);

    out_pulse(1+i,:) = 0.5*eps0*c0*abs(pulse).^2;
end

figure
if (NUM_TRANSVERSE > 1)
    surf((t-t(1))/ps,transverse_grid_y/um,out_pulse*cm*cm/1e6,'edgecolor','none')
    
    ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
    view(67,16)
    ylabel('y [\mum]')
    zlabel('I [MW/cm^2]')
else
    plot((t-t(1))/ps, out_pulse*cm*cm/1e6, 'b-')
    ylabel('I [MW/cm^2]')
end

grid on
xlabel('t [ps]')
xlim([min(t-t(1)), max(t-t(1))]/ps)



%% Plot QW density and temperature
qw_density = zeros(NUM_TRANSVERSE_DEVICE,length(t));
qw_temp    = zeros(NUM_TRANSVERSE_DEVICE,length(t));
abs_density = zeros(NUM_TRANSVERSE_DEVICE,length(t));
abs_temp    = zeros(NUM_TRANSVERSE_DEVICE,length(t));
for i = 0:(NUM_TRANSVERSE_DEVICE-1)
    
    Ne = loadD([outKey,num2str(plot_num-1),'_Nsum_e_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
    Te = loadD([outKey,num2str(plot_num-1),'_inst_temp_e_QW6_T',num2str(ind_device_y(1+i)),'.dat']);

    qw_density(1+i,:) = Ne;
    qw_temp(1+i,:) = Te;
    
    Ne = loadD([outKey,num2str(plot_num-1),'_Nsum_e_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
    Te = loadD([outKey,num2str(plot_num-1),'_inst_temp_e_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);

    abs_density(1+i,:) = Ne;
    abs_temp(1+i,:) = Te;
end

figure
subplot(2,1,1)
plot(transverse_grid_device_y/um,min(qw_density,[],2),'bo-',...
     transverse_grid_device_y/um,max(qw_density,[],2),'ko-')
 grid on
 xlabel('y [um]')
 ylabel('Density [-]')
 legend('min','max')
 
subplot(2,1,2)
plot(transverse_grid_device_y/um,min(qw_temp,[],2),'bo-',...
     transverse_grid_device_y/um,max(qw_temp,[],2),'ko-')
 grid on
 xlabel('y [um]')
 ylabel('Temperature [-]')
 legend('min','max')
 
 
figure
subplot(2,1,1)
plot(transverse_grid_device_y/um,min(abs_density,[],2),'bo-',...
     transverse_grid_device_y/um,max(abs_density,[],2),'ko-')
 grid on
 xlabel('y [um]')
 ylabel('Density [-]')
 legend('min','max')
 
subplot(2,1,2)
plot(transverse_grid_device_y/um,min(abs_temp,[],2),'bo-',...
     transverse_grid_device_y/um,max(abs_temp,[],2),'ko-')
 grid on
 xlabel('y [um]')
 ylabel('Temperature [-]')
 legend('min','max')

 


for i = 1:length(t)
    
    abs_density(:,i) = abs_density(:,i) - abs_density(:,end);
    abs_temp(:,i) = abs_temp(:,i) - abs_temp(:,end);
   
    qw_density(:,i) = qw_density(:,i) - qw_density(:,end);
    qw_temp(:,i) = qw_temp(:,i) - qw_temp(:,end);
end



figure
subplot(2,1,1)
if (NUM_TRANSVERSE>1)
    surf((t-t(1))/ps,transverse_grid_device_y/um,qw_density/1e16,'edgecolor','none')
    ylabel('y [um]')
else
	plot((t-t(1))/ps, qw_density/1e16)
    ylabel('N [10^{16}m^{-2}]')
end
xlabel('t [ps]')
title('Density QW1 [1e16]')
colorbar
grid on

subplot(2,1,2)
if (NUM_TRANSVERSE>1)
    surf((t-t(1))/ps,transverse_grid_device_y/um,qw_temp,'edgecolor','none')
    ylabel('y [um]')
else
	plot((t-t(1))/ps, qw_temp)
    ylabel('T [K]')
end
xlabel('t [ps]')
title('Temperature QW1 [K]')
colorbar
grid on

figure
subplot(2,1,1)
if (NUM_TRANSVERSE>1)
    surf((t-t(1))/ps,transverse_grid_device_y/um,abs_density/1e14,'edgecolor','none')
    ylabel('y [um]')
else
	plot((t-t(1))/ps, abs_density/1e14)
    ylabel('N [10^{14}m^{-2}]')
end
xlabel('t [ps]')
title('Density ABS1 [1e14]')
colorbar
grid on

subplot(2,1,2)
if (NUM_TRANSVERSE>1)
    surf((t-t(1))/ps,transverse_grid_device_y/um,abs_temp,'edgecolor','none')
    ylabel('y [um]')
else
	plot((t-t(1))/ps, abs_temp)
    ylabel('T [K]')
end
xlabel('t [ps]')
title('Temperature ABS1 [K]')
colorbar
grid on



%% Density and temperature
temperature_profile = zeros(NUM_TRANSVERSE, 1);
density_profile = zeros(NUM_TRANSVERSE, 1);
for i = 0:(NUM_TRANSVERSE-1)
    
    tmp = loadD([outKey,'lattice_setup_QW6_T',num2str(1+i),'.dat']);
    density_profile(1+i) = tmp(1);
    temperature_profile(1+i) = tmp(2);
end


figure
[AX,H1,H2] = plotyy(transverse_grid_y/um, density_profile/1e16,transverse_grid_y/um, temperature_profile);
set(AX,{'ycolor'},{'b';'r'})
set(H1,'color','b','LineStyle','-','marker','o')
set(H2,'color','r','LineStyle','-','marker','o')

grid on
xlabel('y [\mum]')

axes(AX(1))
ylim([0, 2.5])
ylabel('Density [10^{16} m^{-2}]')
set(AX(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')

axes(AX(2))
ylim([300, 390])
ylabel('Temperature[K]')
set(AX(2), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')




%% Peak intensity and FWHM (time and transverse dim)
max_I = zeros(NUM_TRANSVERSE, plot_num);
max_It = zeros(1,plot_num);
fwhm_t = zeros(1,plot_num);
fwhm_y = zeros(1,plot_num);
for j = 0:(plot_num-1)
    
    t = loadD([outKey,num2str(j),'_t.dat']);
    
    trans_profile = zeros(NUM_TRANSVERSE,length(t));
    for i = 0:(NUM_TRANSVERSE-1)
    
        pulse_re = loadD([outKey,num2str(j),'_E_re_OUTPUT_T',num2str(i),'.dat']);
        pulse_im = loadD([outKey,num2str(j),'_E_im_OUTPUT_T',num2str(i),'.dat']);
        pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*w0);
        
        max_I(1+i,1+j) = max(0.5*eps0*c0*abs(pulse).^2);
        
        trans_profile(1+i,:) = pulse;
    end
        
    % Find a pulse to focus on
    if (NUM_TRANSVERSE>1)
        [~,ind] = findpeaks(abs(trans_profile(NUM_TRANSVERSE/2,:)).^2,'sortstr','descend');
        
        ind = sort(ind(1:2));


        fwhm_y(1+j) = findFWHM(transverse_grid_y*FOCUS, abs(trans_profile(:,ind(1))) );
        %fwhm_y(1+j) = findSpotSize(transverse_grid_y,abs(trans_profile(:,ind(1))));
        fwhm_t(1+j) = findFWHM(t, abs(trans_profile(NUM_TRANSVERSE/2,:)) );
        
    else
        fwhm_t(1+j) = findFWHM(t, abs(trans_profile) );
    end
    
    
    
    max_It(1+j) = mean(t);
end

yy = max_I*cm*cm/1e6;
%yy = log(yy/max(yy(:)));

figure
subplot(3,1,1)
if (NUM_TRANSVERSE>1)
    surf(max_It/ns,transverse_grid_y*FOCUS/um,yy,'edgecolor','none')
    ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
    ylabel('y [\mum]')
    zlabel('I(t) [MW/cm^2]')
    view(25,22)
else
    plot(max_It/ns, yy)
    ylabel('I(t) [MW/cm^2]')
end

zl = zlim;
hold on
yy2 = (density_profile-min(density_profile))/(max(density_profile)-min(density_profile));
plot3(0*ones(size(transverse_grid_y))/ps, transverse_grid_y/um, zl(1) + yy2*(zl(2)-zl(1)),'ko-')
hold off

xlabel('t [ns]')

grid on



subplot(3,1,2)

if (NUM_TRANSVERSE>1)
    plot(max_It/ns,yy(NUM_TRANSVERSE/2,:),'bo-')
else
    plot(max_It/ns, yy)
end
grid on
xlabel('t [ns]')
ylabel('I(t) [MW/cm^2]')



subplot(3,1,3)
[AX,H1,H2] = plotyy(max_It/ns,fwhm_t/fs, max_It/ns,fwhm_y/um);
set(AX,{'ycolor'},{'b';'r'})
set(H1,'color','b','LineStyle','-','marker','o')
set(H2,'color','r','LineStyle','-','marker','o')

grid on
xlabel('t [ns]')

axes(AX(1))
ylabel('FWHM [fs]')
set(AX(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')

axes(AX(2))
ylabel('FWHM [\mum]')
set(AX(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')




end


function v = loadD(name)
% Read a single double or a list of doubles into v
% v is a row vector

fid = fopen(name,'rb');
v = fread(fid,'double');
fclose(fid);

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

function [k,w,y] = getFFT2(y,dx,dt)
%% Find fourier transform of y(NX,NY) sampled with timestep (dx,dt)
% PAD offers options to zeropad the data for increased resolution and to
% compare with other signals of different length
% 

[n,m] = size(y);
NFFT_t = 2^nextpow2(m);
NFFT_x = 2^nextpow2(n);

L = length(y(:));

y = fftshift(fft2(y,NFFT_x,NFFT_t)/L);
k = (1.0/(2*dx))*linspace(-1,1,NFFT_x);
f = (1.0/(2*dt))*linspace(0,1,NFFT_t/2+1);
y = sqrt(2)*y(:,1:NFFT_t/2+1);


y = y(:, end:-1:1);
w = 2*pi*f; % Go to omega and reverse frequency axis



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

X1 = X1*Xscale;

end

function spot_size = findSpotSize(DT,pulse)
% Find spot of pulse w in exp(-(r/w)^2) 
% Input is a pulse with time difference DT. That is... DT is refrenced from
% 0

[maxp, maxi] = max(pulse); % max pulse

% Find left boundary
i0 = maxi;
while ((i0 > 2)&&(pulse(i0-1) > exp(-1)*maxp))
    
    i0 = i0 -1;
end

% Find right boundary
i1 = maxi;
while ((i1 < length(pulse))&&(pulse(i1+1) > exp(-1)*maxp))
    
    i1 = i1 +1;
end
%fwhm = DT(i1) - DT(i0);


% Interpolate data
ind0 = i0-1:i0+1;
T0 = invertYData(DT(ind0),pulse(ind0),exp(-1)*maxp);

ind1 = i1-1:i1+1;
T1 = invertYData(DT(ind1),pulse(ind1),exp(-1)*maxp);


spot_size = abs(T1 - T0)/2;
end


