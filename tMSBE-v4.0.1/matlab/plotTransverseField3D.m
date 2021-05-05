function plotTransverseField3D()


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
NUM_TRANSVERSE_DEVICE = length(dir([outKey,num2str(0),'_E_re_QW1_T*.dat']))
plot_num = 16;




%% Plot final output
t = loadD([outKey,num2str(plot_num-1),'_t.dat']);
out_pulse = zeros(NUM_TRANSVERSE,length(t));
for i = 0:(NUM_TRANSVERSE-1)

    pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_OUTPUT_T',num2str(i),'.dat']);
    pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_OUTPUT_T',num2str(i),'.dat']);
    pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*w0);

    out_pulse(1+i,:) = 0.5*eps0*c0*abs(pulse).^2;
end

%{
figure
hold on
surf((t-t(1))/ps,transverse_grid_y/um,log(out_pulse*cm*cm/1e6),'edgecolor','none')
hold off
grid on
xlabel('t [ps]')
ylabel('y [\mum]')
zlabel('I [MW/cm^2]')
ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
xlim([min(t-t(1)), max(t-t(1))]/ps)
view(67,16)
%}

%% Isolate pulse
DT = 1800*fs;


[~,ind] = findpeaks(out_pulse(NUM_TRANSVERSE/2,:),'SortStr','descend');
PLOT_R = linspace(min(transverse_grid_y),max(transverse_grid_y),200);

ind_t = abs(t-t(ind(1)))<DT;
t_new = t(ind_t) - t(ind(1));
out_pulse = out_pulse(:,ind_t);

figure
hold on
surf(t_new/ps,transverse_grid_y/um,out_pulse*cm*cm/1e6,'edgecolor','none')
hold off
grid on
xlabel('t [ps]')
ylabel('y [\mum]')
zlabel('I [MW/cm^2]')
ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
%xlim([min(t-t(1)), max(t-t(1))]/ps)
view(67,16)

%% CUT in half
out_pulse = out_pulse((1+NUM_TRANSVERSE/2):end,:);
transverse_grid_y = transverse_grid_y((1+NUM_TRANSVERSE/2):end);

figure
hold on
surf(t_new/ps,transverse_grid_y/um,out_pulse*cm*cm/1e6,'edgecolor','none')
hold off
grid on
xlabel('t [ps]')
ylabel('y [\mum]')
zlabel('I [MW/cm^2]')
ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
%xlim([min(t-t(1)), max(t-t(1))]/ps)
view(67,16)

%% Rotate one slice
out_pulse = out_pulse/max(out_pulse(:));


plot_index = 1:10:length(t_new);

[T,X,Y] = meshgrid(t_new(plot_index)/ps, PLOT_R'/um, PLOT_R'/um);

rho = sqrt(X.^2 + Y.^2);

Z_tot = 0*X;

for i = 1:length(plot_index)
    
    i
    plot_slice = out_pulse(:,plot_index(i));
    Z_tot(:,i,:) = interp1(transverse_grid_y/um, plot_slice,rho(:,i,:));

end

% Create cuts
cutoff = -12.8;
Z_tot = log(Z_tot);

%cutoff = 0.01;

ind_cut = 1:(length(PLOT_R)/2); % Create cut through data

% Create plane cuts through data
Z_tot(ind_cut,:,:) = []; % Cut with plane



X(ind_cut,:,:) = [];
Y(ind_cut,:,:) = [];
T(ind_cut,:,:) = [];

size(X)
size(Y)
size(T)
size(Z_tot)


%% Surface two
figure
p1 = patch(isosurface(T,X,Y,Z_tot,cutoff),'FaceColor','blue','EdgeColor','none');
p2 = patch(isocaps(T,X,Y,Z_tot,cutoff),'FaceColor','interp','EdgeColor','none');



view(24,2)

axis tight
%daspect([1,2,1])
%colormap(gray(100))
camlight(gca,24,12)
camlight
lighting gouraud
isonormals(T,X,Y, Z_tot,p1)
%isonormals(T,X,Y, Z_tot,p2)

%xlabel('t [ps]')
%ylabel('y [\mum]')
%zlabel('x [\mum]')
%grid on


set(gca,'visible','off')


asd



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


