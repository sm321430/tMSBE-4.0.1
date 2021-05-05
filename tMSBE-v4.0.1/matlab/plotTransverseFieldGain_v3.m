%function plotTransverseField()
clear all
close all

%% IOs for various features in code order
IO_gainLoad=1; %Load refSpecQW simulation field and compute GDD
IO_gainPower=0; %Compute power and energy for two peaks
IO_absorberLoad=0; %Load refSpecABS simulation field and compute GDD
IO_initialFields=0; %Plot field raw data [MW/cm^2]
IO_gainQWpeaks=0; %Plot field peaks [MW/cm^2]
IO_reportSpotSizes=0; %Compute/Report spot sizes for two pulses
IO_gainSurf=0; %Surface plot of gain
IO_absSurf=0; %Surface plot of absorption
IO_netGainSurf=0; %Surface plot of net gain
IO_tempPump=1; %Plot temperature and density pump profiles
IO_GDD=0; %GDD surface plots
IO_crossSection=0; %Cross section plots
IO_save=1; %Save all plots generated
IO_break=0; %Break after first positive IO_flag is completed
IO_iterateQW=1; %For loading lattice setup (not fully tested when off)

setupConstants
setupPlot

stride_points=[1;10]; %Stride for plotting [space,time]
XLIM = [1.16, 1.26]; %Frequency plotting width
YLIM = [-300, 300];  %Spatial plotting width
%TLIM = [-0.25,0.25]; %Initial field temporal plot width
ZLIM = [-3,3];     %Surface plot of gain limits
ZLIM_abs = [0.9,2.4]; %Surface plot of absorption limits
ZLIM_netGain = [-3,1.5]; %Surface plot of net gain limits

outKey_base = '../../tMSBE-VCAV36-transModelocking2/run/';
outKey=[outKey_base,'refSpecQW__'];
outKey_abs=[outKey_base,'refSpecABS__'];

%outKey = '../run/refSpecQW__';
location='CAVP'; %Field location for uploading and saving
date='102620';
test='VCAV36-transModelocking2-n2p1-noTrans-noSESAM';
test_folder='test';
saveKey_local='Fall2019-Spring2020/VCAV/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey = [saveKey_fold,'/',test,'-',location,'-'];
user_entry = input(['saveKey=',saveKey,'? (y to continue): '], 's');
if user_entry~='y'
    gafuggle
end
if exist(saveKey_fold,'dir')~=7
    mkdir(saveKey_fold);
end


OC_loss=1.0; %%Output coupling loss within the cavity, added to ABS
w0 = loadD([outKey,'w0.dat']);
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])
W_LIM = w0+e*0.1/hbar*[-1.0,1.0];
w_plot = linspace(W_LIM(1),W_LIM(2),1000);

transverse_grid_y = loadD([outKey,'transverse_grid_y.dat']);

NUM_TRANSVERSE = length(transverse_grid_y);
plot_num = 0;

t = loadD([outKey,num2str(plot_num),'_t.dat']);
Eout = zeros(NUM_TRANSVERSE,length(t));

gain = zeros(NUM_TRANSVERSE,length(w_plot));
loss = zeros(NUM_TRANSVERSE,length(w_plot));
gdd = zeros(NUM_TRANSVERSE,length(w_plot));
temperature_profile = zeros(NUM_TRANSVERSE, 1);
density_profile = zeros(NUM_TRANSVERSE, 1);

if IO_gainLoad==1
    disp('Loading gain QW');
    t_gain=t;
    for i = 0:(NUM_TRANSVERSE-1)
       if(mod(i,10)==0 || i==NUM_TRANSVERSE-1)
           disp(['Loading point: ',num2str(i)]); 
       end

       if IO_iterateQW==1
            tmp = loadD([outKey,'lattice_setup_QW6_T',num2str(1+i),'.dat']);
            density_profile(1+i) = tmp(1);
            temperature_profile(1+i) = tmp(2);
       end


       %% GAIN

       pulse_re = loadD([outKey,num2str(plot_num),'_E_re_',location,'_T',num2str(i),'.dat']);
       pulse_im = loadD([outKey,num2str(plot_num),'_E_im_',location,'_T',num2str(i),'.dat']);

       E_cav = (pulse_re + 1i*pulse_im).*exp(-1i*t_gain*w0);
       Eout(1+i,:) = E_cav;

       % Find peaks
       [pk,pk_i] = findpeaks(abs(E_cav),'sortstr','descend');
       pk(3:end) = [];
       pk_i(3:end) = [];
       if(i==floor(NUM_TRANSVERSE/2) && numel(pk_i)>0)
            pk_gain_center=pk_i(1); %Store largest peak at central point
       end
       pk_i = sort(pk_i);

       if (numel(pk_i) > 0 )
           sizeIl = 400*fs;
           sizeIr = 2000*fs;

           % Focus in on important part only
           % Input to QW
           x0 = t_gain-t_gain(pk_i(1));
           indp = x0>=-sizeIl;
           indm = x0<=sizeIr;
           ind0 = indp==indm;
           t0_qw = t_gain(ind0);
           y0_qw = E_cav(ind0);
           E0_qw(i+1,1:length(E_cav(ind0)))=E_cav(ind0);
       else
           disp('No peaks in GAIN found');
       end
       if (numel(pk_i) == 2 )
           sizeIl = 400*fs;
           sizeIr = 2000*fs;

           % Output from QW
           x1 = t_gain-t_gain(pk_i(2));
           indp = x1>=-sizeIl;
           indm = x1<=sizeIr;
           ind1 = indp==indm;
           t1_qw = t_gain(ind1);
           y1_qw = E_cav(ind1);
           E1_qw(i+1,1:length(E_cav(ind1)))=E_cav(ind1);

           [w_qw,Y1_qw,Y2_qw] = getSpectrums(y0_qw,y1_qw,t0_qw(2)-t0_qw(1),t1_qw(2)-t1_qw(1));

           rat1 = Y2_qw./Y1_qw;

           refSpecQW = (abs(rat1)-1);
           gain(1+i,:) = interp1(w_qw,refSpecQW,w_plot);
       
            
            %% Dispersion: QW
            if IO_GDD==1
                phi = atan2(imag(rat1),real(rat1)); %Phase of complex signal
                phi = unwrap(phi);    %Correct for 2pi phase jump
                %phi = smooth(phi);

                dw = w_qw(2)-w_qw(1);
                gdd_qw  = (phi(3:end) - 2*phi(2:end-1) + phi(1:end-2))/(dw^2);
                gdd_w  = w_qw(2:end-1);

                gdd(1+i,:) = interp1(gdd_w, gdd_qw, w_plot); 
            end
        else
            gain(1+i,:) = 0;
            gdd(1+i,:) = 0;
            disp('No two distinct peaks found');
        end
    end
end

if IO_gainPower==1
    disp('Plotting power curves');
    power=zeros(2,length(E0_qw(1,:)));
    energy=zeros(2,1);
    power(1,:)=trapz(abs(E0_qw));
    energy(1)=trapz(power(1,:)); 
    power(2,1:length(E1_qw(1,:)))=trapz(abs(E1_qw));
    energy(2)=trapz(power(2,:));
    
    tmp_fig=figure('name','Power plots');
    set(tmp_fig,'Position',POS);
    semilogy(power(1,:));
    hold on
    semilogy(power(2,:));
    semilogy(abs(power(2,:)-power(1,:)));
    xlabel('time [a.u.]');
    ylabel('power [a.u.]');
    legend('Peak 1','Peak 2','\Delta');
    if IO_save==1
        saveas(tmp_fig,[saveKey,'power.png']);
    end
end

if IO_absorberLoad==1
    disp('Loading absorbing QW');
    t_abs = loadD([outKey_abs,num2str(plot_num),'_t.dat']);
    Eout2 = zeros(NUM_TRANSVERSE,length(t_abs));

    for i = 0:(NUM_TRANSVERSE-1)
        if(mod(i,10)==0 || i==NUM_TRANSVERSE-1)
            disp(['Loading point: ',num2str(i)]); 
        end
        %% LOSS
        pulse_re = loadD([outKey_abs,num2str(plot_num),'_E_re_CAVP_T',num2str(i),'.dat']);
        pulse_im = loadD([outKey_abs,num2str(plot_num),'_E_im_CAVP_T',num2str(i),'.dat']);

        E_cav = (pulse_re + 1i*pulse_im).*exp(-1i*t_abs*w0);
        Eout2(1+i,:) = E_cav;

        % Find peaks
        [pk,pk_i] = findpeaks(abs(E_cav),'sortstr','descend');
        pk(3:end) = [];
        pk_i(3:end) = [];
        if(i==floor(NUM_TRANSVERSE/2))
            pk_abs_center=pk_i(1);
        end
        pk_i = sort(pk_i);

        if (numel(pk_i) > 0)

            sizeIl = 200*fs;
            sizeIr = 600*fs;

            % Focus in on important part only
            % Input to QW

            x0 = t_abs-t_abs(pk_i(1));
            indp = x0>=-sizeIl;
            indm = x0<=sizeIr;
            ind3 = indp==indm;
            t0_qw = t_abs(ind3);
            y0_qw = E_cav(ind3);
        else
            disp('No peaks in ABSORBER found');
        end
        if (numel(pk_i) == 2)
            
            sizeIl = 200*fs;
            sizeIr = 600*fs;

            % Output from QW
            x0 = t_abs-t_abs(pk_i(2));
            indp = x0>=-sizeIl;
            indm = x0<=sizeIr;
            ind3 = indp==indm;
            t1_qw = t(ind3);
            y1_qw = E_cav(ind3);

            [w_qw,Y1_qw,Y2_qw] = getSpectrums(y0_qw,y1_qw,t0_qw(2)-t0_qw(1),t1_qw(2)-t1_qw(1));

            rat1 = Y2_qw./Y1_qw;

            refSpecABS = (1-abs(rat1));
            loss(1+i,:) = interp1(w_qw,refSpecABS,w_plot)+OC_loss/100;
        else
            loss(1+i,:) = 0;
        end
    end
end   
    
if IO_initialFields==1
   disp('Plotting initial fields');
   if IO_gainLoad==1
        disp('Plotting gain field');
        tmp_fig=figure('name','Input gain field');
        set(tmp_fig,'Position',POS);
        if exist('pk_gain_center','var')
            t_cur=t_gain-t_gain(pk_gain_center);
        else
            t_cur=t_gain;
        end
        contourf(t_cur(1:stride_points(2):end)/ps,...
        transverse_grid_y(1:stride_points(1):end)/um,...
        (abs(Eout(1:stride_points(1):end,1:stride_points(2):end)*cm*cm/1e6)),'edgecolor','none');
        colorbar;
        xlabel('time [ps]');
        ylabel(' y [\mu m] ');
        if (exist('pk_gain_center','var') && exist('YLIM','var'))
            ylim(YLIM);
        end
        if exist('TLIM','var')
            xlim(TLIM);
        end
        if IO_save==1
            saveas(tmp_fig,[saveKey,'E0_gain.png']);
        end
   end
   if IO_absorberLoad==1
        disp('Plotting absorber field');
        tmp_fig=figure('name','Input absorber field');
        set(tmp_fig,'Position',POS);
        contourf(t_abs(1:stride_points(2):end),...
        transverse_grid_y(1:stride_points(1):end)/um,...
        abs(Eout2(1:stride_points(1):end,1:stride_points(2):end)*cm*cm/1e6),'edgecolor','none')
        colorbar
        ylim(YLIM);
        xlabel('time [ps]');
        ylabel(' y [\mu m] ');
        if IO_save==1
            saveas(tmp_fig,[saveKey,'E0_abs.png']);
        end
   end
   if IO_break==1
        disp('Breaking after plotting initial fields');
        asd 
   end
end

if IO_gainQWpeaks==1
    disp('Plotting QW peaks data');
    tmp_fig=figure('name','Input gain field');
    set(tmp_fig,'Position',POS);
    contourf(abs(E0_qw(1:stride_points(1):end,1:stride_points(2):end)*cm*cm/1e6),'edgecolor','none')
    colorbar
    if IO_save==1
        saveas(tmp_fig,[saveKey,'E0_qw.png']);
    end
    
    tmp_fig=figure('name','Reflected gain field');
    set(tmp_fig,'Position',POS);
    contourf(abs(E1_qw(1:stride_points(1):end,1:stride_points(2):end)*cm*cm/1e6),'edgecolor','none')
    colorbar
    if IO_save==1
        saveas(tmp_fig,[saveKey,'E1_qw.png']);
    end
    if IO_break==1
        disp('Breaking after plotting QW peaks');
        asd
    end
end
    
if IO_reportSpotSizes==1
    yy = abs(Eout)*um;

    [~,ind] = findpeaks(yy(NUM_TRANSVERSE/2,:),'sortstr','descend');
    ind = sort(ind(1:2));

    spot = findSpotSize(transverse_grid_y,yy(:,ind(1)));
    fwhm_y = findFWHM(transverse_grid_y,yy(:,ind(1)));

    ind_pulse = abs(t-t(ind(1)))<200*fs;
    fwhm_t = findFWHM(t,yy(NUM_TRANSVERSE/2,ind_pulse));

    disp(['BEFORE: spot = ',num2str(spot/um,'%.2f'),' [um], fwhm = ',num2str(fwhm_y/um,'%.2f'),' [um], ',num2str(fwhm_t/fs,'%.2f'),' [fs]'])

    spot = findSpotSize(transverse_grid_y,yy(:,ind(2)));
    fwhm_y = findFWHM(transverse_grid_y,yy(:,ind(2)));

    ind_pulse = abs(t-t(ind(2)))<200*fs;
    fwhm_t = findFWHM(t,yy(NUM_TRANSVERSE/2,ind_pulse));

    disp(['AFTER: spot = ',num2str(spot/um,'%.2f'),' [um], fwhm = ',num2str(fwhm_y/um,'%.2f'),' [um], ',num2str(fwhm_t/fs,'%.2f'),' [fs]'])
end

%% Plotting gain surface
if IO_gainSurf==1
    disp('Plotting gain as surface');
    tmp_fig=figure('Name', 'Gain total');
    set(tmp_fig,'Position',POS);
    hold on
    surf(hbar*w_plot/e,transverse_grid_y/um,100*gain,'edgecolor','none')
    hold off
    xlabel('Energy [eV]')
    ylabel('y [\mum]')
    zlabel('Amplification [%]')
    grid on
    set(gca,'Ydir','reverse')

    [X,Y] = meshgrid(hbar*w_plot/e,transverse_grid_y/um);
    hold on
    C = contourf(X, Y, 100*gain, [-0.00001 0.00001],'linewidth',1);
    hold off

    zlim(ZLIM)
    xlim(XLIM)
    ylim(YLIM)


    zl = zlim;
    caxis(zl)
    view(-90,90)
    colormap('jet')
    colorbar
    if IO_save==1
        saveas(tmp_fig,[saveKey,'GainSurf.png']);
    end
    if IO_break==1
        asd
    end
end

%% Surface plot of absorption in cavity
if IO_absSurf==1
    disp('Plotting absorption as surface');
    tmp_fig=figure('Name', 'Absorption/Loss total');
    set(tmp_fig,'Position',POS);
    hold on
    surf(hbar*w_plot/e,transverse_grid_y/um,100*loss,'edgecolor','none')
    hold off
    xlabel('Energy [eV]')
    ylabel('y [\mum]')
    zlabel('Loss [%]')
    grid on
    set(gca,'Ydir','reverse')

    [X,Y] = meshgrid(hbar*w_plot/e,transverse_grid_y/um);
    hold on
    C = contourf(X, Y, 100*loss, [OC_loss-0.00001 OC_loss+0.00001],'linewidth',1);
    hold off

    zlim(ZLIM_abs)
    xlim(XLIM)
    ylim(YLIM)


    zl = zlim;
    caxis(zl)
    view(-90,90)
    colormap('jet')
    colorbar
    if IO_save==1
        saveas(tmp_fig,[saveKey,'AbsSurf.png']);
    end
    if IO_break==1
        asd
    end
end

%% PLot net gain subtracting off losses
if IO_netGainSurf==1
    disp('Plotting net gain as surface');
    tmp_fig=figure('Name', 'Net gain');
    set(tmp_fig,'Position',POS);
    hold on
    surf(hbar*w_plot/e,transverse_grid_y/um,100*(gain-loss),'edgecolor','none')
    hold off
    xlabel('Energy [eV]')
    ylabel('y [\mum]')
    zlabel('Net gain [%]')
    grid on
    set(gca,'Ydir','reverse')

    [X,Y] = meshgrid(hbar*w_plot/e,transverse_grid_y/um);
    hold on
    C = contourf(X, Y, 100*(gain-loss), [-0.00001 0.00001],'linewidth',1);
    hold off

    zlim(ZLIM_netGain)
    %zlim([-10,32])
    xlim(XLIM)
    ylim(YLIM)


    zl = zlim;
    caxis(zl)
    view(-90,90)
    colormap('jet')
    colorbar
    if IO_save==1
        saveas(tmp_fig,[saveKey,'NetGainSurf.png']);
    end
    if IO_break==1
        asd
    end
end

    %% Temperature / pump profile
if IO_tempPump==1
    disp('Plotting temperature and pump profiles');
    h = figure('Name','Temperature/Pump profiles');
    set(h,'Position',POS_pump);
    [AX,H1,H2] = plotyy(transverse_grid_y/um, density_profile*1.0e-16, transverse_grid_y/um, temperature_profile);
    
    set(AX,{'ycolor'},{'b';'r'})
    set(H1,'color','b','LineStyle','-','marker','none')
    set(H2,'color','r','LineStyle','-','marker','none')

    grid on
    xlabel('y [\mum]')

    axes(AX(1))
    ylabel('Density [1/m^2 10^{16}]')
    xlim(YLIM)
    set(AX(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')

    yl = ylim;
    
%     bg_SG_degree = 16;
%     bg_SG_FWHM = 0.78*(1050*um);
%     bg_sigma = bg_SG_FWHM/(2*sqrt(2)*log(2)^(1/(2*bg_SG_degree)));
%     bg_shape = exp(-(transverse_grid_y/(sqrt(2)*bg_sigma)).^(2*bg_SG_degree));
%     hold on 
%     plot(transverse_grid_y/um,yl(2)*bg_shape,'k--')
%     hold off
    
    

    axes(AX(2))
    ylabel('Temperature [K]')
    xlim(YLIM)
    set(AX(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
    
    if IO_save==1
        saveas(h,[saveKey,'TempPumpProfiles.png']);
    end
    if IO_break==1
        asd
    end
end

    %% GDD
if IO_GDD==1  
    disp('Plotting GDD');
    % Trim gdd data
    [cx,cy] = convert_contour_to_polygon(C);
    ir = cx < XLIM(1);
    cx(ir) = []; cy(ir) = [];
    ir = cx > XLIM(2);
    cx(ir) = []; cy(ir) = [];
    ir = cy < YLIM(1);
    cx(ir) = []; cy(ir) = [];
    ir = cy > YLIM(2);
    cx(ir) = []; cy(ir) = [];

    
    ind_in = inpolygon(X(:),Y(:),cx,cy);
    gdd(~ind_in) = nan;
    
    
    % Dispersion
    tmp_fig=figure('Name','GDD surf');
    set(tmp_fig,'Position',POS);
    surf(X,Y,gdd/(fs*fs),'edgecolor','none')
    xlabel('Energy [eV]')
    ylabel('y [\mum]')
    zlabel('GDD [fs^2]')
    grid on
    set(gca,'Ydir','reverse')
    
    %zlim([-10,50])
    xlim(XLIM)
    ylim(YLIM)
    
    %zl = zlim;
    %caxis(zl)
    view(-90,90)
    colormap('jet')
    colorbar
    if IO_save==1
        saveas(tmp_fig,[saveKey,'GDDSurf.png']);
    end
    if IO_break==1
        asd
    end
end     
    %% Cross-section gain y = 0
if IO_crossSection==1
    disp('Plotting cross section');
    gain_y0 = 100*gain(NUM_TRANSVERSE/2,:);
    abs_y0 = 100*loss(NUM_TRANSVERSE/2,:);
    
    net_gain_y0 = gain_y0-abs_y0;
    ind_pos = net_gain_y0 > 0;
    ind_rem = (hbar*w_plot/e<XLIM(1)).*(hbar*w_plot/e > XLIM(2)) > 0;
    ind_pos(ind_rem) = 0;
    net_gain_wlim = w_plot(ind_pos);
    net_gain_width = max(net_gain_wlim) - min(net_gain_wlim);
    
    
    tmp_fig=figure('Name', 'Center gain/abs');
    set(tmp_fig,'Position',POS);
    plot(hbar*w_plot/e, gain_y0,'b-',...
         hbar*w_plot/e, abs_y0,'k-')
    grid on
    xlim(XLIM)
    xlabel('Energy [eV]')
    ylabel('Amplification [%]')
    title('Crossection at y = 0')
    %ylim([-10,10])
    
    yl = ylim;
    xl = xlim;
    text(xl(1) + 0.01*(xl(2)-xl(1)), yl(1) + 0.8*(yl(2)-yl(1)),['Net gain width = ',num2str(1000*hbar*net_gain_width/e,'%.1f'),' [meV]'])
    
    hold on
    yl =ylim;
    plot(hbar*w0*[1,1]/e, yl,'k--')
    hold off
    if IO_save==1
        saveas(tmp_fig,[saveKey,'CentralRefl.png']);
    end
    
    tmp_fig=figure('Name', 'Center GDD');
    set(tmp_fig,'Position',POS);
    plot(hbar*w_plot/e, gdd(NUM_TRANSVERSE/2,:)/(fs*fs),'b-')
    grid on
    xlim(XLIM)
    xlabel('Energy [eV]')
    ylabel('GDD [fs^2]')
    %ylim([-100,100])
    
    
    hold on
    yl =ylim;
    plot(hbar*w0*[1,1]/e, yl,'k--')
    hold off
    if IO_save==1
        saveas(tmp_fig,[saveKey,'CentralGDD.png']);
    end
    if IO_break==1
        asd
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

function [Y,X] = cleanUp(x, y, XL)

x0 = XL(1);
x1 = XL(2);

ind_m = x <= x0;
ind_p = x >= x1;
ind = ind_m==ind_p;

X = x(ind);
Y = y(ind);

end

function v = loadD(name)
% Read a single double or a list of doubles into v
% v is a row vector

fid = fopen(name,'rb');
v = fread(fid,'double');
fclose(fid);

end

function [x,y] = convert_contour_to_polygon(C)
%% Transfer contour matrix C to [x,y] pairs

[n,m] = size(C);

ind_split = [];


i= 1;
while i <= m
    
    skip = C(2,i);
    ind_split = [ind_split, i];
    
    i = i + skip +1;
    
end

C(:,ind_split) = [];
x = C(1,:);
y = C(2,:);


end

function fwhm = findSpotSize(DT,pulse)
% Find FWHM of central pulse. (For interfering pulses)
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


fwhm = abs(T1 - T0)/2;
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


