%% Plot transverse field outputs for tMSBE
%General visualization for tMSBE with non-normal incidence capabilities
%June 12, 2020 S.A.Mclaren
%June 26, 2020 - Added additional IOs for superior control (SAM)

%General output functions for tMSBE
clear all
close all

%IOs for each plot located here
IO_convergence=0; %Plot convergence of FWHM, peak intensity, and energy
IO_temp_dens=0; %Plot density and temperature (temp not currently configured)
IO_densityVideo=0; %Video of density evolution (Requires IO_temp_dens flag)
IO_spectrums=0; %Get and plot spectrums for individual pulses
IO_finalOutput=1; %Plot final output
IO_finalOutput_indPulse=0; %Individual plots of pulses at final output
IO_iterateQW=1; %Include QWs if appropriate
IO_power=1; %Plot power overlay
IO_transIntensity=0; %Plot transverse intensity profiles
IO_pointDiff=0; %Plot time differences between peaks
IO_findFWHM=0; %Find FWHM for each peak (saved in pulseFWHM)
IO_fourFields=0; %Contour plots for four fields individual pulses (generates many plots)
IO_save=1; %Save figures
IO_separatePeaks=1; %Ensure peaks are distinct by set amount
%IO_twoArmQW=1;

%% Preliminaries
setupConstants;
setupPlot;

% global um;
% global ps;
% 
% fs = 1.0e-15;
% ps = 1.0e-12;
% um = 1.0e-6;
% ns = 1.0e-9;
% cm = 1.0e-2;
% nm = 1.0e-9;
% 
% hbar = 1.054589e-34;
% e = 1.602189e-19;
% c0   = 2.99792458E+08;
% mu0  = (4.0e-7)*pi;
% eps0 = 1.0/(mu0*c0*c0);
% n=3.4453; %Background refractive index (change as needed)
% 
% set(0,'defaulttextinterpreter','tex') %Default
% %set(0,'defaulttextinterpreter','latex')
% set(0,'defaultAxesFontName', 'Arial')
% set(0,'defaultAxesFontsize', 30)
% set(0,'defaultAxesLinewidth', 5)
% set(0,'defaultTextFontsize', 30)
% set(0,'defaultlineMarkerSize',26)
% set(0,'defaultlinelinewidth',5) %Thin lines
% maps %Load maps file for color scheme
% POS=[1,1,1200,800];
% POS2=[1,1,500,900];
% %POS_plot=[0.14,0.1100,0.53125,0.85];
% %POS_plot=[0.1330,0.1100,0.5930,0.8150];
% my_lineStyle={'-kd','-bo','-r+', '-gx', '-c*','-mo','-k+','-kx','-k*','-b+','-bx','-b*','-g+'};

date='050321';
%outKey ='/Volumes/SAMbackup/tMSBE-dualProp-data-2021/tMSBE-v4.0.1-DualPropCavity2-n3p0/run/out__';
outKey='/Volumes/SAMbackup/tMSBE-dualProp-data-2021/tMSBE-v4.2-KerrLensSML1_untested_DualProp_sepChips/run/out__';
test='tMSBE-dualProp-sepChips1';
%outKey = '../run/out__';
location='OUTPUTBACK'; %Field location for uploading and saving
OC=1;

%Automatic plot discovery
plot_num=1;
 while isfile([char(outKey),num2str(plot_num),'_E_re_',char('OUTPUT'),'_T0.dat'])
       plot_num=plot_num+1; %Check for next plot. If exists, add to counter
 end
 
plot_num=1; %Output number +1


test_folder='test';
saveKey_local='Fall2020-Summer2021/dualProp/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey = [saveKey_fold,'/',test,'-out',num2str(plot_num-1),'-',location,'-'];
user_entry = input(['saveKey=',saveKey,'? (y to continue): '], 's');
if user_entry~='y'
    gafuggle
end
if exist(saveKey_fold,'dir')~=7
    mkdir(saveKey_fold);
end

w0 = loadD([outKey,'w0.dat']);
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])


FOCUS = 1; % Compensate for lens in output: file output is from SESAM side, which has a spot (w/focus) where w is on the GAIN chip spot
round_trip_time = loadD([outKey,'round_trip_time.dat']);
transverse_grid_y = loadD([outKey,'transverse_grid_y.dat']);
transverse_grid_y=FOCUS*transverse_grid_y; %Scale grid based on focus.

t = loadD([outKey,num2str(plot_num-1),'_t.dat']);

%Exclude edges in time
init_ind=1;%find((t-t(1))/ps>=13.8,1);
final_ind=length(t);%find((t-t(1))/ps>=15,1);
t=t(init_ind:final_ind);

%Loading and plotting parameters
num_peaks=3; %Number of peaks for computing peaks
num_peaks_plot=num_peaks;
%tmp_width_pulse=5.0;
PulseTLIM=[-0.2,0.2]; %Output plot temporal width (ps)
PulseTTICK=[-0.2,0,0.2];
%tmp_width_pulse_log=1.5; %Output plot temporal width (ps) for log plot
spc_width=50; %Spacial width of plots (um)
spc_width_density=50; %Spacial width for density plots (um)
spectrum_width=[1.17,1.24]; %Spectrum widths (eV)
PulseLogTLIM=[-5,15.0];
PulseLogTTICK=[-5,0,5,10,15.0]; %Log Pulse Plot temporal tick marks
PulseLogXTICK=[-spc_width,-floor(spc_width/2),0,floor(spc_width/2), spc_width]; %Log Pulse spatial tick marks
PowerTLIM=[-2.0,2.0]; %Power plot  temporal width
PowerTTICK=[-0.25,0.0,0.5]; %Power plot tick marks
tmp_width_power=3;%0.5*round_trip_time/(ps); %Energy integration parameter

NUM_TRANSVERSE = max([length(dir([outKey,num2str(plot_num-1),'_E_re_OUTPUT_T*.dat'])),...
    length(dir([outKey,num2str(plot_num-1),'_E_re_OUTPUTBACK_T*.dat'])),length(dir([outKey,num2str(plot_num-1),'_E_fp_re_QW6_T*.dat']))]);
if IO_iterateQW==1
    if NUM_TRANSVERSE>1
        NUM_TRANSVERSE_DEVICE = max(length(dir([outKey,num2str(plot_num-1),'_E_fp_re_QW6_T*.dat'])),length(dir([outKey,num2str(plot_num-1),'_E_fp_re_QW6_T*.dat'])));
        %NUM_TRANSVERSE_DEVICES_SHORT=length(find(abs(transverse_grid_device_y)<spc_width*um+dx));
        NUM_TRANSVERSE=max(NUM_TRANSVERSE,NUM_TRANSVERSE_DEVICE);
    else
        NUM_TRANSVERSE_DEVICE = 1;
        %NUM_TRANSVERSE_DEVICES_SHORT = 1;
    end
end

counter=0;

if IO_iterateQW==1
    transverse_grid_device_y = loadD([outKey,'transverse_grid_device_y.dat']);
    dx=0.9*abs(transverse_grid_device_y(3)-transverse_grid_device_y(2));
    [~, ind_device_y] = intersect(transverse_grid_y, transverse_grid_device_y);
    NUM_TRANSVERSE_DEVICE=length(ind_device_y);
else
    if NUM_TRANSVERSE>1
        dx=transverse_grid_y(2)-transverse_grid_y(1);
    else
        dx=0;
    end
end

pulse_tmp=zeros(1,length(t));
out_pulse=zeros(1,length(t));
out_pulse_fp=zeros(1,length(t));
out_pulse_fm=zeros(1,length(t));
out_pulse_bp=zeros(1,length(t));
out_pulse_bm=zeros(1,length(t));

switch location
    case 'QW6'
    for i = 0:(NUM_TRANSVERSE_DEVICE-1)
        i
        pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
        pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
        pulse = (pulse_re(init_ind:final_ind) + 1i*pulse_im(init_ind:final_ind)).*exp(-1i*t*w0);
        if abs(transverse_grid_device_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_device_y(1+i);
            counter=counter+1;
        end
        location_short='QW6';
    end
    case 'QW6-twoArm'
    for i = 0:(NUM_TRANSVERSE_DEVICE-1)
            i
            pulse_fp_re = loadD([outKey,num2str(plot_num-1),'_E_fp_re_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_fp_im = loadD([outKey,num2str(plot_num-1),'_E_fp_im_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_fp(i+1,:) = (pulse_fp_re(init_ind:final_ind) + 1i*pulse_fp_im(init_ind:final_ind)).*exp(-1i*t*w0);
            pulse_bp_re = loadD([outKey,num2str(plot_num-1),'_E_bp_re_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_bp_im = loadD([outKey,num2str(plot_num-1),'_E_bp_im_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_bp(i+1,:) = (pulse_bp_re(init_ind:final_ind) + 1i*pulse_bp_im(init_ind:final_ind)).*exp(-1i*t*w0);
            pulse_fm_re = loadD([outKey,num2str(plot_num-1),'_E_fm_re_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_fm_im = loadD([outKey,num2str(plot_num-1),'_E_fm_im_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_fm(i+1,:) = (pulse_fm_re(init_ind:final_ind) + 1i*pulse_fm_im(init_ind:final_ind)).*exp(-1i*t*w0);
            pulse_bm_re = loadD([outKey,num2str(plot_num-1),'_E_bm_re_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_bm_im = loadD([outKey,num2str(plot_num-1),'_E_bm_im_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_bm(i+1,:) = (pulse_bm_re(init_ind:final_ind) + 1i*pulse_bm_im(init_ind:final_ind)).*exp(-1i*t*w0);
            pulse=pulse_fp(i+1,:)+pulse_bp(i+1,:)+pulse_fm(i+1,:)+pulse_bm(i+1,:);
        if abs(transverse_grid_device_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse_fp(1+counter,:)=abs(pulse_fp(i+1,:));
            out_pulse_bp(1+counter,:)=abs(pulse_bp(i+1,:));
            out_pulse_fm(1+counter,:)=abs(pulse_fm(i+1,:));
            out_pulse_bm(1+counter,:)=abs(pulse_bm(i+1,:));
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_device_y(1+i);
            counter=counter+1;
        end
        location_short='QW6';
    end
    case 'ABS1'
    for i = 0:(NUM_TRANSVERSE_DEVICE-1)
        pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
        pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
        pulse = (pulse_re(init_ind:final_ind) + 1i*pulse_im(init_ind:final_ind)).*exp(-1i*t*w0);
        if abs(transverse_grid_device_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_device_y(1+i);
            counter=counter+1;
        end
        location_short='ABS1';
    end
    case 'ABS1-twoArm'
    for i = 0:(NUM_TRANSVERSE_DEVICE-1)
            pulse_fp_re = loadD([outKey,num2str(plot_num-1),'_E_fp_re_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_fp_im = loadD([outKey,num2str(plot_num-1),'_E_fp_im_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_fp = (pulse_fp_re(init_ind:final_ind) + 1i*pulse_fp_im(init_ind:final_ind)).*exp(-1i*t*w0);
            pulse_bp_re = loadD([outKey,num2str(plot_num-1),'_E_bp_re_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_bp_im = loadD([outKey,num2str(plot_num-1),'_E_bp_im_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_bp = (pulse_bp_re(init_ind:final_ind) + 1i*pulse_bp_im(init_ind:final_ind)).*exp(-1i*t*w0);
            pulse_fm_re = loadD([outKey,num2str(plot_num-1),'_E_fm_re_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_fm_im = loadD([outKey,num2str(plot_num-1),'_E_fm_im_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_fm = (pulse_fm_re(init_ind:final_ind) + 1i*pulse_fm_im(init_ind:final_ind)).*exp(-1i*t*w0);
            pulse_bm_re = loadD([outKey,num2str(plot_num-1),'_E_bm_re_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_bm_im = loadD([outKey,num2str(plot_num-1),'_E_bm_im_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
            pulse_bm = (pulse_bm_re(init_ind:final_ind) + 1i*pulse_bm_im(init_ind:final_ind)).*exp(-1i*t*w0);
            pulse=pulse_fp+pulse_fm+pulse_bp+pulse_bm;
        if abs(transverse_grid_device_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_device_y(1+i);
            counter=counter+1;
        end
        location_short='ABS1';
    end
    case 'OUTPUT'
    for i = 0:NUM_TRANSVERSE-1 %issue with last point im
        i
        pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(i),'.dat']);
        pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_',location,'_T',num2str(i),'.dat']);
        pulse = (pulse_re(init_ind:final_ind) + 1i*pulse_im(init_ind:final_ind)).*exp(-1i*t*w0);
        if abs(transverse_grid_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_y(1+i);
            counter=counter+1;
        end
    end
    case 'OUTPUTBACK'
    for i = 0:NUM_TRANSVERSE-1 %issue with last point im
        i
        pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(i),'.dat']);
        pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_',location,'_T',num2str(i),'.dat']);
        pulse = (pulse_re(init_ind:final_ind) + 1i*pulse_im(init_ind:final_ind)).*exp(-1i*t*w0);
        if abs(transverse_grid_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_y(1+i);
            counter=counter+1;
        end
    end
    case 'CAVP'
    for i = 0:NUM_TRANSVERSE-1 %issue with last point im
        i
        pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(i),'.dat']);
pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_',location,'_T',num2str(i),'.dat']);
        pulse = (pulse_re(init_ind:final_ind) + 1i*pulse_im(init_ind:final_ind)).*exp(-1i*t*w0);
        if abs(transverse_grid_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_y(1+i);
            counter=counter+1;
        end
    end
    otherwise
      try
        for i = 0:NUM_TRANSVERSE-1 %issue with last point im
            i
        pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(i),'.dat']);
pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_',location,'_T',num2str(i),'.dat']);
        pulse = (pulse_re(init_ind:final_ind) + 1i*pulse_im(init_ind:final_ind)).*exp(-1i*t*w0);
        if abs(transverse_grid_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_y(1+i);
            counter=counter+1;
        end
        end
        location_short=location;
      catch
        warning('Location is not recognized.');
        asd
      end
end

    NUM_TRANSVERSE_DEVICES_SHORT=counter;
    %transverse_grid_short_y=zeros(NUM_TRANSVERSE_DEVICES_SHORT,1);

%% Find peak_num distinct peaks
t_cur=(t-t(1))/ps; %Current time refernced from t(1) in picoseconds
ind_trans=ones(counter,20); %indices for peaks
peaks_trans=zeros(counter,50); %Values for peaks
points=zeros(num_peaks,1); %time values for peaksinitialized
tmp_pulse=out_pulse; %Temporary pulse variable for manipulation
% Find top peaks, away from edges for all transverse points
for ll=1:counter
    [peaks,time_ind]=findpeaks(tmp_pulse(ll,:),'SortStr','descend');
    ii=min(50,length(peaks));
    ind_trans(ll,1:ii)=time_ind(1:ii);
    peaks_trans(ll,1:ii)=peaks(1:ii);
end

% Out of found peaks, find max separated peaks
jk=1;
while jk<=num_peaks
    [row,col]=find(tmp_pulse==max(max(peaks_trans)),1);
    if IO_separatePeaks==1
        try
        if ((min(abs(points-t_cur(col)))>2 && (t_cur(end)-t_cur(col))>1.0)...
              ||((t_cur(col)<2 && t_cur(col)>0.25) && max(points(points<2))==0))
            points(jk)=t_cur(col);
            jk=jk+1;
            
        end
        catch
        end
    else
        points(jk)=t_cur(col);
        jk=jk+1; 
    end
    tmp_pulse(:,col)=0;
    [row,col]=find(peaks_trans==max(max(peaks_trans)),1);
    peaks_trans(row,col)=0;
end
points=sort(points);
points_diff=points(2:end)-points(1:end-1);

%% Plot point differences
if IO_pointDiff==1
   tmp_fig=figure(203);
   set(tmp_fig,'Name','Final Output Plot');
   set(tmp_fig,'Position',POS);
   plot(points(2:end),points_diff,'-ko');
   xlabel('time [ps]');
   ylabel('\Delta t [ps]');
   xlim([points(2),points(end)]);
   if IO_save==1
    saveas(tmp_fig,[saveKey,'finalOutputPointDifference.png']);   
   end
end

if IO_findFWHM==1
   disp('Finding FWHM');
   pulse_FWHM=zeros(num_peaks,1);
   for j=1:num_peaks
      ind=find(t/ps==points(j),1);
      pulse_FWHM(j)=findFWHM(transverse_grid_y,out_pulse(:,ind));
   end
   pulse_FWHM_ratio=(pulse_FWHM(2:end)-pulse_FWHM(1:end-1))./pulse_FWHM(1:end-1);
end

%% Plot spectrums
if IO_spectrums==1
    disp('Plotting spectrums');
    sizeIl=300*fs;
    sizeIr=2000*fs;
    if num_peaks==1
       points(2)=points(1); 
    end
    if points(1)*ps>sizeIl && (t(end)/ps-points(2)*ps)>sizeIr
       if NUM_TRANSVERSE>1
           for l=1:counter
               [w_qw, Y1_qw, Y2_qw] = getSpectrums_fromPoints(t_cur,sizeIl,sizeIr,points(1:2),squeeze(pulse_tmp(l,:)));
               pulse_spectrum1(l,1:length(Y1_qw))=Y1_qw;
               pulse_spectrum2(l,1:length(Y2_qw))=Y2_qw;
           end
       else
            [w_qw, Y1_qw, Y2_qw] = getSpectrums_fromPoints(t_cur,sizeIl,sizeIr,points(1:2),pulse_tmp);
            pulse_spectrum1=Y1_qw;
            pulse_spectrum2=Y2_qw;
       end  
       
       
       if NUM_TRANSVERSE>1
            pulse_spectrum1_integrated=zeros(length(Y1_qw),1);
            pulse_spectrum2_integrated=zeros(length(Y2_qw),1);
            for l=1:length(Y1_qw)
                pulse_spectrum1_integrated(l)=trapz(pulse_spectrum1(:,l));
            end
            for l=1:length(Y2_qw)
                pulse_spectrum2_integrated(l)=trapz(pulse_spectrum2(:,l));
            end
       else
           pulse_spectrum1_integrated=pulse_spectrum1;
           pulse_spectrum2_integrated=pulse_spectrum2;
       end
    else
        disp('IO_spectrums::pulses too close to end of window');
        asd
    end
    
    save([saveKey,'pulseSpectrum1.mat'],'pulse_spectrum1_integrated','w_qw');
    save([saveKey,'pulseSpectrum2.mat'],'pulse_spectrum2_integrated','w_qw');
    
    if NUM_TRANSVERSE>1
        tmp_fig=figure(707);
        set(tmp_fig,'Name','Pulse spectrum 1');
        contourf(hbar*w_qw/e, transverse_grid_short_y/um,abs(pulse_spectrum1),'edgecolor','none')
        ylim([-spc_width,spc_width]);
        xlim(spectrum_width);
        xlabel('Energy [eV]');
        ylabel('y [\mu m]');
        colorbar
        if IO_save==1
            saveas(tmp_fig,[saveKey,'pulseSpectrum1.png']);   
        end

        tmp_fig=figure(708);
        set(tmp_fig,'Name','Pulse spectrum 2');
        contourf(hbar*w_qw/e, transverse_grid_short_y/um,abs(pulse_spectrum2),'edgecolor','none')
        ylim([-spc_width,spc_width]);
        xlim(spectrum_width);
        xlabel('Energy [eV]');
        ylabel('y [\mu m]');
        colorbar

        if IO_save==1
            saveas(tmp_fig,[saveKey,'pulseSpectrum2.png']);   
        end
    else
        tmp_fig=figure(707);
        set(tmp_fig,'Name','Pulse spectrum 1');
        plot(hbar*w_qw/e, abs(pulse_spectrum1));
        %plot(2*pi*c0./w_qw*10^9,abs(pulse_spectrum1))
        grid on
        %xlim([800,1200]);
        xlim(spectrum_width);
        xlabel('Energy [eV]');
        ylabel('y [\mu m]');
        if IO_save==1
            saveas(tmp_fig,[saveKey,'pulseSpectrum1.png']);   
        end

        tmp_fig=figure(708);
        set(tmp_fig,'Name','Pulse spectrum 2');
        plot(hbar*w_qw/e, abs(pulse_spectrum2))
        grid on
        xlim(spectrum_width);
        xlabel('Energy [eV]');
        ylabel('y [\mu m]');
        if IO_save==1
            saveas(tmp_fig,[saveKey,'pulseSpectrum2.png']);   
        end
    end
end
%% Plot final output
if IO_finalOutput==1
disp('Plotting final output figures')

    %%Different plotting procedures for TMSBE
    if NUM_TRANSVERSE==1
        if (strcmp(location,'QW6-twoArm')||strcmp(location,'ABS1-twoArm'))
            tmp_fig=figure(11);
            set(tmp_fig,'Name','Final Output Plot');
            set(tmp_fig,'Position',POS);
            semilogy((t-t(1))/ps, (abs(pulse_fp)),(t-t(1))/ps, (abs(pulse_fm)),(t-t(1))/ps,...
                (abs(pulse_bp)),(t-t(1))/ps, (abs(pulse_bm)))
            legend('FP','FM','BP','BM');
            xlabel('time [ps]');
            ylabel('|E| [V/m]');
            xlim([min((t-t(1))/ps),max((t-t(1))/ps)]);
            if IO_save==1
                saveas(tmp_fig,[saveKey,'finalOutputPlotLog_fourFields.png']);   
            end

            tmp_fig=figure(110);
            set(tmp_fig,'Name','Final Output Plot');
            set(tmp_fig,'Position',POS);
            plot((t-t(1))/ps, (abs(pulse_fp)),(t-t(1))/ps, (abs(pulse_fm)),(t-t(1))/ps,...
                (abs(pulse_bp)),(t-t(1))/ps, (abs(pulse_bm)))
            legend('FP','FM','BP','BM');
            xlabel('time [ps]');
            ylabel('|E| [V/m]');
            xlim([min((t-t(1))/ps),max((t-t(1))/ps)]);
            if IO_save==1
                saveas(tmp_fig,[saveKey,'finalOutputPlot_fourFields.png']);  
            end
            
            tmp_fig=figure(111);
            set(tmp_fig,'Name','Final Output Plot (real)');
            set(tmp_fig,'Position',POS);
            plot((t-t(1))/ps, real(pulse_fp_re),(t-t(1))/ps,real (pulse_fm_re),(t-t(1))/ps,...
                real(pulse_bp_re),(t-t(1))/ps, real(pulse_bm_re))
            legend('FP','FM','BP','BM');
            xlabel('time [ps]');
            ylabel('real(E) [V/m]');
            xlim([min((t-t(1))/ps),max((t-t(1))/ps)]);
            if IO_save==1
                saveas(tmp_fig,[saveKey,'finalOutputPlotReal_fourFields.png']);  
            end
        end

        tmp_fig=figure(666);
        set(tmp_fig,'Name','Final Output Plot (log)');
        set(tmp_fig,'Position',POS);
        plot((t-t(1))/ps, log(out_pulse*cm*cm/1e6)/OC)
        xlabel('time [ps]');
        ylabel('I [MW/cm2]');
        grid on
        xlim([min((t-t(1))/ps),max((t-t(1))/ps)]);
        if IO_save==1
            saveas(tmp_fig,[saveKey,'finalOutputPlotLog.png']); 
        end

        tmp_fig=figure(22);
        set(tmp_fig,'Name','Final Output Plot');
        set(tmp_fig,'Position',POS);
        plot((t-t(1))/ps, (out_pulse*cm*cm/1e6)/OC)
        xlabel('time [ps]');
        ylabel('I [MW/cm2]');
        grid on
        xlim([min((t-t(1))/ps),max((t-t(1))/ps)]);
        if IO_save==1
            saveas(tmp_fig,[saveKey,'finalOutputPlot.png']);
        end

        tmp_fig=figure(33);
        set(tmp_fig,'Name','Final Output Plot');
        set(tmp_fig,'Position',POS);
        plot((t-t(1))/ps, real(pulse))
        grid on
        xlabel('time [ps]');
        ylabel('real(E)');
        xlim([min((t-t(1))/ps),max((t-t(1))/ps)]);
        if IO_save==1
            saveas(tmp_fig,[saveKey,'finalOutputPlotReal.png']); 
        end
    else
        tmp_fig=figure(233);
        set(tmp_fig,'Name','Final Output (logPlot)');
        set(tmp_fig,'Position',POS);
        contourf((t-t(1))/ps, transverse_grid_short_y/um,log(out_pulse*cm*cm/1e6),'edgecolor','none')
        ylim([-spc_width,spc_width])
        colormap(map3);
        %set(gca,'Linewidth',4);
        xlabel('time [ps]');
        ylabel('y [\mu m]');
        xlim([min((t-t(1))/ps),max((t-t(1))/ps)]);
        colorbar
        if IO_save==1
            saveas(tmp_fig,[saveKey,'finalOutputLog.png']);
        end

        tmp_fig=figure(333);
        set(tmp_fig,'Name','Final Output');
        set(tmp_fig,'Position',POS);
        contourf((t-t(1))/ps, transverse_grid_short_y/um,(out_pulse*cm*cm/1e6),'edgecolor','none')
        ylim([-spc_width,spc_width])
        colormap(map3);
        xlabel('time [ps]');
        ylabel('y [\mu m]');
        %set(gca,'Linewidth',4);
        xlim([min((t-t(1))/ps),max((t-t(1))/ps)]);
        colorbar
        if IO_save==1
            saveas(tmp_fig,[saveKey,'finalOutput.png']);
        end
    end
end


%% Integrate power across transverse dimension
if IO_power==1
disp('Plotting power overlay')
power=zeros(length(t),1);
energy=zeros(length(points),1);
if NUM_TRANSVERSE>1
for j=1:length(t)
power(j)=trapz(transverse_grid_short_y,out_pulse(:,j));
end
else
power=out_pulse;
end
max_power=max(power);
 tmp_fig=figure(222);
  set(tmp_fig,'Position',POS);
for kk=1:num_peaks_plot%length(points)
    t_cur=(t-t(1))/ps-points(kk); %Center t=0 at current peak
%subplot(1,length(points),kk)
% tmp_fig=figure(220+kk);
% set(tmp_fig,'Name',['Power RT',num2str(kk-1)]);
% time=(t-t(1))/ps-points(kk);
    if kk==1
        txt='t_0';
%     elseif kk==2
%         txt='t_0+T_{rt}';
    else
       txt=['t_0+',num2str(round(points(kk)-points(1),1)),'ps'];
%        txt=['t_0+',num2str(kk-1),'T_{rt}'];
    end
 hold on
plot(t_cur(1:length(t)),power(1:length(t))/max_power,my_lineStyle{kk},'DisplayName',txt);
pulse_ind_left=find(t_cur>-tmp_width_power,1);
if isempty(pulse_ind_left)
   pulse_ind_left=1;
end
pulse_ind_right=find(t_cur>tmp_width_power,1);
if isempty(pulse_ind_right)
   pulse_ind_right=length(t_cur);
end
energy(kk)=trapz(t_cur(pulse_ind_left:pulse_ind_right)*ps,power(pulse_ind_left:pulse_ind_right));
energy_ratio=(energy(2:end)-energy(1:end-1))./energy(1:end-1);
%energy(kk)=trapz(t_cur,power);

%ylim([0,1]);

% if kk==1
%     title('t_0','FontSize',24);
% else
%     title(['t_0+',num2str(kk-1),'RT'],'FontSize',24);
% end
 xlim(PowerTLIM);
%     if kk==1
 %   ylabel('Power [a.u]')
%     else
%     set(gca,'YTick', [])
%     end
% end
% hh = mtit(' ');
% set(hh.th,'Visible','off');
% hhx=xlabel(hh.ah,'t [ps]','Fontsize',24);
% set(hhx,'Visible','on');
%xlabel('t[ps]');
end
    set(gca,'YTick', [0,0.25,0.5,0.75,1.0]);
    set(gca,'XTick', PowerTTICK);
      xlabel('time [ps]');
      ylabel('Power [a.u.]');
    hold off
     grid on
     legend show;
     if IO_save==1
        saveas(tmp_fig,[saveKey,'QWpowerOverlay.png']);
     end
end

if IO_transIntensity==1
tmp_fig=figure(100);
set(tmp_fig,'Name','Final Output Pulse Overlay');
hold on
set(tmp_fig,'Position',POS);
for kk=1:num_peaks_plot
    t_cur=(t-t(1))/ps-points(kk);
    pulse_ind_left=find(t_cur>-11,1)-1;
    pulse_ind_right=find(t_cur>11,1);
    if isempty(pulse_ind_left)||pulse_ind_left==0
        pulse_ind_left=1;
    end
    if isempty(pulse_ind_right)
        pulse_ind_right=length(t_cur);
    end
    trans_intensity=zeros(length(transverse_grid_short_y),1);
    for j=1:length(transverse_grid_short_y)
      trans_intensity(j)=trapz(t_cur(pulse_ind_left:pulse_ind_right)*ps,out_pulse(j,pulse_ind_left:pulse_ind_right));
    end
    
    pulse_ind=find(t_cur==0,1);
    if kk==1
        txt='t_0';
%    elseif kk==2
%        txt='t_0+T_{rt}';
    else
       txt=['t_0+',num2str(round(points(kk)-points(1),1)),'ps'];
%       txt=['t_0+',num2str(kk-1),'T_{rt}'];
    end
plot(transverse_grid_short_y/um,trans_intensity*cm*cm*10^3,my_lineStyle{kk},'DisplayName',txt); 
end
xlim([-spc_width,spc_width]);
ylabel('Fluence [mJ/cm^2]');
  xlabel('y [\mu m]');
%ylim([0,round((max(max(out_pulse/10))*cm*cm/1e6)*10)]);
set(gca,'XTick', PulseLogXTICK);
%set(gca,'YTick', [0,round((max(max(out_pulse/10))*cm*cm/1e6)*5),round((max(max(out_pulse/10))*cm*cm/1e6)*10)]);
hold off
grid on
box on
legend show
if IO_save==1
    saveas(tmp_fig,[saveKey,'QWpulseOverlay.png']);
end
end

%% Plots of individual output pulses 
if IO_finalOutput_indPulse==1
for kk=1:length(points)
%A(kk)=subplot(1,length(points),kk);
tmp_fig=figure(1100+kk);
set(tmp_fig,'Name',['Final Output Pulse RT', num2str(kk-1)]);
set(tmp_fig,'Position',POS);
if (NUM_TRANSVERSE > 1)
    t_cur=(t-t(1))/ps-points(kk);
    if max(abs(transverse_grid_short_y))~=0
    %contourf(t_cur, transverse_grid_short_y/um,(out_pulse*cm*cm/1e6),'edgecolor','none');
    contourf(t_cur,transverse_grid_short_y/um,(out_pulse*cm*cm/1e6),'edgecolor','none');
    %ylim([-spc_width,spc_width])
    colormap(map3);
    %caxis([0,2.5]); %WTFF
    %caxis([min(min((out_pulse*cm*cm/1e6))) max(max((out_pulse*cm*cm/1e6)))])
    else
    surf((t-t(1))/ps, transverse_grid_y/um,out_pulse*cm*cm/1e6,'edgecolor','none')
    ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
    end
    xlim(PulseTLIM)
    ylim([-spc_width,spc_width]);
%     if kk==1
%         title('t_0','FontSize',24);
%     else
%         title(['t_0+',num2str(kk-1),'T_{RT}'],'FontSize',24);
%     end
%     if kk==1
    ylabel('y [\mum]')
    xlabel('t [ps]')
%     else
    set(gca,'YTick', [-spc_width,0,spc_width])
%     end
    zlabel('I [MW/cm^2]')
    set(gca,'XTick', PulseTTICK)
    grid on
    B=colorbar;
else
    plot((t-t(1))/ps-points(kk), out_pulse*cm*cm/1e6, 'b-')
    ylabel('I [MW/cm^2]');
    xlim(PulseTLIM);
    xlabel('t [ps]');
end
if IO_save==1
    saveas(tmp_fig,[saveKey,'QWpulseRT',num2str(kk-1),'.png']);
end
end
end

if IO_fourFields==1
    for kk=1:length(points)
        %% FP plots
        tmp_fig=figure(2000+kk);
        set(tmp_fig,'Name',['Final Output Pulse FP RT', num2str(kk-1)]);
        set(tmp_fig,'Position',POS);
        if (NUM_TRANSVERSE > 1)
            t_cur=(t-t(1))/ps-points(kk);
            if max(abs(transverse_grid_short_y))~=0
            contourf(t_cur, transverse_grid_short_y/um,out_pulse_fp,'edgecolor','none');
            %ylim([-spc_width,spc_width])
            colormap(map3);
            %caxis([min(min((out_pulse*cm*cm/1e6))) max(max((out_pulse*cm*cm/1e6)))])
            else
            surf((t-t(1))/ps, transverse_grid_y/um,out_pulse_fp,'edgecolor','none')
            ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
            end
            xlim(PulseTLIM)
        %     if kk==1
        %         title('t_0','FontSize',24);
        %     else
        %         title(['t_0+',num2str(kk-1),'T_{RT}'],'FontSize',24);
        %     end
        %     if kk==1
            ylabel('y [\mum]')
            xlabel('t [ps]')
        %     else
            set(gca,'YTick', [-spc_width,0,spc_width])
        %     end
            zlabel('|E| [V/m]')
            set(gca,'XTick', PulseTTICK)
            grid on
        else
            plot((t-t(1))/ps, abs(pulse_fp), 'b-')
            ylabel('|E| [V/m]');
        end
        B=colorbar;
        if IO_save==1
            saveas(tmp_fig,[saveKey,'QWpulseFP_RT',num2str(kk-1),'.png']);
        end
        
        %% FM plots
        tmp_fig=figure(3000+kk);
        set(tmp_fig,'Name',['Final Output Pulse FM RT', num2str(kk-1)]);
        set(tmp_fig,'Position',POS);
        if (NUM_TRANSVERSE > 1)
            t_cur=(t-t(1))/ps-points(kk);
            if max(abs(transverse_grid_short_y))~=0
            contourf(t_cur, transverse_grid_short_y/um,out_pulse_fm,'edgecolor','none');
            %ylim([-spc_width,spc_width])
            colormap(map3);
            %caxis([min(min((out_pulse*cm*cm/1e6))) max(max((out_pulse*cm*cm/1e6)))])
            else
            surf((t-t(1))/ps, transverse_grid_y/um,out_pulse_fm,'edgecolor','none')
            ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
            end
            xlim(PulseTLIM)
        %     if kk==1
        %         title('t_0','FontSize',24);
        %     else
        %         title(['t_0+',num2str(kk-1),'T_{RT}'],'FontSize',24);
        %     end
        %     if kk==1
            ylabel('y [\mum]')
            xlabel('t [ps]')
        %     else
            set(gca,'YTick', [-spc_width,0,spc_width])
        %     end
            zlabel('|E| [V/m]')
            set(gca,'XTick', PulseTTICK)
            grid on
        else
            plot((t-t(1))/ps, abs(pulse_fm), 'b-')
            ylabel('|E| [V/m]');
        end
        B=colorbar;
        if IO_save==1
            saveas(tmp_fig,[saveKey,'QWpulseFM_RT',num2str(kk-1),'.png']);
        end
        
        %% BP plots
        tmp_fig=figure(4000+kk);
        set(tmp_fig,'Name',['Final Output Pulse BP RT', num2str(kk-1)]);
        set(tmp_fig,'Position',POS);
        if (NUM_TRANSVERSE > 1)
            t_cur=(t-t(1))/ps-points(kk);
            if max(abs(transverse_grid_short_y))~=0
            contourf(t_cur, transverse_grid_short_y/um,out_pulse_bp,'edgecolor','none');
            %ylim([-spc_width,spc_width])
            colormap(map3);
            %caxis([min(min((out_pulse*cm*cm/1e6))) max(max((out_pulse*cm*cm/1e6)))])
            else
            surf((t-t(1))/ps, transverse_grid_y/um,out_pulse_bp,'edgecolor','none')
            ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
            end
            xlim(PulseTLIM)
        %     if kk==1
        %         title('t_0','FontSize',24);
        %     else
        %         title(['t_0+',num2str(kk-1),'T_{RT}'],'FontSize',24);
        %     end
        %     if kk==1
            ylabel('y [\mum]')
            xlabel('t [ps]')
        %     else
            set(gca,'YTick', [-spc_width,0,spc_width])
        %     end
            zlabel('|E| [V/m]')
            set(gca,'XTick', PulseTTICK)
            grid on
        else
            plot((t-t(1))/ps, abs(pulse_bp), 'b-')
            ylabel('|E| [V/m]');
        end
        B=colorbar;
        if IO_save==1
            saveas(tmp_fig,[saveKey,'QWpulseBP_RT',num2str(kk-1),'.png']);
        end
        
        %% BM plots
        tmp_fig=figure(2000+kk);
        set(tmp_fig,'Name',['Final Output Pulse BM RT', num2str(kk-1)]);
        set(tmp_fig,'Position',POS);
        if (NUM_TRANSVERSE > 1)
            t_cur=(t-t(1))/ps-points(kk);
            if max(abs(transverse_grid_short_y))~=0
            contourf(t_cur, transverse_grid_short_y/um,out_pulse_bm,'edgecolor','none');
            %ylim([-spc_width,spc_width])
            colormap(map3);
            %caxis([min(min((out_pulse*cm*cm/1e6))) max(max((out_pulse*cm*cm/1e6)))])
            else
            surf((t-t(1))/ps, transverse_grid_y/um,out_pulse_bm,'edgecolor','none')
            ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
            end
            xlim(PulseTLIM)
        %     if kk==1
        %         title('t_0','FontSize',24);
        %     else
        %         title(['t_0+',num2str(kk-1),'T_{RT}'],'FontSize',24);
        %     end
        %     if kk==1
            ylabel('y [\mum]')
            xlabel('t [ps]')
        %     else
            set(gca,'YTick', [-spc_width,0,spc_width])
        %     end
            zlabel('|E| [V/m]')
            set(gca,'XTick', PulseTTICK)
            grid on
        else
            plot((t-t(1))/ps, abs(pulse_bm), 'b-')
            ylabel('|E| [V/m]');
        end
        B=colorbar;
        if IO_save==1
            saveas(tmp_fig,[saveKey,'QWpulseBM_RT',num2str(kk-1),'.png']);
        end
    end
end
% 
% for kk=1:length(points)
%     tmp_fig=figure(2000+kk);
% set(tmp_fig,'Name',['Final Output Pulses (Log) RT',num2str(kk)]);
% set(tmp_fig,'Position',POS);
% %A(kk)=subplot(1,length(points),kk);
% if (NUM_TRANSVERSE > 1)
%     if max(abs(transverse_grid_short_y))~=0
%     contourf((t-t(1))/ps-points(kk), transverse_grid_short_y/um,log(out_pulse*cm*cm/1e6),'edgecolor','none');
%     ylim([-spc_width,spc_width])
%     colormap(map3);
%     caxis([min(min(log(out_pulse*cm*cm/1e6))) max(max(log(out_pulse*cm*cm/1e6)))])
%     else
%     surf((t-t(1))/ps, transverse_grid_y/um,out_pulse*cm*cm/1e6,'edgecolor','none')
%     ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
%     end
%     xlim(PulseLogTLIM)
% %     if kk==1
% %         title('t_0','FontSize',24);
% %     else
% %         title(['t_0+',num2str(kk-1),'RT'],'FontSize',24);
% %     end
% %     if kk==1
% %     ylabel('y [\mum]','Fontsize', 24)
% %     else
% %    set(gca,'YTick', [])
% %     end
%     zlabel('I [MW/cm^2]')
%     set(gca,'XTick', PulseLogTTICK);
%     set(gca,'YTick', PulseLogXTICK);
% else
%     plot((t-t(1))/ps, out_pulse*cm*cm/1e6, 'b-')
%     ylabel('I [MW/cm^2]');
% end
% %   hh = mtit(' ');
% % set(hh.th,'Visible','off');
% % hhx=xlabel(hh.ah,'t [ps]','Fontsize',24);
% % set(hhx,'Visible','on');
% B=colorbar;
% % set(B, 'Position', [.89 .11 .025 .8150])
% % for i=1:length(points)
% % pos=get(A(i), 'Position'); 
% % axes(A(i))
% % set(A(i), 'Position', [pos(1) pos(2) .09 pos(4)])
% % end
% saveas(tmp_fig,[saveKey,'Qwpulses_logRT',num2str(kk-1),'.png']);
% end
%% Peak intensity and FWHM (time and transverse dim) and Pulse energy
if IO_convergence==1
    if ~(strcmp(location,'OUTPUT')||strcmp(location,'OUTPUTBACK'))
       error('Location must be at OUTPUT/OUTPUTBACK for convergence study')
    end
disp('Plotting convergence figures')
max_I = zeros(counter, plot_num);
max_It = zeros(1,plot_num);
fwhm_t = zeros(1,plot_num);
fwhm_y = zeros(1,plot_num);
energy_evol=zeros(1,plot_num);
%temperature_profile = zeros(NUM_TRANSVERSE_DEVICES_SHORT, 1);
%density_profile = zeros(NUM_TRANSVERSE_DEVICES_SHORT, 1);
%transverse_grid_short_y=zeros(NUM_TRANSVERSE_DEVICES_SHORT,1);
power=zeros(length(t),1);
out_pulse=0;

for j = 0:(plot_num-1)
    current_plot=j;
    t = loadD([outKey,num2str(j),'_t.dat']);
    trans_profile = zeros(NUM_TRANSVERSE,length(t));
    counter=0;
    NUM_TRANSVERSE
    for i = 0:(NUM_TRANSVERSE-1)
     if abs(transverse_grid_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
        i
        pulse_re = loadD([outKey,num2str(j),'_E_re_',location,'_T',num2str(i),'.dat']);
        pulse_im = loadD([outKey,num2str(j),'_E_im_',location,'_T',num2str(i),'.dat']);
        %tmp = loadD([outKey,'lattice_setup_','QW6_T',num2str(ind_device_y(1+i)),'.dat']);
        %density_profile(1+counter) = tmp(1);
        %temperature_profile(1+counter) = tmp(2);
        pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*w0);
        pulse_tmp(1+counter,1:length(pulse))= pulse;
        out_pulse(1+counter,1:length(pulse)) = 0.5*n*eps0*c0*abs(pulse).^2;
        out_pulse(1+counter,length(pulse)+1:end)=0;
       %transverse_grid_short_y(1+counter)=transverse_grid_device_y(1+i);
        max_I(1+counter,1+j) = max(0.5*n*eps0*c0*abs(pulse).^2);
        counter=counter+1;
     end
    end
% Find a pulse to focus on
        [~,ind] = findpeaks(out_pulse(ceil(counter/2),:),'sortstr','descend');
        ind = sort(ind(1:min(2,length(ind))));        
        pulse_ind_left=find((t-t(ind(1)))/ps>-tmp_width_power,1);
        if isempty(pulse_ind_left)
            pulse_ind_left=1;
        end
        pulse_ind_right=find((t-t(ind(1)))/ps>tmp_width_power,1);
        if isempty(pulse_ind_right)
            pulse_ind_right=length(t);
        end
        if (NUM_TRANSVERSE>1)
            fwhm_y(1+j) = findFWHM(transverse_grid_y, out_pulse(:,ind(1)) );
            %fwhm_y(1+j) = findSpotSize(transverse_grid_y,abs(trans_profile(:,ind(1))));
            fwhm_t(1+j) = findFWHM(t, out_pulse(ceil(counter/2),:) );
            power=zeros(pulse_ind_right-pulse_ind_left+1,1);
            for jj=pulse_ind_left:pulse_ind_right
                power(1+jj-pulse_ind_left)=trapz(transverse_grid_short_y,out_pulse(:,jj));
            end
            energy_evol(j+1)=trapz(t(pulse_ind_left:pulse_ind_right),power);
        else
            fwhm_t(1+j) = findFWHM(t, out_pulse );
            energy_evol(j+1)=trapz(t(pulse_ind_left:pulse_ind_right),out_pulse(pulse_ind_left:pulse_ind_right));
        end
        max_It(1+j) = mean(t);
end

yy = max_I*cm*cm/1e6;
%yy = log(yy/max(yy(:)));

tmp_fig=figure(101);
set(tmp_fig,'Name','Convergence-Profile');
set(tmp_fig,'Position',POS);
%subplot(4,1,1)
if (NUM_TRANSVERSE>1)
    surf(max_It/ns,transverse_grid_short_y/um,yy,'edgecolor','none')
    ylim([min(transverse_grid_short_y),max(transverse_grid_short_y)]/um)
%     ylabel('y [\mum]','Fontsize',24)
%     zlabel('I(t) [MW/cm^2]','Fontsize',24)
    view(25,22)
else
    plot(max_It/ns, yy,'color','b','Linestyle','-', 'marker', 'o');
    xlabel('time [ns]');
    ylabel('Intensity [MW/cm^{2}]');
%     ylabel('I(t) [MW/cm^2]','Fontsize',24)
end
%set(gca,'FontSize', 24)

zl = zlim;
%hold on
%yy2 = (density_profile-min(density_profile))/(max(density_profile)-min(density_profile));
%plot3(0*ones(size(transverse_grid_short_y))/ps, transverse_grid_short_y/um, zl(1) + yy2*(zl(2)-zl(1)),'ko-')
%hold off

%xlabel('t [ns]')
grid on
if IO_save==1
    saveas(tmp_fig,[saveKey,'QWconvergence-Profile.png']);
end

%Plots central intensity
tmp_fig=figure(102);
set(tmp_fig,'Name','Convergence-PeakIntensity');
set(tmp_fig,'Position',POS);
if (NUM_TRANSVERSE>1)
    plot(max_It/ns,yy(ceil(counter/2),:),'bo-');
    xlabel('time [ns]');
    ylabel('Intensity [MW/cm^{2}]');
else
    plot(max_It/ns, yy,'color','b','Linestyle','-', 'marker', 'o');
    xlabel('time [ns]');
    ylabel('Intensity [MW/cm^{2}]');
end
grid on
% xlabel('t [ns]','Fontsize',24)
% ylabel('I(t) [MW/cm^2]','Fontsize',24)
if IO_save==1
    saveas(tmp_fig,[saveKey,'QWconvergence-PeakIntensity.png']);
end

%Plots FWHM in time and space
if NUM_TRANSVERSE>1
    tmp_fig=figure(103);
    set(tmp_fig,'Name','FWHM');
    set(tmp_fig,'Position',POS);
    [AX,H1,H2] = plotyy(max_It/ns,fwhm_t/fs, max_It/ns,fwhm_y/um);
    set(AX,{'ycolor'},{'b';'r'})
    set(H1,'color','b','LineStyle','-','marker','o')
    set(H2,'color','r','LineStyle','-','marker','o')
    grid on
    axes(AX(1))
    xlabel('t [ns]')
    ylabel('FWHM [fs]')
    axes(AX(2))
    ylabel('FWHM [\mum]')
    
    
else
    tmp_fig=figure(103);
    set(tmp_fig,'Name','FWHM');
    set(tmp_fig,'Position',POS);
    plot(max_It/ns,fwhm_t/fs,'color','b','Linestyle','-', 'marker', 'o');
    xlabel('time [ns]');
    ylabel('FWHM [fs]');
    %set(AX,{'ycolor'},{'b'})
    %set(H1,'color','b','LineStyle','-','marker','o')
    %set(H2,'color','r','LineStyle','-','marker','o')
    grid on
end

if IO_save==1
    saveas(tmp_fig,[saveKey,'QWconvergence-FWHM.png']);
end

tmp_fig=figure(104);
set(tmp_fig,'Name','Convergence-Energy');
set(tmp_fig,'Position',POS);
hold on
plot(max_It/ns,energy_evol/(10^6),'bo-');
%plot(max_It(17)/ns,energy_evol(17)/(10^6),'ro-');
%plot(max_It(20)/ns,energy_evol(20)/(10^6),'mo-');
%plot(max_It(23)/ns,energy_evol(23)/(10^6),'ko-');
hold off
xlabel('time [ns]');
ylabel('Energy [J/cm^2]');
grid on
% xlabel('t [ns]')
% ylabel('Energy [a.u.]')
% set(gca,'FontSize', 24)
if IO_save==1
    saveas(tmp_fig,[saveKey,'QWconvergence-Energy.png']);
end
end

%% Plot QW density and temperature
if IO_temp_dens==1
    disp('Plotting density figures')
    qw_density  = zeros(NUM_TRANSVERSE_DEVICES_SHORT,length(t));
    qw_temp     = zeros(NUM_TRANSVERSE_DEVICES_SHORT,length(t));
    abs_density = zeros(NUM_TRANSVERSE_DEVICES_SHORT,length(t));
    abs_temp    = zeros(NUM_TRANSVERSE_DEVICES_SHORT,length(t));
    counter     = 1;

    for i = 0:(NUM_TRANSVERSE_DEVICE-1)
        if abs(transverse_grid_device_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
          Ne = loadD([outKey,num2str(plot_num-1),'_Nsum_e_',location_short,'_T',num2str(ind_device_y(1+i)),'.dat']);
        %  Te = loadD([outKey,num2str(plot_num-1),'_inst_temp_e_',location_short,'_T',num2str(ind_device_y(1+i)),'.dat']);
          Ne=Ne(init_ind:final_ind);
        %  Te=Te(init_ind:final_ind);

          qw_density(counter,:) = Ne;
         % qw_temp(counter,:) = Te;    
          counter=counter+1;
        end
    end

    if NUM_TRANSVERSE>1
        tmp_fig=figure(111);
        set(tmp_fig,'Name','Density');
        set(tmp_fig,'Position',POS);
        x_ind_left=find(transverse_grid_short_y/um>=-1.05*spc_width_density,1);
        x_ind_right=length(transverse_grid_short_y);
        contourf((t-t(1)-points(1)*ps)/ps,transverse_grid_short_y(x_ind_left:x_ind_right)/um,...
        qw_density(x_ind_left:x_ind_right,:)/10^16,'edgecolor','none')
        %caxis([min(min(qw_density(x_ind_left:x_ind_right,:))), max(max(qw_density(x_ind_left:x_ind_right,:)))]);
        %caxis([1.9,2.0]);
        ylim([-spc_width_density,spc_width_density]);
        %xlim([-0.5,points(end)-points(1)+1]);
        set(gca,'YTick', [-spc_width_density,0,spc_width_density]);
        tmp=(t(end)-t(1))/ps-points(1);
        ind=find((t-t(1)+tmp*ps)>=round_trip_time/2,1);
        time_ticks=[(t(ind)-t(1))/ps-points(1);round(points-points(1),1);(t(end)-t(1))/ps-points(1)];
        set(gca, 'XTick',time_ticks, 'xticklabel',{'t_0','t_1','t_2','t_0+T_{rt}/2'})
        xlim([(t(ind)-t(1))/ps-points(1),((t(end)-t(1))/ps-points(1))]);
        xlabel('time [a.u.]');
        ylabel('y [\mu m]');
        %colormap('hot')
        colorbar
    else
        tmp_fig=figure(123);
        set(tmp_fig,'Name','Density');
        set(tmp_fig,'Position',POS); 
        plot(t_cur,Ne*(10^(-16)))
        xlabel('time [ps]');
        xlim([t_cur(1),t_cur(end)]);
        ylabel('Electron Density [10^{16}m^{-2}]');
    end
    if IO_save==1
        saveas(tmp_fig,[saveKey,'QWdensity.png']);
    end
    
    if IO_densityVideo==1
        if NUM_TRANSVERSE>1
            tmp_fig=figure(1010);
            set(tmp_fig,'Name','Density Video Plot');
            cnt=0;
            if t_cur(1)>=0
               t_cur=t_cur-points(1);  %Center on first pulse
            end
            for j=1750:floor(length(t_cur)/1000):3250 %200 frames
                plot(transverse_grid_short_y/um,qw_density(:,j)/1.0E16);
                hold on
                tmp=out_pulse(:,j);
                maxIntens=max(max(out_pulse));
                normedIntes=tmp/maxIntens;
                scale=2.375-2.25;
                shift=2.25;
                tmpScaled=normedIntes*scale*0.95+shift;
                plot(transverse_grid_short_y/um,tmpScaled);
                hold off
                %xlim([min(transverse_grid_short_y),max(transverse_grid_short_y)]/um)
                xlim([-175,175]);
                ylim([2.25,2.375]);
                ylabel(' Density [10^{16} m^{-2}]');
                xlabel('y [\mu m]');
                txt=text(50,2.3675,['time=',num2str(t_cur(j)),'ps']);
                drawnow
                cnt=cnt+1;
                ax = gca;
                ax.Units = 'pixels';
                pos = ax.Position;
                ti = ax.TightInset;
                rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
                ax.Units = 'normalized';
                F(cnt)=getframe(gcf);
            end
            if IO_save==1
                myVideo = VideoWriter([saveKey,'INVvideo.avi']);
                myVideo.FrameRate = 20;
                open(myVideo);
                writeVideo(myVideo, F);
                close(myVideo);
            end
        else
            disp('Density video in 1D does not work')
        end
    end

    % % 
    % figure
    % subplot(2,1,1)
    % plot(transverse_grid_device_y/um,min(qw_density,[],2),'bo-',...
    %      transverse_grid_device_y/um,max(qw_density,[],2),'ko-')
    %  grid on
    %  xlabel('y [um]')
    %  ylabel('Density [-]')
    %  legend('min','max')
    %  
    % subplot(2,1,2)
    % plot(transverse_grid_device_y/um,min(qw_temp,[],2),'bo-',...
    %      transverse_grid_device_y/um,max(qw_temp,[],2),'ko-')
    %  grid on
    %  xlabel('y [um]')
    %  ylabel('Temperature [-]')
    %  legend('min','max')
    % %  
    %  
    % figure
    % subplot(2,1,1)
    % plot(transverse_grid_device_y/um,min(abs_density,[],2),'bo-',...
    %      transverse_grid_device_y/um,max(abs_density,[],2),'ko-')
    %  grid on
    %  xlabel('y [um]')
    %  ylabel('Density [-]')
    %  legend('min','max')
    %  
    % subplot(2,1,2)
    % plot(transverse_grid_device_y/um,min(abs_temp,[],2),'bo-',...
    %      transverse_grid_device_y/um,max(abs_temp,[],2),'ko-')
    %  grid on
    %  xlabel('y [um]')
    %  ylabel('Temperature [-]')
    %  legend('min','max')
    % 
    %  
    % 
    % 
    % for i = 1:length(t)
    %     
    %     abs_density(:,i) = abs_density(:,i) - abs_density(:,end);
    %     abs_temp(:,i) = abs_temp(:,i) - abs_temp(:,end);
    %    
    %     qw_density(:,i) = qw_density(:,i) - qw_density(:,end);
    %     qw_temp(:,i) = qw_temp(:,i) - qw_temp(:,end);
    % end
    % 
    % figure
    % subplot(2,1,1)
    % if (NUM_TRANSVERSE>1)
    %     surf((t-t(1))/ps,transverse_grid_device_y/um,qw_density/1e16,'edgecolor','none')
    %     ylabel('y [um]')
    % else
    %       plot((t-t(1))/ps, qw_density/1e16)
    %     ylabel('N [10^{16}m^{-2}]')
    % end
    % xlabel('t [ps]')
    % title('Density QW1 [1e16]')
    % colorbar
    % grid on
    % 
    % subplot(2,1,2)
    % if (NUM_TRANSVERSE>1)
    %     surf((t-t(1))/ps,transverse_grid_device_y/um,qw_temp,'edgecolor','none')
    %     ylabel('y [um]')
    % else
    %       plot((t-t(1))/ps, qw_temp)
    %     ylabel('T [K]')
    % end
    % xlabel('t [ps]')
    % title('Temperature QW1 [K]')
    % colorbar
    % grid on
    % 
    % figure
    % subplot(2,1,1)
    % if (NUM_TRANSVERSE>1)
    %     surf((t-t(1))/ps,transverse_grid_device_y/um,abs_density/1e14,'edgecolor','none')
    %     ylabel('y [um]')
    % else
    %       plot((t-t(1))/ps, abs_density/1e14)
    %     ylabel('N [10^{14}m^{-2}]')
    % end
    % xlabel('t [ps]')
    % title('Density QW6 [1e14]')
    % colorbar
    % grid on
    % 
    % subplot(2,1,2)
    % if (NUM_TRANSVERSE>1)
    %     surf((t-t(1))/ps,transverse_grid_device_y/um,abs_temp,'edgecolor','none')
    %     ylabel('y [um]')
    % else
    %       plot((t-t(1))/ps, abs_temp)
    %     ylabel('T [K]')
    % end
    % xlabel('t [ps]')
    % title('Temperature QW6 [K]')
    % colorbar
    % grid on
    % end

    %FFT transform
    % dy=transverse_grid_short_y(2)-transverse_grid_short_y(1);
    % pad=50;
    % [w,y]=getFFT(pulse_tmp(:,1000),dy,pad);
    % output_FFT=zeros(length(y),length(t));
    % for j=1:length(t)
    % [w,y]=getFFT(pulse_tmp(:,j),dy,pad);
    % output_FFT(:,j)=y;
    % end
    % 
    % figure
    % surf((t-t(1))/ps,w,abs(output_FFT),'edgecolor','none')
    % xlim([0,(t(end)-t(1))/ps]);
    % xlabel('time (ps)','Fontsize',16);
    % ylim([min(w),max(w)/3]);
    % ylabel('Frequency (a.u.)','Fontsize',16);
    % title('FFT of output .','Fontsize', 18);
    % asd
    % 

    % figure
    % subplot(2,1,1)
    % surf((t-t(1))/ps, transverse_grid_short_y/um,real(pulse_tmp),'edgecolor','none')
    %     ylim([min( transverse_grid_short_y),max(transverse_grid_short_y)]/um)
    % grid on
    % xlabel('t [ps]')
    % xlim([min(t-t(1)), max(t-t(1))]/ps)
    % 
    % subplot(2,1,2)
    % surf((t-t(1))/ps, transverse_grid_short_y/um,imag(pulse_tmp),'edgecolor','none')
    %     ylim([min( transverse_grid_short_y),max(transverse_grid_short_y)]/um)
    % grid on
    % xlabel('t [ps]')
    % xlim([min(t-t(1)), max(t-t(1))]/ps)
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
try
[maxp, maxi] = max(pulse); % max pulse

% Find left boundary
i0 = maxi;
while ((i0 > 2)&&(pulse(i0-1) > 0.5*maxp))

    i0 = i0 -1;
end

% Find right boundary
i1 = maxi;
while ((i1 < length(pulse)-1)&&(pulse(i1+1) > 0.5*maxp))

    i1 = i1 +1;
end
%fwhm = DT(i1) - DT(i0);


% Interpolate data
ind0 = i0-1:i0+1;
T0 = invertYData(DT(ind0),pulse(ind0),0.5*maxp);

ind1 = i1-1:i1+1;
T1 = invertYData(DT(ind1),pulse(ind1),0.5*maxp);


fwhm = abs(T1 - T0);
catch
    fwhm=0.0;
end
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

