%%Comparisons for outputs/evolution of TMSBE modelocking simulations
%S.A.Mclaren (April 27, 2021)
% close all
% clear all

%IOs for each plot located here
IO_CWandCCW=1; %Not fully configured
IO_convergence=1; %Plot convergence of FWHM, peak intensity, and energy
IO_spectrums=0; % Plot spectrums of particular peaks. Not currently configure
%IO_temp_dens=0; %Plot density and temperature (Not currently configured)
%IO_finalOutput=0; %Plot final output (Not currently configured)
%IO_finalOutput_indPulse=1; %Individual plots of pulses at final output (Not currently configured)
IO_power=1; %Plot power overlay for latest output
%IO_transIntensity=0; %Plot transverse intensity profiles (Not currently configured)
IO_pointDiff=0; %Plot time differences between peaks for consecutive simulations (RCAV setup). Requires IO_convergence
%IO_findFWHM=0; %Find FWHM for each peak (saved in pulseFWHM) (Not currently configured)
IO_save=1; %Save figures
IO_separatePeaks=0; %Ensure peaks are distinct by set amount

%Setup
setupPlot
setupConstants
setupSave

width=[2.5,2.5]; %Time on either side of peak in ps
width_short=[-4.0,10.0];
width_spectrum=[1.8,2.2];
YLIM_intensity=[10^(-13);10^(1)];
XLIM_SPECTRUM=[1.1,1.3];
PowerTLIM=[-0.5,0.5]; %Power plot  temporal width
PowerTTICK=[-0.5,0.0,0.5]; %Power plot tick marks
sizeIl=500*fs;
sizeIr=2000*fs;

if IO_CWandCCW==1
    lineStyle_mod=2;
    lineStyle_shift=1;
else
    lineStyle_mod=1;
    lineStyle_shift=1;
end

%% All load locations
% VCAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV55-control-equalArms-n2p0-noThresh/run/out__'};
% VCAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV62-control-equalArms-n2p0-noThresh-Ny288/run/out__'};
% VCAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV57-totalThresh-equalArms-n2p0/run/out__'};
% VCAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV58-fieldFreeThresh-equalArms-n2p0/run/out__'};
% VCAVsim(5).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV56-fullThresh-equalArms-n2p0/run/out__'};
% VCAVsim(6).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV64-equalArms-n2p0-fullThresh-2em2-Ny288-theta2/run/out__'};
% VCAVsim(7).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV60-fullThreshOnly-Thresh-2e-8-equalArms-n2p0/run/out__'};
% 
% VCAVsim(1).legend={'Ny=144, No thresh., VCAV55'};
% VCAVsim(2).legend={'Ny=288, No thresh., VCAV62'};
% VCAVsim(3).legend={'Ny=144, Coulomb-free thresh., VCAV57'};
% VCAVsim(4).legend={'Ny=144, Field-free and Coulomb-free, VCAV58'};
% VCAVsim(5).legend={'Ny=144, All thresh., T_0=2, VCAV56'};
% VCAVsim(6).legend={'Ny=288, All thresh., T_0=2e-2, VCAV64'};
% VCAVsim(7).legend={'Ny=144, Only Analytic, VCAV60'};

% VCAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV59-control-equalArms-n2p0-noThresh-theta4/run/out__'};
% VCAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV59-control-equalArms-n2p0-noThresh-theta4-redo/run/out__'};
% VCAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV61-fullThresh-equalArms-n2p0-theta4/run/out__'};
% VCAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV65-control-equalArms-n2p0-Ny288-theta4/run/out__'};
% VCAVsim(5).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV63-equalArms-n2p0-fullThresh-2em2-Ny288-theta4/run/out__'};
% 
% VCAVsim(1).legend={'Ny=144, No thresh., VCAV59'};
% VCAVsim(2).legend={'Ny=144, No thresh., redo, VCAV59A'};
% VCAVsim(3).legend={'Ny=144, All thresh., T_0=2, VCAV61'};
% VCAVsim(4).legend={'Ny=288, No thresh., VCAV65'};
% VCAVsim(5).legend={'Ny=288, All thresh., T_0=2e-2, VCAV63'};

% VCAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV64-equalArms-n2p0-fullThresh-2em2-Ny288-theta2/run/out__'};
% VCAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV63-equalArms-n2p0-fullThresh-2em2-Ny288-theta4/run/out__'};
% VCAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV68-equalArms-n2p0-fullThresh-2em2-Ny288-theta8/run/out__'};
% VCAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV75-equalArms-n2p0-fullThresh-2em2-Ny288-theta12/run/out__'};
% VCAVsim(5).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV69-equalArms-n2p0-fullThresh-2em2-Ny288-theta16/run/out__'};
% VCAVsim(6).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV76-equalArms-n2p0-fullThresh-2em2-Ny288-theta20/run/out__'};
%VCAVsim(7).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV66-Arms_3to1-n2p0-fullThresh-2em2-Ny288-theta4/run/out__'};
%VCAVsim(7).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV67-Arms7to1-n2p0-fullThresh-2em2-Ny288-theta4/run/out__'};


% VCAVsim(1).legend={'\theta=2^\circ'};
% VCAVsim(2).legend={'\theta=4^\circ'};
% VCAVsim(3).legend={'\theta=8^\circ'};
% VCAVsim(4).legend={'\theta=12^\circ'};
% VCAVsim(5).legend={'\theta=16^\circ'};
% VCAVsim(6).legend={'\theta=20^\circ'};
%VCAVsim(7).legend={'Arms 3-1, theta=4, VCAV66'};
%VCAVsim(7).legend={'Arms 7-1, \theta=4^\circ, VCAV67'};

% VCAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV63-equalArms-n2p0-fullThresh-2em2-Ny288-theta4/run/out__'};
% VCAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV67-Arms7to1-n2p0-fullThresh-2em2-Ny288-theta4/run/out__'};
% VCAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV86-equalArms-6400lambda-n2p0-fullThresh-2em2-Ny288-theta4/run/out__'};
% %VCAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV89-ArmsTest1-SESAM-5600-n2p0-fullThresh-2em2-Ny288-theta4/run/out__'};
% VCAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV85-Arms7to1-6400lambda-n2p0-fullThresh-2em2-Ny288-theta4/run/out__'};
% % 
% VCAVsim(1).legend={'3.4mm, 1:1'};
% VCAVsim(2).legend={'3.4mm, 7:1'};
% VCAVsim(3).legend={'6.7mm, 1:1'};
% VCAVsim(4).legend={'6.7mm, 7:1'};
%VCAVsim(5).legend={'800\lambda, 5600\lambda'};
% 
% VCAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV66-Arms_3to1-n2p0-fullThresh-2em2-Ny288-theta4/run/out__'};
% VCAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-VCAV66-Arms_3to1-n2p0-fullThresh-2em2-Ny288-theta4-redoFixedOutputFreq-restart/run/out__'};
% VCAVsim(1).legend={'Unstable simulation'};
% VCAVsim(2).legend={'Stable simulation'};

% VCAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV37-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5/run/out__'};
% VCAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV37-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5/run/out__'};
% VCAVsim(1).legend={'CW'};
% VCAVsim(2).legend={'CCW'};
% VCAVsim(1).location={'OUTPUT'};
% VCAVsim(2).location={'OUTPUTBACK'};

% VCAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-VCAV103-2D-n2p35-allThresh-theta2-3200lam-noSpont-noExpSBE-focus5-lessCores/run/out__'};
% VCAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-VCAV104-2D-n2p35-noThresh-theta2-3200lam-noSpont-noExpSBE-focus5-lessCores/run/out__'};
% VCAVsim(1).legend={'NO THRESH'};
% VCAVsim(2).legend={'THRESH'};
% VCAVsim(1).location={'OUTPUT'};
% VCAVsim(2).location={'OUTPUT'};

%VCAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV37-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5/run/out__'};
%VCAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV37-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5/run/out__'};
%VCAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV95-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5-seedPulse-2e5/run/out__'};
%VCAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV95-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5-seedPulse-2e5/run/out__'};
VCAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV99-2D-n2p5-fullThresh-theta2-1700-1500-3200lam-noSpontEmis-wExpSBE-focus5-seedPulse-2e5/run/out__'};
VCAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV99-2D-n2p5-fullThresh-theta2-1700-1500-3200lam-noSpontEmis-wExpSBE-focus5-seedPulse-2e5/run/out__'};
%VCAVsim(7).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV110-2D-n2p5-fullThresh-theta2-5100-4500-9600lam-noSpontEmis-wExpSBE-focus5-seedPulse-2e5/run/out__'};
%VCAVsim(8).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV110-2D-n2p5-fullThresh-theta2-5100-4500-9600lam-noSpontEmis-wExpSBE-focus5-seedPulse-2e5/run/out__'};
VCAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV122-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5-seedPulse-2e5-highRes/run/out__'};
VCAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV122-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5-seedPulse-2e5-highRes/run/out__'};


%VCAVsim(1).legend={'SPONT CW'};
%VCAVsim(2).legend={'CCW'};
%VCAVsim(3).legend={'SEED CW'};
%VCAVsim(4).legend={'CCW'};
VCAVsim(1).legend={'N_y=216, CW'};
VCAVsim(2).legend={'CCW'};
VCAVsim(3).legend={'N_y=336, CW'};
VCAVsim(4).legend={'CCW'};
VCAVsim(1).location={'OUTPUT'};
VCAVsim(2).location={'OUTPUTBACK'};
VCAVsim(3).location={'OUTPUT'};
VCAVsim(4).location={'OUTPUTBACK'};



for j=1:4
  VCAVsim(1).type={'RCAV'};  
end

num_sims=length(VCAVsim);

for j = 1 : num_sims
   %VCAVsim(j).location=location; %Location for output. Should be the same for all
   tmp_cnt=1; %temporary counter for plot num
   while isfile([char(VCAVsim(j).outKey),num2str(tmp_cnt-1),'_E_re_',char(VCAVsim(j).location),'_T1.dat'])
          tmp_cnt=tmp_cnt+1; %Check for next plot. If exists, add to counter
   end
   VCAVsim(j).plot_num=tmp_cnt-1; %number of output plots for convergence
   %VCAVsim(j).plot_point=2; %which peak to plot
   VCAVsim(j).spc_width=250; %Spatial width for loading in microns
   VCAVsim(j).NUM_TRANSVERSE = length(dir([char(VCAVsim(j).outKey),num2str(VCAVsim(j).plot_num-1),'_E_re_OUTPUT_T*.dat']));
   VCAVsim(j).round_trip_time = loadD([char(VCAVsim(j).outKey),'round_trip_time.dat']);
   if strcmp(VCAVsim(j).type,'RCAV')
     VCAVsim(j).round_trip_time=VCAVsim(j).round_trip_time/2;  
   end
   VCAVsim(j).transverse_grid_y = loadD([char(VCAVsim(j).outKey),'transverse_grid_y.dat']);
   VCAVsim(j).dx=VCAVsim(j).transverse_grid_y(2)-VCAVsim(j).transverse_grid_y(1);
   VCAVsim(j).tmp_width_power=1.0;%0.1*VCAVsim(j).round_trip_time/(ps);
   VCAVsim(j).w0 = loadD([char(VCAVsim(j).outKey),'w0.dat']);
end   

if IO_convergence==1
    % Load and find peaks for all outKeys
    for ll = 1:num_sims
        ll
        for j = 0:(VCAVsim(ll).plot_num-1)
            j
            current_plot=j;
            t = loadD([char(VCAVsim(ll).outKey),num2str(j),'_t.dat']);
            VCAVsim(ll).counter=0;
            for i = 0:(VCAVsim(ll).NUM_TRANSVERSE-1)
                if abs(VCAVsim(ll).transverse_grid_y(1+i))<VCAVsim(ll).spc_width*um+VCAVsim(ll).dx %max(abs(transverse_grid_y))/8
                    pulse_re = loadD([char(VCAVsim(ll).outKey),num2str(j),'_E_re_',char(VCAVsim(ll).location),'_T',num2str(i),'.dat']);
                    pulse_im = loadD([char(VCAVsim(ll).outKey),num2str(j),'_E_im_',char(VCAVsim(ll).location),'_T',num2str(i),'.dat']);
                    VCAVsim(ll).pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*VCAVsim(ll).w0);
                    VCAVsim(ll).out_pulse(1+VCAVsim(ll).counter,1:length(VCAVsim(ll).pulse)) = 0.5*n*eps0*c0*abs(VCAVsim(ll).pulse).^2;
                    VCAVsim(ll).max_I(1+VCAVsim(ll).counter,1+j) = max(0.5*n*eps0*c0*abs(VCAVsim(ll).pulse).^2);
                    VCAVsim(ll).counter=VCAVsim(ll).counter+1;
                    VCAVsim(ll).transverse_grid_short_y(VCAVsim(ll).counter)=VCAVsim(ll).transverse_grid_y(i);
                end
            end
            % Find a pulse to focus on
            [~,ind] = findpeaks(VCAVsim(ll).out_pulse(ceil(VCAVsim(ll).counter/2),:),'sortstr','descend');
            ind = sort(ind(1:min(2,length(ind))));        
            VCAVsim(ll).pulse_ind_left=find((t-t(ind(1)))/ps>-VCAVsim(ll).tmp_width_power,1);
            VCAVsim(ll).pulse_ind_right=find((t-t(ind(1)))/ps>VCAVsim(ll).tmp_width_power,1);
            peakPoint=t(ind(1));
            if (isempty(VCAVsim(ll).pulse_ind_left)||isempty(VCAVsim(ll).pulse_ind_right))
                VCAVsim(ll).pulse_ind_left=find((t-t(ind(2)))/ps>-VCAVsim(ll).tmp_width_power,1);
                VCAVsim(ll).pulse_ind_right=find((t-t(ind(2)))/ps>VCAVsim(ll).tmp_width_power,1);
                peakPoint=t(ind(2));
            end
            VCAVsim(ll).peakPoint(1+j)=peakPoint;
            VCAVsim(ll).fwhm_y(1+j) = findFWHM(VCAVsim(ll).transverse_grid_y, VCAVsim(ll).out_pulse(:,ind(1)) );
            VCAVsim(ll).fwhm_t(1+j) = findFWHM(t, VCAVsim(ll).out_pulse(ceil(VCAVsim(ll).counter/2),:) );
            VCAVsim(ll).power=zeros(VCAVsim(ll).pulse_ind_right-VCAVsim(ll).pulse_ind_left+1,1);
            if ll==1
               max_power=0;
            else
                if max(VCAVsim(ll).power)>max_power
                    max_power=max(VCAVsim(ll).power);
                end
            end
            VCAVsim(ll).t=t(VCAVsim(ll).pulse_ind_left:VCAVsim(ll).pulse_ind_right)-peakPoint;
            for jj=VCAVsim(ll).pulse_ind_left:VCAVsim(ll).pulse_ind_right
                VCAVsim(ll).power(1+jj-VCAVsim(ll).pulse_ind_left)=trapz(VCAVsim(ll).transverse_grid_short_y,VCAVsim(ll).out_pulse(:,jj));
            end
            VCAVsim(ll).energy_evol(j+1)=trapz(t(VCAVsim(ll).pulse_ind_left:VCAVsim(ll).pulse_ind_right),VCAVsim(ll).power);
            VCAVsim(ll).max_It(1+j) = mean(t);         
        end
        VCAVsim(ll).scaledInt = VCAVsim(ll).max_I*cm*cm/1e6;
    end
    
    if IO_pointDiff==1 %compute evolution of peaks difference betw consecutive sims.
        for ll=1:2:num_sims %run through all sims by pairs
            VCAVsim(ll).round_trip_time=VCAVsim(ll).round_trip_time/2;
            VCAVsim(ll).deltaT=0;%(2*5.5*1030*nm+(max(VCAVsim(ll).transverse_grid_y)-min(VCAVsim(ll).transverse_grid_y))*tan(2*pi/180))/c0; %Delay to correct for different lengths for IO_pointDiff
            for j=0:(min(VCAVsim(ll).plot_num,VCAVsim(ll).plot_num)-1) %compare all peaks up to shorter sim
               VCAVsim(ll).peakDiff(j+1) = VCAVsim(ll).peakPoint(j+1)-VCAVsim(ll+1).peakPoint(j+1);
               VCAVsim(ll).peakDiff(j+1) = VCAVsim(ll).peakDiff(j+1)-VCAVsim(ll).deltaT;
               %Adjust difference for different peak arrivals, e.g. peak 1 vs. peak 2
               if VCAVsim(ll).peakDiff(j+1)<0
                   VCAVsim(ll).peakDiff(j+1)=mod(VCAVsim(ll).peakDiff(j+1),VCAVsim(ll).round_trip_time)-VCAVsim(ll).round_trip_time;
               elseif VCAVsim(ll).peakDiff(j+1)>0
                   VCAVsim(ll).peakDiff(j+1)=mod(VCAVsim(ll).peakDiff(j+1),VCAVsim(ll).round_trip_time);
               end
            end
        end
        
        tmp_fig=figure(1002); %Plot peak intensity convergence
        set(tmp_fig,'Name','Convergence-PointDiff');
        set(tmp_fig,'Position',POS);
        hold on
        for ll=1:2:num_sims
            plot(VCAVsim(ll).max_It/ns,VCAVsim(ll).peakDiff/fs, 'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', [char(VCAVsim(ll).legend),'-',char(VCAVsim(ll+1).legend)]);
        end
        hold off
        grid on
        lgd=legend('show');
        set(lgd,'Location','best');
        %ylim([-500,500]);
        xlabel('time [ns]');
        ylabel('Arrival Time Difference [fs]');
        if IO_save==1
            saveas(tmp_fig,[saveKey,'QWconvergence-PointDiff.png']);
        end
    end

    
    tmp_fig=figure(102); %Plot peak intensity convergence
    set(tmp_fig,'Name','Convergence-PeakIntensity');
    set(tmp_fig,'Position',POS);
    hold on
    for ll=1:num_sims
        plot(VCAVsim(ll).max_It/ns,VCAVsim(ll).scaledInt(ceil(VCAVsim(ll).counter/2),:), my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(VCAVsim(ll).legend));
    end
    hold off
    grid on
    lgd=legend('show');
    set(lgd,'Location','best');
    xlabel('time [ns]');
    ylabel('Intensity [MW/cm^{2}]');
    if IO_save==1
        saveas(tmp_fig,[saveKey,'QWconvergence-PeakIntensity.png']);
    end

    tmp_fig=figure(103);
    set(tmp_fig,'Name','FWHM Time');
    set(tmp_fig,'Position',POS);
    hold on
    for ll=1:num_sims
        plot(VCAVsim(ll).max_It/ns,VCAVsim(ll).fwhm_t/fs, my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(VCAVsim(ll).legend));
    end
    hold off
    grid on
    lgd=legend('show');
    set(lgd,'Location','best');
    xlabel('time [ns]')
    ylabel('FWHM [fs]')
    if IO_save==1
        saveas(tmp_fig,[saveKey,'QWconvergence-FWHM_time.png']);
    end

    tmp_fig=figure(104);
    set(tmp_fig,'Name','FWHM Time');
    set(tmp_fig,'Position',POS);
    hold on
    for ll=1:num_sims
        plot(VCAVsim(ll).max_It/ns,VCAVsim(ll).fwhm_y/um, my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(VCAVsim(ll).legend));
    end
    hold off
    grid on
    lgd=legend('show');
    set(lgd,'Location','best');
    xlabel('time [ns]')
    ylabel('FWHM [\mum]')
    if IO_save==1
        saveas(tmp_fig,[saveKey,'QWconvergence-FWHM_space.png']);
    end

    tmp_fig=figure(105);
    set(tmp_fig,'Name','Convergence-Energy');
    set(tmp_fig,'Position',POS);
    hold on
    for ll=1:num_sims
        plot(VCAVsim(ll).max_It/ns,VCAVsim(ll).energy_evol/(10^6), my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(VCAVsim(ll).legend));
    end
    hold off
    grid on
    lgd=legend('show');
    set(lgd,'Location','best');
    xlabel('time [ns]');
    ylabel('Energy [J/cm^2]');
    grid on
    if IO_save==1
        saveas(tmp_fig,[saveKey,'QWconvergence-Energy.png']);
    end
end

%% Integrate power across transverse dimension
if IO_power==1
    disp('Plotting power overlay')
    for ll = 1:num_sims
        ll
        current_plot=VCAVsim(ll).plot_num-1;
        t = loadD([char(VCAVsim(ll).outKey),num2str(current_plot),'_t.dat']);
        VCAVsim(ll).counter=0;
        for i = 0:(VCAVsim(ll).NUM_TRANSVERSE-1)
            if abs(VCAVsim(ll).transverse_grid_y(1+i))<VCAVsim(ll).spc_width*um+VCAVsim(ll).dx %max(abs(transverse_grid_y))/8
                pulse_re = loadD([char(VCAVsim(ll).outKey),num2str(current_plot),'_E_re_',char(VCAVsim(ll).location),'_T',num2str(i),'.dat']);
                pulse_im = loadD([char(VCAVsim(ll).outKey),num2str(current_plot),'_E_im_',char(VCAVsim(ll).location),'_T',num2str(i),'.dat']);
                VCAVsim(ll).pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*VCAVsim(ll).w0);
                VCAVsim(ll).out_pulse(1+VCAVsim(ll).counter,1:length(VCAVsim(ll).pulse)) = 0.5*n*eps0*c0*abs(VCAVsim(ll).pulse).^2;
                VCAVsim(ll).counter=VCAVsim(ll).counter+1;
                VCAVsim(ll).transverse_grid_short_y(VCAVsim(ll).counter)=VCAVsim(ll).transverse_grid_y(i);
            end
        end
        % Find a pulse to focus on
        [~,ind] = findpeaks(VCAVsim(ll).out_pulse(ceil(VCAVsim(ll).counter/2),:),'sortstr','descend');
        ind = sort(ind(1:min(2,length(ind))));        
        VCAVsim(ll).pulse_ind_left=find((t-t(ind(1)))/ps>-VCAVsim(ll).tmp_width_power,1);
        VCAVsim(ll).pulse_ind_right=find((t-t(ind(1)))/ps>VCAVsim(ll).tmp_width_power,1);
        peakPoint=t(ind(1));
        if (isempty(VCAVsim(ll).pulse_ind_left)||isempty(VCAVsim(ll).pulse_ind_right))
            VCAVsim(ll).pulse_ind_left=find((t-t(ind(2)))/ps>-VCAVsim(ll).tmp_width_power,1);
            VCAVsim(ll).pulse_ind_right=find((t-t(ind(2)))/ps>VCAVsim(ll).tmp_width_power,1);
            peakPoint=t(ind(2));
        end
        VCAVsim(ll).fwhm_y(1+j) = findFWHM(VCAVsim(ll).transverse_grid_y, VCAVsim(ll).out_pulse(:,ind(1)) );
        VCAVsim(ll).fwhm_t(1+j) = findFWHM(t, VCAVsim(ll).out_pulse(ceil(VCAVsim(ll).counter/2),:) );
        VCAVsim(ll).power=zeros(VCAVsim(ll).pulse_ind_right-VCAVsim(ll).pulse_ind_left+1,1);
        VCAVsim(ll).t=t(VCAVsim(ll).pulse_ind_left:VCAVsim(ll).pulse_ind_right)-peakPoint;
        tmp=0;
        for jj=VCAVsim(ll).pulse_ind_left:VCAVsim(ll).pulse_ind_right
            tmp=tmp+1;
            VCAVsim(ll).power(tmp)=trapz(VCAVsim(ll).transverse_grid_short_y,VCAVsim(ll).out_pulse(:,jj));
        end
        if ll==1
           max_power=max(VCAVsim(ll).power);
        else
            if max(VCAVsim(ll).power)>max_power
                max_power=max(VCAVsim(ll).power);
            end
        end
        max_power;
    end
    tmp_fig=figure(222);
    set(tmp_fig,'Position',POS);
    hold on
    for ll = 1:num_sims
        plot(VCAVsim(ll).t/ps,VCAVsim(ll).power/max_power, my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(VCAVsim(ll).legend));
    end
    xlim(PowerTLIM);
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