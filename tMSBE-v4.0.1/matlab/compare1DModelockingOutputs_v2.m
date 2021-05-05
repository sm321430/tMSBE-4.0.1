%%Comparisons for outputs/evolution of 1D MSBE modelocking simulations
%S.A.Mclaren (December 03, 2020)
close all
clear all

%IOs for each plot located here
IO_CWandCCW=0; %doubles up simulation numbers and includes CCW direction
IO_convergence=0; %Plot convergence of FWHM, peak intensity, and energy
IO_pointDiff=0; %Plot convergence difference in time between peaks. Requires IO_convergence. Not configured
IO_spectrums=0; % Plot spectrums of particular peaks. Not currently configure
%IO_temp_dens=0; %Plot density and temperature (Not currently configured)
IO_finalOutput=1; %Plot final output
IO_finalOutputPlot_roundTripTime=0; %if IO_finalOutput==1 plots line for round trip time for each sim.
%IO_finalOutput_indPulse=1; %Individual plots of pulses at final output (Not currently configured)
IO_power=1; %Plot power overlay for latest output
%IO_transIntensity=0; %Plot transverse intensity profiles (Not currently configured)
%IO_pointDiff=0; %Plot time differences between peaks (Not currently configured)
%IO_findFWHM=0; %Find FWHM for each peak (saved in pulseFWHM) (Not currently configured)
IO_phaseSpace=0; %Plot converged values in parameter space given by labeled values and axes 
IO_save=1; %Save figures
IO_separatePeaks=1; %Ensure peaks are distinct by set amount

%Setup
setupPlot
setupConstants
setupSave

width=[20.0,20.0]; %Time on either side of peak in ps
PowerTLIM=[-0.25,0.25]; %Power plot  temporal width
PowerTTICK=[-0.25,0.0,0.25]; %Power plot tick marks
spectrum_width=[1.17,1.24];
sizeIl=500*fs;
sizeIr=2000*fs;


%% All load locations
setupComparisons
if IO_CWandCCW==1
    for j=num_sims:-1:1
        CAVsim(2*j)=CAVsim(j);
        CAVsim(2*j-1)=CAVsim(j);
        CAVsim(2*j).location='OUTPUTBACK';
        CAVsim(2*j).legend='CCW';
        j
    end
    num_sims=num_sims*2;
end

round_trip_time_max=0;
for j = 1 : num_sims
   tmp_cnt=1; %temporary counter for plot num
   while isfile([char(CAVsim(j).outKey),num2str(tmp_cnt-1),'_E_re_',char(CAVsim(j).location),'_T0.dat'])
          tmp_cnt=tmp_cnt+1; %Check for next plot. If exists, add to counter
   end
   
   CAVsim(j).plot_num=tmp_cnt-1; %number of output plots for convergence
   %CAVsim(j).plot_point=2; %which peak to plot
   CAVsim(j).round_trip_time = loadD([char(CAVsim(j).outKey),'round_trip_time.dat'])/(1+IO_CWandCCW);
   %CAVsim(j).round_trip_time = CAVsim(j).round_trip_time/2;
   CAVsim(j).tmp_width_power=CAVsim(j).round_trip_time/(3*ps); %max(1.25*PowerTLIM);% %0.5; %F%0.1*CAVsim(j).round_trip_time/(ps);
   CAVsim(j).w0 = loadD([char(CAVsim(j).outKey),'w0.dat']);
   CAVsim(j).OC=0.01;
   if CAVsim(j).round_trip_time/ps>round_trip_time_max
        round_trip_time_max=CAVsim(j).round_trip_time/ps;
   end
end 

if IO_CWandCCW==1
    lineStyle_mod=2;
    lineStyle_shift=1;
else
    lineStyle_mod=1;
    lineStyle_shift=2;
end

if IO_convergence==1
    % Load and find peaks for all outKeys
    for ll = 1:num_sims
        ll
        for j = 0:(CAVsim(ll).plot_num-1)
            j
            current_plot=j;
            t = loadD([char(CAVsim(ll).outKey),num2str(j),'_t.dat']);
            CAVsim(ll).counter=0;
            pulse_re = loadD([char(CAVsim(ll).outKey),num2str(j),'_E_re_',char(CAVsim(ll).location),'_T',num2str(0),'.dat']);
            pulse_im = loadD([char(CAVsim(ll).outKey),num2str(j),'_E_im_',char(CAVsim(ll).location),'_T',num2str(0),'.dat']);
            CAVsim(ll).pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*CAVsim(ll).w0);
            CAVsim(ll).out_pulse = 0.5*n*eps0*c0*abs(CAVsim(ll).pulse).^2;
            CAVsim(ll).max_I(1,1+j) = max(CAVsim(ll).out_pulse);
            CAVsim(ll).round_trip_time_mod=CAVsim(ll).round_trip_time/2;
           
            % Find a pulse to focus on
            [~,ind] = findpeaks(CAVsim(ll).out_pulse,'sortstr','descend');
            
            CAVsim(ll).pulse_ind_left=find((t-t(ind(1)))/ps>-CAVsim(ll).tmp_width_power,1);
            CAVsim(ll).pulse_ind_right=find((t-t(ind(1)))/ps>CAVsim(ll).tmp_width_power,1);
            peakPoint=t(ind(1));
            
            if (isempty(CAVsim(ll).pulse_ind_left)||isempty(CAVsim(ll).pulse_ind_right)||CAVsim(ll).pulse_ind_left==1)
                CAVsim(ll).pulse_ind_left=find((t-t(ind(2)))/ps>-CAVsim(ll).tmp_width_power,1);
                CAVsim(ll).pulse_ind_right=find((t-t(ind(2)))/ps>CAVsim(ll).tmp_width_power,1);
                peakPoint=t(ind(2));
            end
            
            CAVsim(ll).fwhm_t(1+j) = findFWHM(t, CAVsim(ll).out_pulse);
            CAVsim(ll).power=CAVsim(ll).out_pulse(CAVsim(ll).pulse_ind_left:CAVsim(ll).pulse_ind_right);
            CAVsim(ll).peakPoint(1+j)=peakPoint;
            
            ind=sort(ind(1:floor((t(end)-t(1))/CAVsim(ll).round_trip_time_mod)));
            CAVsim(ll).peakPoint_first(j+1)=t(ind(1));
            CAVsim(ll).peakPoint_second(j+1)=t(ind(2));
            if ll==1
               max_power=0;
            else
                if max(CAVsim(ll).power)>max_power
                    max_power=max(CAVsim(ll).power);
                end
            end
            CAVsim(ll).t=t(CAVsim(ll).pulse_ind_left:CAVsim(ll).pulse_ind_right)-peakPoint;
            CAVsim(ll).energy_evol(j+1)=trapz(t(CAVsim(ll).pulse_ind_left:CAVsim(ll).pulse_ind_right),CAVsim(ll).power);
            %CAVsim(ll).energy_evol(j+1)=trapz(t,CAVsim(ll).out_pulse)/((t(end)-t(1))/(CAVsim(ll).round_trip_time/2));
            CAVsim(ll).max_It(1+j) = mean(t);         
        end
        CAVsim(ll).scaledInt = CAVsim(ll).max_I*cm*cm/1e6;
    end
    
    if IO_pointDiff==1 && IO_CWandCCW==1 %compute evolution of peaks difference betw consecutive sims.
        for ll=1:2:num_sims %run through all sims by pairs
            CAVsim(ll).deltaT=0;%(2*5.5*1030*nm+(max(CAVsim(ll).transverse_grid_y)-min(CAVsim(ll).transverse_grid_y))*tan(2*pi/180))/c0; %Delay to correct for different lengths for IO_pointDiff
            for j=0:(min(CAVsim(ll).plot_num,CAVsim(ll).plot_num)-1) %compare all peaks up to shorter sim
                %Ensure difference is always taken with leading CCW pulse for consistancy
                CAVsim(ll).peakDiff(j+1)=CAVsim(ll).peakPoint_first(j+1)-CAVsim(ll+1).peakPoint_first(j+1);
                if CAVsim(ll).peakDiff(j+1) < 0
                   CAVsim(ll).peakDiff(j+1)=CAVsim(ll).peakPoint_second(j+1)-CAVsim(ll+1).peakPoint_first(j+1); 
                end
            end
        end
        
        tmp_fig=figure(1002); %Plot peak intensity convergence
        set(tmp_fig,'Name','Convergence-PointDiff');
        set(tmp_fig,'Position',POS);
        hold on
        for ll=1:2:num_sims
            plot(CAVsim(ll).max_It/ns,CAVsim(ll).peakDiff/fs, 'DisplayName', [char(CAVsim(ll).legend),'-',char(CAVsim(ll+1).legend)]);
        end
        hold off
        grid on
        lgd=legend('show');
        set(lgd,'Location','best');
        ylim([-1000,1000]);
        xlabel('time [ns]');
        ylabel('Arrival Time Difference [fs]');
        if IO_save==1
            saveas(tmp_fig,[saveKey,'QWconvergence-PointDiff.png']);
        end
    end
    
    tmp_fig=figure(102);
    set(tmp_fig,'Name','Convergence-PeakIntensity');
    set(tmp_fig,'Position',POS);
    hold on
    for ll=1:num_sims
        plot(CAVsim(ll).max_It/ns,CAVsim(ll).scaledInt(1,:), my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(CAVsim(ll).legend));
    end
    hold off
    grid on
    lgd=legend('show');
    if(num_sims>=6)
        lgd.FontSize=30;
    end
    xlim([0,250]);
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
        plot(CAVsim(ll).max_It/ns,CAVsim(ll).fwhm_t/fs,  my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(CAVsim(ll).legend));
    end
    hold off
    grid on
    lgd=legend('show');
    if(num_sims>=6)
        lgd.FontSize=30;
    end
    set(lgd,'Location','best');
    xlabel('time [ns]')
    ylabel('FWHM [fs]')
    if IO_save==1
        saveas(tmp_fig,[saveKey,'QWconvergence-FWHM_time.png']);
    end
    
    tmp_fig=figure(105);
    set(tmp_fig,'Name','Convergence-Energy');
    set(tmp_fig,'Position',POS);
    hold on
    for ll=1:num_sims
        plot(CAVsim(ll).max_It/ns,CAVsim(ll).energy_evol/(10^6), my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(CAVsim(ll).legend));
    end
    hold off
    grid on
    lgd=legend('show');
    if(num_sims>=6)
        lgd.FontSize=30;
    end
    set(lgd,'Location','best');
    xlabel('time [ns]');
    ylabel('Energy [J/cm^2]');
    grid on
    if IO_save==1
        saveas(tmp_fig,[saveKey,'QWconvergence-Energy.png']);
    end
end

%% Plot final output
if IO_finalOutput==1
    disp('Plotting final output');
    for ll=1:num_sims
        ll
        current_plot=CAVsim(ll).plot_num-1;
        t = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_t.dat']);
        CAVsim(ll).t=t;
        pulse_re = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_E_re_',char(CAVsim(ll).location),'_T',num2str(0),'.dat']);
        pulse_im = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_E_im_',char(CAVsim(ll).location),'_T',num2str(0),'.dat']);
        
        CAVsim(ll).pulse = (pulse_re + 1i*pulse_im).*exp(-1i*CAVsim(ll).t*CAVsim(ll).w0);
        CAVsim(ll).out_pulse = abs(CAVsim(ll).pulse);%WTFF 0.5*n*eps0*c0*abs(CAVsim(ll).pulse).^2;
        % Find a pulse to focus on
        [~,ind] = findpeaks(CAVsim(ll).out_pulse,'sortstr','descend');
        ind = sort(ind(1:min(2,length(ind))));        
        if(~isempty(ind))
            CAVsim(ll).pulse_ind_left=find((t-t(ind(1)))/ps>-CAVsim(ll).tmp_width_power,1);
            CAVsim(ll).pulse_ind_right=find((t-t(ind(1)))/ps>CAVsim(ll).tmp_width_power,1);
            peakPoint=t(ind(1));
            if (isempty(CAVsim(ll).pulse_ind_left)||isempty(CAVsim(ll).pulse_ind_right)||CAVsim(ll).pulse_ind_left==1)
                CAVsim(ll).pulse_ind_left=find((t-t(ind(2)))/ps>-CAVsim(ll).tmp_width_power,1);
                CAVsim(ll).pulse_ind_right=find((t-t(ind(2)))/ps>CAVsim(ll).tmp_width_power,1);
                peakPoint=t(ind(2));
            end
        else
            CAVsim(ll).pulse_ind_left=1;
            CAVsim(ll).pulse_ind_right=length(t);
            peakPoint=1;
        end
        CAVsim(ll).peakPoint=peakPoint;
        CAVsim(ll).t=CAVsim(ll).t-peakPoint;
        %CAVsim(ll).t=CAVsim(ll).t-CAVsim(ll).t(1);
    end
    delay=zeros(num_sims,1); %variable for artificially shifting data
    
    tmp_fig=figure(22);
    set(tmp_fig,'Name','Final Output Plot');
    set(tmp_fig,'Position',POS);
    hold on
    %delay=[-4.37*ps,-4.37*ps,0,0];
    for ll=1:num_sims
        plot((CAVsim(ll).t+delay(ll))/ps, (CAVsim(ll).out_pulse*cm*cm/1e6), my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(CAVsim(ll).legend));
        if IO_finalOutputPlot_roundTripTime==1
            go=plot(CAVsim(ll).round_trip_time*[1,1]/ps,[min(CAVsim(ll).out_pulse*cm*cm/1e6),max(CAVsim(ll).out_pulse*cm*cm/1e6)],my_plot_style{mod(ll,lineStyle_mod)+2},'color', my_plot_color{ceil(ll/lineStyle_mod)});
            go.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
    hold off
    grid on
    %xlim([-5*round_trip_time_max,5*round_trip_time_max*1.1]);
    xlim([-1,round_trip_time_max*0.6]);
    %xlim([-3,200]);
    lgd=legend('show');
    if(num_sims>=6)
        lgd.FontSize=30;
    end
    set(lgd,'Location','best');  
    xlabel('time [ps]');
    ylabel('Intensity [MW/cm^2]');
    grid on
    if IO_save==1
        saveas(tmp_fig,[saveKey,'finalOutputPlot.png']);
    end
  
    tmp_fig=figure(222);
    set(tmp_fig,'Name','Final Output Plot (Log)');
    set(tmp_fig,'Position',POS);
    %delay=[-4.37*ps,-4.37*ps,0,0];
    hold on
    for ll=1:num_sims
        plot((CAVsim(ll).t+delay(ll))/ps, log(CAVsim(ll).out_pulse*cm*cm/1e6), my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(CAVsim(ll).legend));
        if IO_finalOutputPlot_roundTripTime==1
            go=plot(CAVsim(ll).round_trip_time*[1,1]/ps,[min(log(CAVsim(ll).out_pulse*cm*cm/1e6)),max(log(CAVsim(ll).out_pulse*cm*cm/1e6))],my_plot_style{mod(ll,lineStyle_mod)+2},'color', my_plot_color{ceil(ll/lineStyle_mod)});
            go.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
    hold off
    grid on
    %xlim([-5*round_trip_time_max,5*round_trip_time_max*1.1]);
    xlim([-1,round_trip_time_max*0.6]);
    %xlim([-3,40]);
    lgd=legend('show');
    if(num_sims>=6)
        lgd.FontSize=30;
    end
    set(lgd,'Location','best');  
    xlabel('time [ps]');
    ylabel('log(Int.) [MW/cm^2]');
    %ylim([-40,1]);
    grid on
    if IO_save==1
        saveas(tmp_fig,[saveKey,'finalOutputPlotLog.png']);
    end
end

%% Plot spectrums
if IO_spectrums==1
    disp('Plotting spectrums');
    for ll = 1:num_sims
        ll
        current_plot=CAVsim(ll).plot_num-1;
        t = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_t.dat']);
        t=t-t(1);
        CAVsim(ll).counter=0;
        pulse_re = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_E_re_',char(CAVsim(ll).location),'_T',num2str(0),'.dat']);
        pulse_im = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_E_im_',char(CAVsim(ll).location),'_T',num2str(0),'.dat']);
        CAVsim(ll).pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*CAVsim(ll).w0);
        CAVsim(ll).out_pulse = 0.5*n*eps0*c0*abs(CAVsim(ll).pulse).^2;
        % Find a pulse to focus on
        [~,ind] = findpeaks(CAVsim(ll).out_pulse,'sortstr','descend');
        ind = sort(ind(1:min(2,length(ind))));  
        CAVsim(ll).pulse_ind_left=find((t-t(ind(1)))>-sizeIl,1);
        CAVsim(ll).pulse_ind_right=find((t-t(ind(1)))>sizeIr,1);
        peakPoint=t(ind(1));
        if (isempty(CAVsim(ll).pulse_ind_left)||isempty(CAVsim(ll).pulse_ind_right))
            CAVsim(ll).pulse_ind_left=find((t-t(ind(2)))>-sizeIl,1);
            CAVsim(ll).pulse_ind_right=find((t-t(ind(2)))>sizeIr,1);
            peakPoint=t(ind(2));
        end
        
        [w_qw, Y1_qw, Y2_qw] = getSpectrums_fromPoints(t/ps,sizeIl,sizeIr,t(ind(1:2))/ps,transpose(CAVsim(ll).pulse));
        CAVsim(ll).pulse_spectrum=Y1_qw;
        CAVsim(ll).w_qw=w_qw;
    end
  
    tmp_fig=figure(2009);
    set(tmp_fig,'Name','Pulse spectrums');
    set(tmp_fig,'Position',POS);
    hold on
    for ll = 1:num_sims
        plot(hbar*CAVsim(ll).w_qw/e,abs(CAVsim(ll).pulse_spectrum), my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(CAVsim(ll).legend));
    end
    hold off
    grid on
    xlim(spectrum_width);
    lgd=legend('show');
    if(num_sims>=6)
        lgd.FontSize=30;
    end
    xlabel('Energy [eV]');
    ylabel('Spectrum [a.u.]');
    if IO_save==1
        saveas(tmp_fig,[saveKey,'QWspectrumOverlay.png']);
    end 
end

%% Integrate power across transverse dimension
if IO_power==1
    disp('Plotting power overlay')
    for ll = 1:num_sims
        ll
        current_plot=CAVsim(ll).plot_num-1;
        t = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_t.dat']);
        CAVsim(ll).counter=0;
        pulse_re = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_E_re_',char(CAVsim(ll).location),'_T',num2str(0),'.dat']);
        pulse_im = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_E_im_',char(CAVsim(ll).location),'_T',num2str(0),'.dat']);
        CAVsim(ll).pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*CAVsim(ll).w0);
        CAVsim(ll).out_pulse = 0.5*n*eps0*c0*abs(CAVsim(ll).pulse).^2;
        % Find a pulse to focus on
        [~,ind] = findpeaks(CAVsim(ll).out_pulse,'sortstr','descend');
        ind = sort(ind(1:min(2,length(ind))));        
        CAVsim(ll).pulse_ind_left=find((t-t(ind(1)))/ps>-CAVsim(ll).tmp_width_power,1);
        CAVsim(ll).pulse_ind_right=find((t-t(ind(1)))/ps>CAVsim(ll).tmp_width_power,1);
        peakPoint=t(ind(1));
        if (isempty(CAVsim(ll).pulse_ind_left)||isempty(CAVsim(ll).pulse_ind_right)||CAVsim(ll).pulse_ind_left==1)
            CAVsim(ll).pulse_ind_left=find((t-t(ind(2)))/ps>-CAVsim(ll).tmp_width_power,1);
            CAVsim(ll).pulse_ind_right=find((t-t(ind(2)))/ps>CAVsim(ll).tmp_width_power,1);
            peakPoint=t(ind(2));
        end
        CAVsim(ll).fwhm_t(1+j) = findFWHM(t, CAVsim(ll).out_pulse );
        CAVsim(ll).power=CAVsim(ll).out_pulse(CAVsim(ll).pulse_ind_left:CAVsim(ll).pulse_ind_right);
        CAVsim(ll).t=t(CAVsim(ll).pulse_ind_left:CAVsim(ll).pulse_ind_right)-peakPoint;
        CAVsim(ll).peakPoint=peakPoint;
        if ll==1
           max_power=max(CAVsim(ll).power);
        else
            if max(CAVsim(ll).power)>max_power
                max_power=max(CAVsim(ll).power);
            end
        end
    end
    tmp_fig=figure(221);
    set(tmp_fig,'Position',POS);
    hold on
    for ll = 1:num_sims
        plot(CAVsim(ll).t/ps,CAVsim(ll).power*cm*cm/1e6, my_plot_style{mod(ll,lineStyle_mod)+lineStyle_shift},'color', my_plot_color{ceil(ll/lineStyle_mod)}, 'DisplayName', char(CAVsim(ll).legend));
    end
    xlim(PowerTLIM);
    %set(gca,'YTick', [0,0.25,0.5,0.75,1.0]);
    set(gca,'XTick', PowerTTICK);
    xlabel('time [ps]');
    ylabel('Intensity [MW/cm^2]');
    hold off
    grid on
    lgd=legend('show');
    if(num_sims>=6)
        lgd.FontSize=30;
    end
    if IO_save==1
        saveas(tmp_fig,[saveKey,'QWpowerOverlay.png']);
    end
end

%% Plot latest results in phase space
if IO_phaseSpace==1
    if CAVsim(1).cavType=='RCAV'
        location={'OUTPUT','OUTPUTBACK'};
            for ll=1:num_sims  
            ll
                for kk=1:2 %CW,CCW
                  %Load latest output
                  %Pulse or no pulse
                  %find peak intensity
                  %find FWHM time
                current_plot=CAVsim(ll).plot_num-1;
                t = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_t.dat']);
                CAVsim(ll).counter=0;
                pulse_re = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_E_re_',char(location{kk}),'_T',num2str(0),'.dat']);
                pulse_im = loadD([char(CAVsim(ll).outKey),num2str(current_plot),'_E_im_',char(location{kk}),'_T',num2str(0),'.dat']);
                CAVsim(ll).pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*CAVsim(ll).w0);
                CAVsim(ll).out_pulse = 0.5*n*eps0*c0*abs(CAVsim(ll).pulse).^2;
                CAVsim(ll).max_I(kk) = max(CAVsim(ll).out_pulse);

                % Find a pulse to focus on
                [~,ind] = findpeaks(CAVsim(ll).out_pulse,'sortstr','descend');
                ind = sort(ind(1:min(2,length(ind))));        
                CAVsim(ll).pulse_ind_left=find((t-t(ind(1)))/ps>-CAVsim(ll).tmp_width_power,1);
                CAVsim(ll).pulse_ind_right=find((t-t(ind(1)))/ps>CAVsim(ll).tmp_width_power,1);
                peakPoint=t(ind(1));
                if (isempty(CAVsim(ll).pulse_ind_left)||isempty(CAVsim(ll).pulse_ind_right))
                    CAVsim(ll).pulse_ind_left=find((t-t(ind(2)))/ps>-CAVsim(ll).tmp_width_power,1);
                    CAVsim(ll).pulse_ind_right=find((t-t(ind(2)))/ps>CAVsim(ll).tmp_width_power,1);
                    peakPoint=t(ind(2));
                end
                CAVsim(ll).fwhm_t(kk) = findFWHM(t, CAVsim(ll).out_pulse);
                CAVsim(ll).power=CAVsim(ll).out_pulse(CAVsim(ll).pulse_ind_left:CAVsim(ll).pulse_ind_right);
                CAVsim(ll).peakPoint(1+j)=peakPoint;
                if ll==1
                   max_power=0;
                else
                    if max(CAVsim(ll).power)>max_power
                        max_power=max(CAVsim(ll).power);
                    end
                end
                CAVsim(ll).t=t(CAVsim(ll).pulse_ind_left:CAVsim(ll).pulse_ind_right)-peakPoint;
                CAVsim(ll).energy(kk)=trapz(t(CAVsim(ll).pulse_ind_left:CAVsim(ll).pulse_ind_right),CAVsim(ll).power);        
                CAVsim(ll).scaledInt(kk) = CAVsim(ll).max_I(kk)*cm*cm/1e6;
                CAVsim(ll).meanInt(kk)=mean(CAVsim(ll).out_pulse);

                    if CAVsim(ll).fwhm_t(kk)>0 && CAVsim(ll).scaledInt(kk)>0.75
                        CAVsim(ll).isPulse(kk)=1;
                    else
                        CAVsim(ll).isPulse(kk)=0;
                    end
                end
                phaseSpacePosX(ll)=CAVsim(ll).phaseSpacePos(1);
                phaseSpacePosY(ll)=CAVsim(ll).phaseSpacePos(2);
               
                CW_isPulse(ll)=CAVsim(ll).isPulse(1);
                if CW_isPulse(ll)
                    CW_fwhm_t(ll)=CAVsim(ll).fwhm_t(1);
                    CW_scaledInt(ll)=CAVsim(ll).scaledInt(1);
                else
                    CW_fwhm_t(ll)=0.0;
                    CW_scaledInt(ll)=0.0;
                end
                
                CCW_isPulse(ll)=CAVsim(ll).isPulse(2);
                if CCW_isPulse(ll)
                    CCW_fwhm_t(ll)=CAVsim(ll).fwhm_t(2);
                    CCW_scaledInt(ll)=CAVsim(ll).scaledInt(2);
                else
                    CCW_fwhm_t(ll)=0.0;
                    CCW_scaledInt(ll)=0.0;
                end
            end
            dt = delaunayTriangulation(transpose(phaseSpacePosX),transpose(phaseSpacePosY));
            tri = dt.ConnectivityList; 
            
            tmp_fig=figure('name','CW FWHM');
            trisurf(tri,phaseSpacePosX,phaseSpacePosY,CW_fwhm_t/fs);
            xlabel(CAVsim(1).phaseSpaceLabel(1));
            ylabel(CAVsim(1).phaseSpaceLabel(2));
            xlim([min(phaseSpacePosX),max(phaseSpacePosX)]);
            ylim([min(phaseSpacePosY),max(phaseSpacePosY)]);
            view(0,90);
            colorbar;
            if IO_save==1
                saveas(tmp_fig,[saveKey,'CWFWHM.png']);
            end

            tmp_fig=figure('name','CCW FWHM');
            trisurf(tri,phaseSpacePosX,phaseSpacePosY,CCW_fwhm_t/fs);  
            xlabel(CAVsim(1).phaseSpaceLabel(1));
            ylabel(CAVsim(1).phaseSpaceLabel(2));
            xlim([min(phaseSpacePosX),max(phaseSpacePosX)]);
            ylim([min(phaseSpacePosY),max(phaseSpacePosY)]);
            view(0,90);
            colorbar; 
            if IO_save==1
                saveas(tmp_fig,[saveKey,'CCWFWHM.png']);
            end            
            
            tmp_fig=figure('name','CW max intensity');
            trisurf(tri,phaseSpacePosX,phaseSpacePosY,CW_scaledInt);
            xlabel(CAVsim(1).phaseSpaceLabel(1));
            ylabel(CAVsim(1).phaseSpaceLabel(2));
            xlim([min(phaseSpacePosX),max(phaseSpacePosX)]);
            ylim([min(phaseSpacePosY),max(phaseSpacePosY)]);
            view(0,90);
            colorbar;
            if IO_save==1
                saveas(tmp_fig,[saveKey,'CWMaxIntensity.png']);
            end
            
            tmp_fig=figure('name','CCW max intensity');
            trisurf(tri,phaseSpacePosX,phaseSpacePosY,CCW_scaledInt);
            xlabel(CAVsim(1).phaseSpaceLabel(1));
            ylabel(CAVsim(1).phaseSpaceLabel(2));
            xlim([min(phaseSpacePosX),max(phaseSpacePosX)]);
            ylim([min(phaseSpacePosY),max(phaseSpacePosY)]);            
            view(0,90);
            colorbar;
            if IO_save==1
                saveas(tmp_fig,[saveKey,'CCWMaxIntensity.png']);
            end           
            
            tmp_fig=figure('name','CW pulse existence');
            trisurf(tri,phaseSpacePosX,phaseSpacePosY,CW_isPulse);
            xlabel(CAVsim(1).phaseSpaceLabel(1));
            ylabel(CAVsim(1).phaseSpaceLabel(2));
            xlim([min(phaseSpacePosX),max(phaseSpacePosX)]);
            ylim([min(phaseSpacePosY),max(phaseSpacePosY)]);
            view(0,90);
            colorbar;            
            if IO_save==1
                saveas(tmp_fig,[saveKey,'CWPulseExist.png']);
            end            
            
            tmp_fig=figure('name','CCW pulse existence');
            trisurf(tri,phaseSpacePosX,phaseSpacePosY,CCW_isPulse);
            xlabel(CAVsim(1).phaseSpaceLabel(1));
            ylabel(CAVsim(1).phaseSpaceLabel(2));
            xlim([min(phaseSpacePosX),max(phaseSpacePosX)]);
            ylim([min(phaseSpacePosY),max(phaseSpacePosY)]);
            view(0,90);
            colorbar;
            if IO_save==1
                saveas(tmp_fig,[saveKey,'CCWPulseExist.png']);
            end            
            
            
            tmp_fig=figure('name','Pulse number');
            trisurf(tri,phaseSpacePosX,phaseSpacePosY,CW_isPulse+CCW_isPulse);
            xlabel(CAVsim(1).phaseSpaceLabel(1));
            ylabel(CAVsim(1).phaseSpaceLabel(2));
            xlim([min(phaseSpacePosX),max(phaseSpacePosX)]);
            ylim([min(phaseSpacePosY),max(phaseSpacePosY)]);
            view(0,90);
            colorbar;
            if IO_save==1
                saveas(tmp_fig,[saveKey,'PulseNumber.png']);
            end            
            
    else
        disp("IO_PhaseSpace only setup for ring cav currently");
        asd
    end
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
if i1==length(pulse)
   i1=i1-1; 
end

while ((i1 < (length(pulse)-1))&&(pulse(i1+1) > 0.5*maxp))
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