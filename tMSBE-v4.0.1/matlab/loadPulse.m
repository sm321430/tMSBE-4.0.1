function [t_cur,points,pulse,out_pulse]=loadPulse(outKey,plot_num,num_peaks,location)
    setupConstants
    
    t = loadD([outKey,num2str(plot_num-1),'_t.dat']);
    w0 = loadD([outKey,'w0.dat']);
    pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T0.dat']);
    pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_',location,'_T0.dat']);
    pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*w0);
    out_pulse = 0.5*n*eps0*c0*abs(pulse).^2;
    t_cur=(t-t(1))/ps; %Current time refernced from t(1) in picoseconds
    ind_trans=ones(1,20); %indices for peaks
    peaks_trans=zeros(1,20); %Values for peaks
    points=zeros(num_peaks,1); %time values for peaksinitialized
    
    % Find top peaks, away from edges for all transverse points
    [peaks,time_ind]=findpeaks(out_pulse,'SortStr','descend');
    ii=min(50,length(peaks));
    ind_trans(1,1:ii)=time_ind(1:ii);
    peaks_trans(1,1:ii)=peaks(1:ii);
    % Out of found peaks, find max separated peaks
    jk=1;
    out_pulse_tmp=out_pulse;
    while jk<=num_peaks
        col=find(out_pulse_tmp==max(max(peaks_trans)),1);
        if (min(abs(points-t_cur(col)))>1)||((t_cur(col)<1) && max(points(points<1))==0)
            points(jk)=t_cur(col);
            jk=jk+1;
        end
        out_pulse_tmp(col)=0;
        col=find(peaks_trans==max(max(peaks_trans)),1);
        peaks_trans(col)=0;
    end
    points=sort(points);
end