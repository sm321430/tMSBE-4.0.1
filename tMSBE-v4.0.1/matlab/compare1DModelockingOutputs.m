%%Comparisons for outputs of modelocking simulations
%S.A.Mclaren (June 2020)

clear all
close all

IO_finalPlotLog=1;
IO_finalPlot=1;
IO_spectrumsPlot=1;
IO_spectrumsPlot2=1;
IO_point_diff=0;
IO_dens=0;
IO_save=1;

setupPlot
setupConstants

plot_num=[1,1,1,1,1]; %Output number +1
num_peaks=2;
width=[2.5,2.5];
plot_point=[2,2,2,2,2,2];
width_short=[-4.0,10.0];
width_spectrum=[2.0,2.0];
YLIM_intensity=[10^(-13);10^(1)];
XLIM_SPECTRUM=[1.1,1.3];
location='CAVOC'; %Field location for uploading and saving
date='081020';
test='VCAVXX-redone2';
test_folder='test';
saveKey_local='Fall2019-Spring2020/VCAV/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey = [saveKey_fold,'/',test,'-out',num2str(plot_num(1)-1),'-',location,'-'];
user_entry = input(['saveKey=',saveKey,'? (y to continue): '], 's');
if user_entry~='y'
    gafuggle
end
if exist(saveKey_fold,'dir')~=7
    mkdir(saveKey_fold);
end

outKey = '../../tMSBE-v3.7-threshold1/run/out__'; 
[t_cur,points,pulse,out_pulse]=loadPulse(outKey,plot_num(1),num_peaks,location1);
[w_qw_dens1, spectrum_dens1, spectrum_dens1_2] = getSpectrums_fromPoints(t_cur,width_spectrum(1),width_spectrum(2),points,pulse);
t_cur_ind_left=find(t_cur>(points(plot_point(1))-width(1)),1);
t_cur_ind_right=find(t_cur>(points(plot_point(1))+width(2)),1);
if isempty(t_cur_ind_left)
   t_cur_ind_left=1; 
end
if isempty(t_cur_ind_right)
   t_cur_ind_right=length(t_cur); 
end
t_cur_dens1=t_cur(t_cur_ind_left:t_cur_ind_right);
t_cur_dens1=t_cur_dens1-t_cur_dens1(1)-width(1);
pulse_dens1=pulse(t_cur_ind_left:t_cur_ind_right);
out_pulse_dens1=out_pulse(t_cur_ind_left:t_cur_ind_right);
%points_diff_dens1=points(2:end)-points(1:end-1);
%Ne = loadD([outKey,num2str(plot_point(1)-1),'_Nsum_e_',location,'_T1.dat']);
%Te = loadD([outKey,num2str(plot_point(1)-1),'_inst_temp_e_',location,'_T1.dat']);
%Ne_dens1=Ne(t_cur_ind_left:t_cur_ind_right);
%Te_dens1=Te(t_cur_ind_left:t_cur_ind_right);

outKey = '../../tMSBE-v3.7-threshold2/run/out__';
[t_cur,points,pulse,out_pulse]=loadPulse(outKey,plot_num(2),num_peaks,location2);
[w_qw_dens2, spectrum_dens2, spectrum_dens2_2] = getSpectrums_fromPoints(t_cur,width_spectrum(1),width_spectrum(2),points,pulse);
t_cur_ind_left=find(t_cur>(points(plot_point(2))-width(1)),1);
t_cur_ind_right=find(t_cur>(points(plot_point(2))+width(2)),1);
if isempty(t_cur_ind_left)
   t_cur_ind_left=1; 
end
if isempty(t_cur_ind_right)
   t_cur_ind_right=length(t_cur); 
end
t_cur_dens2=t_cur(t_cur_ind_left:t_cur_ind_right);
t_cur_dens2=t_cur_dens2-t_cur_dens2(1)-width(1);
pulse_dens2=pulse(t_cur_ind_left:t_cur_ind_right);
out_pulse_dens2=out_pulse(t_cur_ind_left:t_cur_ind_right);
%points_diff_dens2=points(2:end)-points(1:end-1);
%Ne = loadD([outKey,num2str(plot_point(2)-1),'_Nsum_e_',location,'_T1.dat']);
%Te = loadD([outKey,num2str(plot_point(2)-1),'_inst_temp_e_',location,'_T1.dat']);
%Ne_dens2=Ne(t_cur_ind_left:t_cur_ind_right);
%Te_dens2=Te(t_cur_ind_left:t_cur_ind_right);
% 
outKey = '../../tMSBE-v3.7-threshold3/run/out__';
[t_cur,points,pulse,out_pulse]=loadPulse(outKey,plot_num(3),num_peaks,location);
[w_qw_dens3, spectrum_dens3, spectrum_dens3_2] = getSpectrums_fromPoints(t_cur,width_spectrum(1),width_spectrum(2),points,pulse);
t_cur_ind_left=find(t_cur>(points(plot_point(3))-width(1)),1);
t_cur_ind_right=find(t_cur>(points(plot_point(3))+width(2)),1);
if isempty(t_cur_ind_left)
   t_cur_ind_left=1; 
end
if isempty(t_cur_ind_right)
   t_cur_ind_right=length(t_cur); 
end
t_cur_dens3=t_cur(t_cur_ind_left:t_cur_ind_right);
t_cur_dens3=t_cur_dens3-t_cur_dens3(1)-width(1);
pulse_dens3=pulse(t_cur_ind_left:t_cur_ind_right);
out_pulse_dens3=out_pulse(t_cur_ind_left:t_cur_ind_right);
%points_diff_dens3=points(2:end)-points(1:end-1);
% %Ne = loadD([outKey,num2str(plot_point(3)-1),'_Nsum_e_',location,'_T1.dat']);
% %Te = loadD([outKey,num2str(plot_point(3)-1),'_inst_temp_e_',location,'_T1.dat']);
% %Ne_dens3=Ne(t_cur_ind_left:t_cur_ind_right);
% %Te_dens3=Te(t_cur_ind_left:t_cur_ind_right);

outKey = '../../tMSBE-v3.7-threshold4/run/out__';
[t_cur,points,pulse,out_pulse]=loadPulse(outKey,plot_num(4),num_peaks,location);
[w_qw_dens4, spectrum_dens4, spectrum_dens4_2] = getSpectrums_fromPoints(t_cur,width_spectrum(1),width_spectrum(2),points,pulse);
t_cur_ind_left=find(t_cur>(points(plot_point(4))-width(1)),1);
t_cur_ind_right=find(t_cur>(points(plot_point(4))+width(2)),1);
if isempty(t_cur_ind_left)
   t_cur_ind_left=1; 
end
if isempty(t_cur_ind_right)
   t_cur_ind_right=length(t_cur); 
end
t_cur_dens4=t_cur(t_cur_ind_left:t_cur_ind_right);
t_cur_dens4=t_cur_dens4-t_cur_dens4(1)-width(1);
pulse_dens4=pulse(t_cur_ind_left:t_cur_ind_right);
out_pulse_dens4=out_pulse(t_cur_ind_left:t_cur_ind_right);
%points_diff_dens4=points(2:end)-points(1:end-1);
% %Ne = loadD([outKey,num2str(plot_point(4)-1),'_Nsum_e_',location,'_T1.dat']);
% %Te = loadD([outKey,num2str(plot_point(4)-1),'_inst_temp_e_',location,'_T1.dat']);
% %Ne_dens4=Ne(t_cur_ind_left:t_cur_ind_right);
% %Te_dens4=Te(t_cur_ind_left:t_cur_ind_right);
% 
% outKey = '../../tMSBE-v3.6-optimization-analytic/run/out__';
% [t_cur,points,pulse,out_pulse]=loadPulse(outKey,plot_num(5),num_peaks,location);
% [w_qw_dens5, spectrum_dens5, spectrum_dens5_2] = getSpectrums_fromPoints(t_cur,width_spectrum(1),width_spectrum(2),points,pulse);
% t_cur_ind_left=find(t_cur>(points(plot_point(5))-width(1)),1);
% t_cur_ind_right=find(t_cur>(points(plot_point(5))+width(2)),1);
% if isempty(t_cur_ind_left)
%    t_cur_ind_left=1; 
% end
% if isempty(t_cur_ind_right)
%    t_cur_ind_right=length(t_cur); 
% end
% t_cur_dens5=t_cur(t_cur_ind_left:t_cur_ind_right);
% t_cur_dens5=t_cur_dens5-t_cur_dens5(1)-width(1);
% pulse_dens5=pulse(t_cur_ind_left:t_cur_ind_right);
% out_pulse_dens5=out_pulse(t_cur_ind_left:t_cur_ind_right);
% %points_diff_dens4=points(2:end)-points(1:end-1);
% % %Ne = loadD([outKey,num2str(plot_point(5)-1),'_Nsum_e_',location,'_T1.dat']);
% % %Te = loadD([outKey,num2str(plot_point(5)-1),'_inst_temp_e_',location,'_T1.dat']);
% % %Ne_dens5=Ne(t_cur_ind_left:t_cur_ind_right);
% % %Te_dens5=Te(t_cur_ind_left:t_cur_ind_right);



round_trip_time = loadD([outKey,'round_trip_time.dat']);
%round_trip_time= (points(2)-points(1))*ps;

if IO_finalPlotLog==1
    tmp_fig=figure(666);
    set(tmp_fig,'Name','Final Output Plot (Log)');
    set(tmp_fig,'Position',POS);
    semilogy(t_cur_dens1, out_pulse_dens1*cm*cm/1e6,t_cur_dens2, out_pulse_dens2*cm*cm/1e6,...
        t_cur_dens3, out_pulse_dens3*cm*cm/1e6, t_cur_dens4, out_pulse_dens4*cm*cm/1e6,...
        t_cur_dens2, abs(out_pulse_dens1-out_pulse_dens2)*cm*cm/1e6,...
        t_cur_dens3, abs(out_pulse_dens1-out_pulse_dens3)*cm*cm/1e6,...
        t_cur_dens4, abs(out_pulse_dens1-out_pulse_dens4)*cm*cm/1e6,'-g');
    hold on
    % for j=1:num_peaks
    %     semilogy((-j*round_trip_time/(2*ps))*[1,1],YLIM_intensity,'k--',...
    %     j*(round_trip_time/(2*ps))*[1,1],YLIM_intensity,'k--');
    % end
     %  semilogy(0*[1,1],YLIM_intensity,'k--');
    hold off

    xlim([-width(1),width(2)]);
    %ylim(YLIM_intensity)
    xlabel('time [ps]');
    ylabel('I [MW/cm2]');
    grid on
     legend('No thresholding','Total field thresholding',...
       'Field free Thresholding',...
         'Analytic field update','Error Total',...
            'Error Field-Free','Error Analytic');
    if IO_save==1
        saveas(tmp_fig,[saveKey,'finalOutputPlotLog.png']); 
    end
end

if IO_finalPlot==1
    tmp_fig=figure(667);
    set(tmp_fig,'Name','Final Output Plot');
    set(tmp_fig,'Position',POS);
    plot(t_cur_dens1, out_pulse_dens1*cm*cm/1e6,t_cur_dens2, out_pulse_dens2*cm*cm/1e6,...
        t_cur_dens3, out_pulse_dens3*cm*cm/1e6,...
        t_cur_dens4, out_pulse_dens4*cm*cm/1e6,...
        t_cur_dens2, abs(out_pulse_dens1-out_pulse_dens2)*cm*cm/1e6,...
        t_cur_dens3, abs(out_pulse_dens1-out_pulse_dens3)*cm*cm/1e6,...
        t_cur_dens4, abs(out_pulse_dens1-out_pulse_dens4)*cm*cm/1e6,'-g');
    %xlim(width_short);
    %ylim([0,5])
    xlabel('time [ps]');
    ylabel('I [MW/cm2]');
    grid on
     legend('No thresholding','Total field thresholding',...
       'Field free Thresholding',...
         'Analytic field update','Error Total',...
            'Error Field-Free','Error Analytic');
    if IO_save==1   
        saveas(tmp_fig,[saveKey,'finalOutputPlot.png']); 
    end
end

if IO_spectrumsPlot==1
    tmp_fig=figure(668);
    set(tmp_fig,'Name','Final Output Spectrums');
    set(tmp_fig,'Position',POS);
    semilogy(hbar*w_qw_dens1/e,abs(spectrum_dens1),hbar*w_qw_dens2/e,abs(spectrum_dens2),...
        hbar*w_qw_dens3/e,abs(spectrum_dens3), hbar*w_qw_dens4/e,abs(spectrum_dens4),...
    hbar*w_qw_dens1/e,abs(spectrum_dens1-spectrum_dens2),hbar*w_qw_dens1/e,abs(spectrum_dens1-spectrum_dens3),...
    hbar*w_qw_dens1/e,abs(spectrum_dens1-spectrum_dens4),'-g');
    xlim(XLIM_SPECTRUM);
    xlabel('Energy [eV]');
    ylabel('Spectrum [a.u.]');
    grid on
     legend('No thresholding','Total field thresholding',...
       'Field free Thresholding',...
         'Analytic field update','Error Total',...
            'Error Field-Free','Error Analytic');
    if IO_save==1
        saveas(tmp_fig,[saveKey,'finalOutputPlotSpectrum.png']); 
    end
end

if IO_spectrumsPlot2==1
    %Normalizing to original spectrum
    spectrum_dens1_norm=spectrum_dens1_2./spectrum_dens1;
    spectrum_dens2_norm=spectrum_dens2_2./spectrum_dens2;
    spectrum_dens3_norm=spectrum_dens3_2./spectrum_dens3;
    spectrum_dens4_norm=spectrum_dens4_2./spectrum_dens4;
    %spectrum_dens5_norm=spectrum_dens5_2./spectrum_dens5;
    
    tmp_fig=figure(669);
    set(tmp_fig,'Name','Final Output Spectrums (Secondary Pulse)');
    set(tmp_fig,'Position',POS);
       semilogy(hbar*w_qw_dens1/e,abs(spectrum_dens1_2),...
           hbar*w_qw_dens2/e,abs(spectrum_dens2_2),...
        hbar*w_qw_dens3/e,abs(spectrum_dens3_2),...
        hbar*w_qw_dens4/e,abs(spectrum_dens4_2),...
    hbar*w_qw_dens1/e,abs(spectrum_dens1_norm-spectrum_dens2_norm),...
    hbar*w_qw_dens1/e,abs(spectrum_dens1_norm-spectrum_dens3_norm),...
    hbar*w_qw_dens1/e,abs(spectrum_dens1_norm-spectrum_dens4_norm),'-g');
    xlim(XLIM_SPECTRUM);
    xlabel('Energy [eV]');
    ylabel('Spectrum [a.u.]');
    grid on
     legend('No thresholding','Total field thresholding',...
       'Field free Thresholding',...
         'Analytic field update','Error Total',...
            'Error Field-Free','Error Analytic');
    if IO_save==1
        saveas(tmp_fig,[saveKey,'finalOutputPlotSpectrum2.png']); 
    end
end

if IO_point_diff==1
    point_diff_min=[min(points_diff_dens1);min(points_diff_dens2);...
        min(points_diff_dens3);min(points_diff_dens4);min(points_diff_dens5);...
        min(points_diff_dens6)];
    point_diff_max=[max(points_diff_dens1);max(points_diff_dens2);...
        max(points_diff_dens3);max(points_diff_dens4);max(points_diff_dens5);...
        max(points_diff_dens6)];

    tmp_fig=figure(662);
    set(tmp_fig,'Name','Point Difference Overlay');
    set(tmp_fig,'Position',POS);
    plot(1:(num_peaks-1),points_diff_dens1,1:(num_peaks-1),points_diff_dens2,...
        1:(num_peaks-1),points_diff_dens3,1:(num_peaks-1),points_diff_dens4,...
        1:(num_peaks-1),points_diff_dens5,1:(num_peaks-1),points_diff_dens6);
    xlabel('time [ns]');
    ylabel('\Delta t [ps]');
    grid on
    legend(['t_0=',num2str((plot_num(1)-1)*5),'ns'],...
     ['t_0=',num2str((plot_num(2)-1)*5),'ns'],...
     ['t_0=',num2str((plot_num(3)-1)*5),'ns'],...
     ['t_0=',num2str((plot_num(4)-1)*5),'ns'],...
     ['t_0=',num2str((plot_num(5)-1)*5),'ns'],...
     ['t_0=',num2str((plot_num(6)-1)*5),'ns']);
    saveas(tmp_fig,[saveKey,'pointDifference_overlay.png']);

    tmp_fig=figure(663);
    set(tmp_fig,'Name','Point Difference Min');
    set(tmp_fig,'Position',POS);
    plot((plot_num-1)*5,point_diff_min);
    xlabel('time [ns]');
    ylabel('\Delta t [ps]');
    grid on
    saveas(tmp_fig,[saveKey,'pointDifference_min.png']);

    tmp_fig=figure(664);
    set(tmp_fig,'Name','Point Difference Max');
    set(tmp_fig,'Position',POS);
    plot((plot_num-1)*5,point_diff_max);
    xlabel('time [ns]');
    ylabel('\Delta t [ps]');
    grid on
    saveas(tmp_fig,[saveKey,'pointDifference_max.png']);

    tmp_fig=figure(665);
    point_diff_shift=(point_diff_max-point_diff_min)/2;
    set(tmp_fig,'Name','Point Difference Shift');
    set(tmp_fig,'Position',POS);
    plot((plot_num-1)*5,point_diff_shift);
    xlabel('time [ns]');
    ylabel('\Delta t [ps]');
    grid on
    saveas(tmp_fig,[saveKey,'pointDifference_shift.png']);
end
 
if IO_dens==1
    density_min=min([Ne_dens1;Ne_dens2;Ne_dens3]);
    density_max=max([Ne_dens1;Ne_dens2;Ne_dens3]);
    dens_lim=[density_min;density_max]*10^(-16);
    tmp_fig=figure(131);  
    set(tmp_fig,'Name','Density plot');
    set(tmp_fig,'Position',POS);
    plot(t_cur_dens1, Ne_dens1*10^(-16),t_cur_dens2, Ne_dens2*10^(-16),...
     t_cur_dens3, Ne_dens3*10^(-16),t_cur_dens4, Ne_dens4*10^(-16)...
    ,t_cur_dens5, Ne_dens5*10^(-16), t_cur_dens6, Ne_dens6*10^(-16));
hold on
for j=1:num_peaks
    plot((-j*round_trip_time/(2*ps))*[1,1],dens_lim,'k--',...
    j*(round_trip_time/(2*ps))*[1,1],dens_lim,'k--');
end
plot(0*[1,1],dens_lim,'k--');
hold off

ylim(dens_lim);
xlim([-width(1),width(2)]);
xlabel('time [ps]');
   yticks([1.83,1.85,1.865]);
ylabel('Density [10^{16}m^{-2}]');
grid on
 legend(['t_0=',num2str((plot_num(1)-1)*5),'ns'],...
     ['t_0=',num2str((plot_num(2)-1)*5),'ns'],...
     ['t_0=',num2str((plot_num(3)-1)*5),'ns'],...
     ['t_0=',num2str((plot_num(4)-1)*5),'ns'],...
     ['t_0=',num2str((plot_num(5)-1)*5),'ns'],...
     ['t_0=',num2str((plot_num(6)-1)*5),'ns'],...
     'T_{rt}/2');
   % legend('L_{OC}=L_{SESAM}','3L_{OC}=L_{SESAM}',...
   %     '7L_{OC}=L_{SESAM}','T_{rt}/2');
  % legend('L_{OC}=L_{SESAM}','3L_{OC}=L_{SESAM}','L_{OC}=3L_{SESAM}',...
  %      '7L_{OC}=L_{SESAM}','L_{OC}=7L_{SESAM}','T_{rt}');
saveas(tmp_fig,[saveKey,'DensityPlot.png']);
end
