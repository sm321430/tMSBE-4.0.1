%% plotInversion takes output from tMSBE file
%% 03/18/19: v1 created with basic implementaiton
%% 03/20/19: Added final output plotting
%% 06/24/19: v2: Added transparent surface plot for single figure output
%% 10/10/20: Enables ABS1 inversion as well. Not setup for traditional devices
%% 02/17/21: Compares Kinetic holes from two separate peaks

clear all
close all

global um;
global ps;
global fs;
maps

%Physical dimensions and constants
fs = 1.0e-15;
ps = 1.0e-12;
um = 1.0e-6;
ns = 1.0e-9;
cm = 1.0e-2;
nm = 1.0e-9;

POS=[1,1,1200,800];
POS2=[1,1,750,1200];
maps %Load maps file for color scheme

setupPlot
setupConstants

%IOs
IO_inversion=1; %Loads occupation numbers (all orders, all peaks)
IO_inversion_plot=1; %Plots occupation numbers around found peaks
IO_video=1; %Make and save video to desktop (needs IO_inversion to be on)
IO_inv_surf=0; %(UNTESTED) Surface plot with guidelines (needs IO_inversion to be on) 
IO_finalOutput=1; %Plot final output
IO_pol=0; %Plot microscopic polarizations (all orders, all peaks)
IO_inversionKinHoleComp=0; %Compare Kinetic holes for ring cavity pulses

max_k=5;
XLIM=[-250,250]; %Transverse plotting width (um)
KLIM=[0,max_k];
num_peaks=1;


outKey ='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV84-1D-n2p9-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus10/run/out__';
date='022321';
test='tMSBE-RCAV84-1D-n2p9-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus10';
%outKey = '../run/out__';
location='ABS1'; %Field location for uploading and saving

%Automatic plot discovery
 plot_num=0;
 while isfile([char(outKey),num2str(plot_num),'_E_re_',char('OUTPUT'),'_T0.dat'])
       plot_num=plot_num+1; %Check for next plot. If exists, add to counter
 end
%plot_num=13; %Output number +1

test_folder='test';
saveKey_local='Fall2020-Summer2021/RingCAV/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey = [saveKey_fold,'/',test,'-out',num2str(plot_num-1),'-',location,'-'];
user_entry = input(['saveKey=',saveKey,'? (y to continue): '], 's');
if user_entry~='y'
    gafuggle
end
if exist(saveKey_fold,'dir')~=7
    mkdir(saveKey_fold);
end

%Physical constants
hbar = 1.054589e-34;
e = 1.602189e-19;
c0   = 2.99792458E+08;
mu0  = (4.0e-7)*pi;
eps0 = 1.0/(mu0*c0*c0);
a0   = 1.062146e-08;


%Output key
w0 = loadD([outKey,'w0.dat']);
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])
round_trip_time = loadD([outKey,'round_trip_time.dat']);

t = loadD([outKey,num2str(plot_num-1),'_t.dat']); %Load time
Nt=length(t); %number of time points
%Exclude edges in time
init_ind=find((t-t(1))/ps>=22.0,1);
final_ind=find((t-t(1))/ps>=23.5,1);
t=t(init_ind:final_ind);

t_0=t(1);
k_QW6=loadD([outKey,'K_QW6_T1.dat']); %k array: Future versions will load separate k-files 
KLIM(1)=min(k_QW6);

FOCUS = 5; % Compensate for lens in output: file output is from SESAM side, which has a spot (w/focus) where w is on the GAIN chip spot
t=t-t(1); %Shift time to start at zero
dt=(t(2)-t(1));
%t_init=42.6*ps;     %20.3*ps; %Start time
%t_max=44.0*ps;      %21.8*ps; %Max time
t_freq=5.0*fs; %output slice frequency
t_fast_freq=10*fs; %Fast slice frequency
SLICE_FREQ=round(t_freq/dt); %Frequency of frame slices
FAST_SLICE_FREQ=round(t_fast_freq/dt); %Frequency of fast frame slices
Nk=length(k_QW6);

%Cut off k values above a certain point with a check for appropriate bounds
k_ind_max=find(k_QW6>=max_k,1); %Plot only k values below max_k
if k_ind_max==0
   k_ind_max=length(k_QW6);
end

%% Load field data for peak finding

%  p_re_qw=loadDA([outKey,num2str(plot_num-1),'_p_re_QW6_T1.dat'],[Nt,Nk]);
%  p_im_qw=loadDA([outKey,num2str(plot_num-1),'_p_im_QW6_T1.dat'],[Nt,Nk]);
%  p_qw=p_re_qw+1i*p_im_qw;
%  
%  pulse_re_qw=loadD([outKey,num2str(plot_num-1),'_E_re_QW6_T1.dat']);
%  pulse_im_qw=loadD([outKey,num2str(plot_num-1),'_E_im_QW6_T1.dat']);
%  pulse_qw=pulse_re_qw+1i*pulse_im_qw;
%  out_pulse=0.5*eps0*c0*abs(pulse_qw);
%  

pulse_fp_re = loadD([outKey,num2str(plot_num-1),'_E_fp_re_',location,'_T1.dat']);
pulse_fp_re=pulse_fp_re(init_ind:final_ind);
pulse_fp_im = loadD([outKey,num2str(plot_num-1),'_E_fp_im_',location,'_T1.dat']);
pulse_fp_im=pulse_fp_im(init_ind:final_ind);
pulse_fp = (pulse_fp_re + 1i*pulse_fp_im).*exp(-1i*t*w0);

pulse_fm_re = loadD([outKey,num2str(plot_num-1),'_E_fm_re_',location,'_T1.dat']);
pulse_fm_im = loadD([outKey,num2str(plot_num-1),'_E_fm_im_',location,'_T1.dat']);
pulse_fm_re=pulse_fm_re(init_ind:final_ind);
pulse_fm_im=pulse_fm_im(init_ind:final_ind);
pulse_fm = (pulse_fm_re + 1i*pulse_fm_im).*exp(-1i*t*w0);

pulse_bp_re = loadD([outKey,num2str(plot_num-1),'_E_bp_re_',location,'_T1.dat']);
pulse_bp_im = loadD([outKey,num2str(plot_num-1),'_E_bp_im_',location,'_T1.dat']);
pulse_bp_re=pulse_bp_re(init_ind:final_ind);
pulse_bp_im=pulse_bp_im(init_ind:final_ind);
pulse_bp = (pulse_bp_re + 1i*pulse_bp_im).*exp(-1i*t*w0);

pulse_bm_re = loadD([outKey,num2str(plot_num-1),'_E_bm_re_',location,'_T1.dat']);
pulse_bm_im = loadD([outKey,num2str(plot_num-1),'_E_bm_im_',location,'_T1.dat']);
pulse_bm_re=pulse_bm_re(init_ind:final_ind);
pulse_bm_im=pulse_bm_im(init_ind:final_ind);
pulse_bm = (pulse_bm_re + 1i*pulse_bm_im).*exp(-1i*t*w0);
pulse=pulse_fp+pulse_fm+pulse_bp+pulse_bm;
out_pulse = 0.5*eps0*c0*abs(pulse_fp+pulse_fm+pulse_bp+pulse_bm).^2;


if IO_pol==1
p_fp1_re_qw=loadDA([outKey,num2str(plot_num-1),'_p_fp1_re_',location,'_T1.dat'],[Nt,Nk]);
p_fm1_re_qw=loadDA([outKey,num2str(plot_num-1),'_p_fm1_re_',location,'_T1.dat'],[Nt,Nk]);
p_fp3_re_qw=loadDA([outKey,num2str(plot_num-1),'_p_fp3_re_',location,'_T1.dat'],[Nt,Nk]);
p_fm3_re_qw=loadDA([outKey,num2str(plot_num-1),'_p_fm3_re_',location,'_T1.dat'],[Nt,Nk]);
p_bp1_re_qw=loadDA([outKey,num2str(plot_num-1),'_p_bp1_re_',location,'_T1.dat'],[Nt,Nk]);
p_bm1_re_qw=loadDA([outKey,num2str(plot_num-1),'_p_bm1_re_',location,'_T1.dat'],[Nt,Nk]);
p_bp3_re_qw=loadDA([outKey,num2str(plot_num-1),'_p_bp3_re_',location,'_T1.dat'],[Nt,Nk]);
p_bm3_re_qw=loadDA([outKey,num2str(plot_num-1),'_p_bm3_re_',location,'_T1.dat'],[Nt,Nk]);
p_fp1_im_qw=loadDA([outKey,num2str(plot_num-1),'_p_fp1_im_',location,'_T1.dat'],[Nt,Nk]);
p_fm1_im_qw=loadDA([outKey,num2str(plot_num-1),'_p_fm1_im_',location,'_T1.dat'],[Nt,Nk]);
p_fp3_im_qw=loadDA([outKey,num2str(plot_num-1),'_p_fp3_im_',location,'_T1.dat'],[Nt,Nk]);
p_fm3_im_qw=loadDA([outKey,num2str(plot_num-1),'_p_fm3_im_',location,'_T1.dat'],[Nt,Nk]);
p_bp1_im_qw=loadDA([outKey,num2str(plot_num-1),'_p_bp1_im_',location,'_T1.dat'],[Nt,Nk]);
p_bm1_im_qw=loadDA([outKey,num2str(plot_num-1),'_p_bm1_im_',location,'_T1.dat'],[Nt,Nk]);
p_bp3_im_qw=loadDA([outKey,num2str(plot_num-1),'_p_bp3_im_',location,'_T1.dat'],[Nt,Nk]);
p_bm3_im_qw=loadDA([outKey,num2str(plot_num-1),'_p_bm3_im_',location,'_T1.dat'],[Nt,Nk]);
p_fp1_qw=p_fp1_re_qw+1i*p_fp1_im_qw;
p_fm1_qw=p_fm1_re_qw+1i*p_fm1_im_qw;
p_fp3_qw=p_fp3_re_qw+1i*p_fp3_im_qw;
p_fm3_qw=p_fm3_re_qw+1i*p_fm3_im_qw;
p_bp1_qw=p_bp1_re_qw+1i*p_bp1_im_qw;
p_bm1_qw=p_bm1_re_qw+1i*p_bm1_im_qw;
p_bp3_qw=p_bp3_re_qw+1i*p_bp3_im_qw;
p_bm3_qw=p_bm3_re_qw+1i*p_bm3_im_qw;

p_fp1_qw=p_fp1_qw(init_ind:final_ind,:);
p_fm1_qw=p_fm1_qw(init_ind:final_ind,:);
p_fp3_qw=p_fp3_qw(init_ind:final_ind,:);
p_fm3_qw=p_fm3_qw(init_ind:final_ind,:);
p_bp1_qw=p_bp1_qw(init_ind:final_ind,:);
p_bm1_qw=p_bm1_qw(init_ind:final_ind,:);
p_bp3_qw=p_bp3_qw(init_ind:final_ind,:);
p_bm3_qw=p_bm3_qw(init_ind:final_ind,:);

renorm_pfp1_re=loadDA([outKey,num2str(plot_num-1),'_renormalized_pfp1_re_',location,'_T1.dat'],[Nt,Nk]);
renorm_pfm1_re=loadDA([outKey,num2str(plot_num-1),'_renormalized_pfm1_re_',location,'_T1.dat'],[Nt,Nk]);
renorm_pfp3_re=loadDA([outKey,num2str(plot_num-1),'_renormalized_pfp3_re_',location,'_T1.dat'],[Nt,Nk]);
renorm_pfm3_re=loadDA([outKey,num2str(plot_num-1),'_renormalized_pfm3_re_',location,'_T1.dat'],[Nt,Nk]);
renorm_pfp1_im=loadDA([outKey,num2str(plot_num-1),'_renormalized_pfp1_im_',location,'_T1.dat'],[Nt,Nk]);
renorm_pfm1_im=loadDA([outKey,num2str(plot_num-1),'_renormalized_pfm1_im_',location,'_T1.dat'],[Nt,Nk]);
renorm_pfp3_im=loadDA([outKey,num2str(plot_num-1),'_renormalized_pfp3_im_',location,'_T1.dat'],[Nt,Nk]);
renorm_pfm3_im=loadDA([outKey,num2str(plot_num-1),'_renormalized_pfm3_im_',location,'_T1.dat'],[Nt,Nk]);
renorm_pfp1=renorm_pfp1_re+1i*renorm_pfp1_im;
renorm_pfm1=renorm_pfm1_re+1i*renorm_pfm1_im;
renorm_pfp3=renorm_pfp3_re+1i*renorm_pfp3_im;
renorm_pfm3=renorm_pfm3_re+1i*renorm_pfm3_im;  
end

if IO_inversion==1
   %% Loading variables: Requires knowledge of lengths [Nt,Nk]
Ne_00_qw=loadDA([outKey,num2str(plot_num-1),'_ne_00_',location,'_T1.dat'],[Nt,Nk]);
Ne_p2_qw_real=loadDA([outKey,num2str(plot_num-1),'_ne_p2_re_',location,'_T1.dat'],[Nt,Nk]);
%Ne_m2_qw_real=loadDA([outKey,num2str(plot_num-1),'_ne_m2_re_',location,'_T1.dat'],[Nt,Nk]);
%Ne_p2_qw_imag=loadDA([outKey,num2str(plot_num-1),'_ne_p2_im_',location,'_T1.dat'],[Nt,Nk]);
%Ne_m2_qw_imag=loadDA([outKey,num2str(plot_num-1),'_ne_m2_im_',location,'_T1.dat'],[Nt,Nk]);
Nh_00_qw=loadDA([outKey,num2str(plot_num-1),'_nh_00_',location,'_T1.dat'],[Nt,Nk]);
Nh_p2_qw_real=loadDA([outKey,num2str(plot_num-1),'_nh_p2_re_',location,'_T1.dat'],[Nt,Nk]);
%Nh_m2_qw_real=loadDA([outKey,num2str(plot_num-1),'_nh_m2_re_',location,'_T1.dat'],[Nt,Nk]);
%Nh_p2_qw_imag=loadDA([outKey,num2str(plot_num-1),'_nh_p2_im_',location,'_T1.dat'],[Nt,Nk]);
%Nh_m2_qw_imag=loadDA([outKey,num2str(plot_num-1),'_nh_m2_im_',location,'_T1.dat'],[Nt,Nk]);
Ne_p2_qw=2.0*(Ne_p2_qw_real);%+1i*Ne_p2_qw_imag;
Nh_p2_qw=2.0*(Nh_p2_qw_real);%+1i*Nh_p2_qw_imag;
%Ne_m2_qw=Ne_m2_qw_real+1i*Ne_m2_qw_imag;
%Nh_m2_qw=Nh_m2_qw_real+1i*Nh_m2_qw_imag;

Ne_00_qw=Ne_00_qw(init_ind:final_ind,:);
Nh_00_qw=Nh_00_qw(init_ind:final_ind,:);
Ne_p2_qw=Ne_p2_qw(init_ind:final_ind,:);
Nh_p2_qw=Nh_p2_qw(init_ind:final_ind,:);
%Ne_m2_qw=Ne_m2_qw(init_ind:final_ind,:);
%Nh_m2_qw=Nh_m2_qw(init_ind:final_ind,:);

if strcmp(location,'QW6')
    Inv_00_qw=Ne_00_qw+Nh_00_qw-1;
else
    Inv_00_qw=Ne_00_qw+Nh_00_qw;
end
Inv_p2_qw=Ne_p2_qw+Nh_p2_qw;
%Inv_m2_qw=Ne_m2_qw+Nh_m2_qw;
Inv_qw=Inv_00_qw+2.0*real(Inv_p2_qw);
end

%% Find peak_num distinct peaks
t_cur=(t-t(1))/ps; %Current time refernced from t(1) in picoseconds
min_time=1;%find(t_cur>3,1); %First index for time>2ps
max_time=length(t);%find(t_cur>t(end)-3,1); %First index for max>2ps
points=0;%-round_trip_time/ps*ones(num_peaks,1); %time values for peaksinitialized
tmp_pulse=out_pulse(min_time:max_time); %Temporary pulse variable for manipulation
t_cur=t_cur(min_time:max_time);

[peaks,time_ind]=findpeaks(tmp_pulse,'SortStr','descend');
ii=min(20,length(peaks));
peaks=peaks(1:ii);

% Out of found peaks, find max separated peaks
jk=1;
while jk<=num_peaks
    tmp_ind=find(tmp_pulse==max(max(peaks)),1);
    if (min(abs(points-t_cur(tmp_ind)))>1)||((t_cur(tmp_ind)<1) && max(points(points<1))==0)
        points(jk)=t_cur(tmp_ind);
        jk=jk+1;
    end
    tmp_pulse(tmp_ind)=0;
    peaks_ind=find(peaks==max(max(peaks)),1);
    peaks(peaks_ind)=0;
end
points=sort(points);

%% Plot final output
if IO_finalOutput==1
    tmp_fig=figure('Name', 'Field intensity');
    set(tmp_fig,'Position',POS);
    plot(t/ps, out_pulse*cm*cm/1e6);
    xlim([t(1),t(end)]/ps);
    xlabel('t (ps)');
    ylabel('Intensity (MW/cm^2)');
    saveas(tmp_fig,[saveKey,'total_field.png']);

    tmp_fig=figure('Name', 'Real(E)');
    set(tmp_fig,'Position',POS);
    plot(t/ps,real(pulse_fp),t/ps,real(pulse_fm),t/ps,real(pulse_bp),t/ps,real(pulse_bm))
    xlim([t(1),t(end)]/ps);
    lgd=legend('FP','FM','BP','BM');
    legend('show');
    set(lgd,'Location','northeast');
    legend
    xlabel('t (ps)');
    ylabel('real(E) (a.u.)');
    saveas(tmp_fig,[saveKey,'real_ind_field.png']);


    tmp_fig=figure('Name', 'Imag(E)');
    set(tmp_fig,'Position',POS);
    plot(t/ps,imag(pulse_fp),t/ps,imag(pulse_fm),t/ps,imag(pulse_bp),t/ps,imag(pulse_bm))
    xlim([t(1),t(end)]/ps);
    lgd=legend('FP','FM','BP','BM');
    legend('show');
    set(lgd,'Location','northeast');
    xlabel('t (ps)');
    ylabel('imag(E) (a.u.)');
    saveas(tmp_fig,[saveKey,'imag_ind_field.png']);
end

if IO_pol==1
for m=1:num_peaks
    t_init=points(m)-0.3;
    t_max=points(m)+0.5;
    t_mid_ind=find(t_cur>=points(m)-0.01,1);
    t_min_ind=find(t_cur>=t_init,1);
    t_max_ind=find(t_cur>=t_max,1);
    t_diff=t_max_ind-t_min_ind;
      
    tmp_fig=figure('Name','P_fp1_k');
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(p_fp1_qw(t_min_ind:t_max_ind,1:k_ind_max)));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'p_fp1_',num2str(m),'.png']);
    
        tmp_fig=figure('Name','P_fm1_k');
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(p_fm1_qw(t_min_ind:t_max_ind,1:k_ind_max)));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'p_fm1_',num2str(m),'.png']);
    
        tmp_fig=figure('Name','P_fp3_k');
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(p_fp3_qw(t_min_ind:t_max_ind,1:k_ind_max)));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'p_fp3_',num2str(m),'.png']);
    
        tmp_fig=figure('Name','P_fm3_k');
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(p_fm3_qw(t_min_ind:t_max_ind,1:k_ind_max)));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'p_fm3_',num2str(m),'.png']);
    
        tmp_fig=figure('Name','P_bp1_k');
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(p_bp1_qw(t_min_ind:t_max_ind,1:k_ind_max)));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'p_bp1_',num2str(m),'.png']);
    
        tmp_fig=figure('Name','P_bm1_k');
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(p_bm1_qw(t_min_ind:t_max_ind,1:k_ind_max)));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'p_bm1_',num2str(m),'.png']);
    
        tmp_fig=figure('Name','P_bp3_k');
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(p_bp3_qw(t_min_ind:t_max_ind,1:k_ind_max)));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'p_bp3_',num2str(m),'.png']);
    
        tmp_fig=figure('Name','P_bm3_k');
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(p_bm3_qw(t_min_ind:t_max_ind,1:k_ind_max)));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'p_bm3_',num2str(m),'.png']); 
%     
%         tmp_fig=figure('Name','renorm_P_fp1_k');
%     set(tmp_fig,'Position',POS);
%     contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(renorm_pfp1(t_min_ind:t_max_ind,1:k_ind_max)));
%     xlabel('k [1/a_0]');
%     ylabel('t [ps]');
%     colorbar
%     saveas(tmp_fig,[saveKey,'renorm_pfp1_',num2str(m),'.png']);
%     
%             tmp_fig=figure('Name','renorm_P_fm1_k');
%     set(tmp_fig,'Position',POS);
%     contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(renorm_pfm1(t_min_ind:t_max_ind,1:k_ind_max)));
%     xlabel('k [1/a_0]');
%     ylabel('t [ps]');
%     colorbar
%     saveas(tmp_fig,[saveKey,'renorm_pfm1_',num2str(m),'.png']);
%     
%             tmp_fig=figure('Name','renorm_P_fp3_k');
%     set(tmp_fig,'Position',POS);
%     contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(renorm_pfp3(t_min_ind:t_max_ind,1:k_ind_max)));
%     xlabel('k [1/a_0]');
%     ylabel('t [ps]');
%     colorbar
%     saveas(tmp_fig,[saveKey,'renorm_pfp3_',num2str(m),'.png']);
%     
%         tmp_fig=figure('Name','renorm_P_fm3_k');
%     set(tmp_fig,'Position',POS);
%     contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(renorm_pfm3(t_min_ind:t_max_ind,1:k_ind_max)));
%     xlabel('k [1/a_0]');
%     ylabel('t [ps]');
%     colorbar
%     saveas(tmp_fig,[saveKey,'renorm_pfm3_',num2str(m),'.png']);
    
end   
    
end

%% Plot kinetic hole comparisons for ring cavity pulses
if IO_inversionKinHoleComp==1
    points_pre=points-0.2;
    points_post=points+0.07;
    t_ind_pre=[find(t/ps>=points_pre(1),1),find(t/ps>=points_pre(2),1)];
    t_ind_post=[find(t/ps>=points_post(1),1),find(t/ps>=points_post(2),1)];
    Inv_slice_pre=Inv_qw(t_ind_pre,:);
    Inv_slice_post=Inv_qw(t_ind_post,:);
    Inv_slices=[Inv_slice_pre;Inv_slice_post];
    tmp_fig=figure;
    plot(k_QW6,Inv_slices(2,:),'-b',k_QW6,Inv_slices(1,:),'-r', ...
        k_QW6,Inv_slices(4,:),':b',k_QW6,Inv_slices(3,:),':r');
    grid on
    xlim(KLIM);
    ylim([-0.3,0.5]); 
    xlabel('k [1/a_0]');
    ylabel('Population Inversion'); 
    lgd=legend('Before CW','Before CCW', 'After CW', 'After CCW');
    set(lgd,'Location','Best');
    saveas(tmp_fig,[saveKey,'inversionSlices.png']);
end

if IO_inversion_plot==1
%Inv_abs=Ne_abs-Nh_abs-1;

    for m=1:num_peaks
    t_init=points(m)-0.5;
    t_max=points(m)+1.0;
    t_mid_ind=find(t_cur>=points(m)-0.01,1);
    t_min_ind=find(t_cur>=t_init,1);
    t_max_ind=find(t_cur>=t_max,1);
    t_diff=t_max_ind-t_min_ind;
      
    if strcmp(location,'QW6')
        tmp_fig=figure('Name','Ne_00+Nh_00-1');
    else
        tmp_fig=figure('Name','Ne_00+Nh_00');
    end
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),Inv_00_qw(t_min_ind:t_max_ind,1:k_ind_max));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'Inversion_00_',num2str(m),'.png']);
    
    tmp_fig=figure('Name','Ne_p2+Nh_p2')
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(Inv_p2_qw(t_min_ind:t_max_ind,1:k_ind_max)));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'Inversion_p2_',num2str(m),'.png']);
%     
%     tmp_fig=figure('Name','Ne_m2+Nh_m2')
%     set(tmp_fig,'Position',POS);
%     contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(Inv_m2_qw(t_min_ind:t_max_ind,1:k_ind_max)));
%     xlabel('k [1/a_0]');
%     ylabel('t [ps]');
%     colorbar
%     saveas(tmp_fig,[saveKey,'Inversion_m2_',num2str(m),'.png']);
    end

    if IO_inv_surf==1
        %Set static axis meshes
        [X_line,K_line]=meshgrid(x_line,k_QW6(1:k_ind_max));
        [X,K]=meshgrid(transverse_grid_device_y(DEVICE_EXC+1:end-DEVICE_EXC)/um,k_QW6(1:k_ind_max));
        
        %Plot transparent surfaces
        for j=1:TOTAL_SLICES_SURF
            tmp=squeeze(Inv_qw(SURF_SLICE(j),:));
            surf(X,K,tmp,'FaceColor',surf_colormap(j,:), 'FaceAlpha',0.1, 'EdgeColor','none');
            set(gca,'Ydir','reverse');
            set(gca,'Xdir','reverse');
            if j==1
                hold on
            end
        end
        
        %Plot guidelines
        for j=1:TOTAL_LINES_SLICES_SURF
            tmp=squeeze(Inv_qw(LINES_SLICE(j),:));
            for l=1:TOTAL_LINES_SURF
              line_mesh(:,l)=tmp(:,line_ind(l));
            end
            h=surf(transpose(X_line),transpose(K_line),transpose(line_mesh),'EdgeColor',line_colormap(j,:),'FaceColor','none','MeshStyle','row');
            %set(gca,'Ydir','reverse');
        end  
        xlim(XLIM);
        ylim(KLIM);
        zlim(ZLIM);
        xticks([XLIM(1),round(XLIM(1)/2),0,round(XLIM(2)/2),XLIM(2)]);
        dKLIM=round((KLIM(2)-KLIM(1))/4);
        yticks([KLIM(1),KLIM(1)+dKLIM,KLIM(1)+dKLIM*2,KLIM(1)+dKLIM*3,KLIM(2)]);
        dZLIM=round((10*ZLIM(2)-10*ZLIM(1))/4);
        %zticks([ZLIM(1),ZLIM(1)+dZLIM/10,ZLIM(1)+dZLIM/5,ZLIM(1)+dZLIM/10*3,ZLIM(2)]);
        zticks([-0.2,0,0.2,0.45]);
        grid on
        box on
        hold off
    %tmp_x=xlabel('   y(\mu m)','FontName','Arial','FontSize',30,'rotation',3,...
       %     'VerticalAlignment','mid','HorizontalAlignment','left');
       % tmp_z=zlabel('n_e+n_h-1','FontName','Arial','FontSize',30,...
       %     'rotation',90,'VerticalAlignment','baseline','HorizontalAlignment','Center');
      %  tmp_y=ylabel('    k(a^{-1}_0)','FontName','Arial','FontSize',28,'rotation',0,...
       %     'HorizontalAlignment','left','VerticalAlignment','mid');
        az=165.7000;
        el=30.8000;
        view(az,el);
      %  tmp=title([num2str(m-1),'RT'],'FontSize',20);
      %  tmp_pos=get(tmp,'Position');
       % set(tmp,'Position',[tmp_pos(1)+25,tmp_pos(2)-2.5,tmp_pos(3)])
       saveas(tmp_fig,[saveKey,'Inversion-',num2str(m-1),'RT.png']);
    end
end

if IO_video==1
    points=points*ps
    t_init=points(1)-380*fs;%t(1);%points(1)-250*fs;
    if length(points)>1
        t_init2=points(2)-380*fs;
    else
        t_init2=t_init;
    end
    t_end=points(1)+0.51*ps;%t(end)-600*fs;%points(4)+390*fs;
    t_fast_lhs=200*fs;
    t_fast_rhs=200*fs;
    t_min_ind=find(t>=t_init,1);
    t_min_ind2=find(t>=t_init2,1);
    t_max_ind=find(t>=t_end,1);
    myVideo = VideoWriter([saveKey,'INVvideo.avi']);
    myVideo.FrameRate = 10;
    open(myVideo);
    fig=figure('Name', 'Population video');
    slice=t_min_ind;
    slice2=t_min_ind2;
    save_slice=[0;0;0];
    j=1;
    peakInv=max(max(Inv_qw));
    while t(slice)<t_end
        tmp=Inv_qw(slice,:);
        %tmp2=Inv_qw(slice2,:);
        %plot(k_QW6,tmp,k_QW6,tmp2);
        plot(k_QW6,tmp);
        %ylim(XLIM);
        %set(gca,'YTICK',[XLIM(1),0,XLIM(2)]);
        xlim(KLIM);
        set(gca,'FontSize',50)
        grid on
        %zlim(ZLIM);
        title(' ')
        %set(gca,'Xdir','reverse');
        xlabel('k [1/a_0]');
        if strcmp(location,'QW6')
            ylabel('Population Inversion');
            ylim([-0.3,0.5]);   
        else
           ylabel('Population Inversion'); 
           ylim([-0.5,1.5]);
        end
        %zlabel('n_e+n_h-1');
        %colormap(map3)
        %caxis(color_axis);
        %colorbar;
        
        %Various views
        %az=60; %13.7;
        %el=35; %36.4;
        %view([az,el]);
    
        %Save specific slices
    ind=find(abs(points-t(slice))==min(abs(points-t(slice))),1);
    if t(slice)>=points(ind)-t_fast_lhs && t(slice)<=points(ind)+t_fast_rhs
        if t(slice)>=points(ind)-t_fast_lhs && save_slice(1)==0 %Before peak
        saveas(gcf,[saveKey,'InversionPrePeak-',num2str(ind-1),'RT.png']);
        save_slice(1)=1;
        end
        if t(slice)>=points(ind) && save_slice(2)==0 %Peak
        saveas(gcf,[saveKey,'InversionPeak-',num2str(ind-1),'RT.png']);
        save_slice(2)=1;
        end
        if t(slice)>=points(ind)+100*fs && save_slice(3)==0 %100fs after peak
        saveas(gcf,[saveKey,'InversionPostPeak-',num2str(ind-1),'RT.png']);
        save_slice(3)=1;
        end
        slice=slice+SLICE_FREQ;
        slice2=slice2+SLICE_FREQ;
    else
        slice=slice+FAST_SLICE_FREQ;
        slice2=slice2+FAST_SLICE_FREQ;
        save_slice=[0;0;0];
    end
        %title(['QW Inversion: t=',num2str(floor(t_0/ns)),'ns+',num2str(floor((t(slice)-t(1))/ps*100)/100),'ps  ']);
        txt=['t=',num2str(floor(t_0/ns)),'ns+',num2str(floor((t(slice)-t(1))/ps*100)/100),'ps'];
        txt2=['t=',num2str(floor(t_0/ns)),'ns+',num2str(floor((t(slice2)-t(1))/ps*100)/100),'ps'];
        lgd=legend(txt,txt2);
        set(lgd,'Location','northeast');
        %text(0.6*KLIM(2),0.8*peakInv,txt);
        %drawnow
        
    %Video creation/saving tools
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    F = getframe(gca);
    ax.Units = 'normalized';
    M(j)=getframe(gcf);
    
    %Update counter
    j=j+1;
    
    
 %  movie(fig,M);
    end
    writeVideo(myVideo, M);
    close(myVideo);
end

% 4.3 7.3
    %11.5 7.3
    %4.3 12.8
% 11.5 12.8

%1726,1241

