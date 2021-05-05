close all
clear all

fs = 1.0e-15;
ps = 1.0e-12;
um = 1.0e-6;
nm = 1.0e-9;

global hbar;
global e;
global c0;

hbar = 1.054589e-34;
e = 1.602189e-19;
c0   = 2.99792458E+08;
eps0 = 8.854187817620E-12;

m0   = 9.109389754e-31;
me   = 0.06085*m0;
mh   = 0.234*m0;
a0   = 1.062146e-08;
mr = me*mh/(me+mh);
Eb = (hbar^2 / (2*mr*a0*a0));

set(0,'defaulttextinterpreter','tex') %Default
%set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontsize', 40)
set(0,'defaultAxesLinewidth', 4)
set(0,'defaultTextFontsize', 40)
set(0,'defaultlineMarkerSize',40)
set(0,'defaultlinelinewidth', 6) %Thin lines
maps %Load maps file for color scheme
POS=[1,1,1200,800];
POS2=[1,1,500,900];
%POS_plot=[0.14,0.1100,0.53125,0.85];
%POS_plot=[0.1330,0.1100,0.5930,0.8150];
my_lineStyle={'-k','-b','-r', '-g', '-m','-c','-k+'};

outKey = '../run/refSpecQW__';
location1='CAVOC'; %Field location for uploading and saving
location2='CAVOC'; %Field location for uploading and saving
location3=location2;
location4=location2;
location5=location2;
location6=location2;
plot_num =4; %Output number +1
date='091120';
test1='theta=2';
test2='theta=4';
test3='theta=6';
%test4='VCAV29-theta20';
%test5='VCAV32-theta25';
%test6='VCAV30-theta30';

testSave='boundaryOrder_com';
test_folder='test';
saveKey_local='Fall2019-Spring2020/VCAV/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey = [saveKey_fold,'/',testSave,'-out',num2str(plot_num-1),'-',location1,'-',location2,'-'];
user_entry = input(['saveKey=',saveKey,'? (y to continue): '], 's');
if user_entry~='y'
    gafuggle
end
if exist(saveKey_fold,'dir')~=7
    mkdir(saveKey_fold);
end


load([saveKey_fold,'/',test1,'-out3-',location1,'-gdd_pulsed_passive.mat'])
load([saveKey_fold,'/',test1,'-out3-',location1,'-gdd_stuff.mat'])
ml_E_qw=E_cav_qw;
ml_gdd_qw=gdd_qw;
ml_gdd_x=gdd_x;
ml_gdd_abs=gdd_abs;
ml_refSpecQW=refSpecQW;
ml_refSpecAbs=refSpecABS;
ml_w=w_gain;

load([saveKey_fold,'/',test2,'-out3-',location2,'-gdd_pulsed_passive.mat'])
load([saveKey_fold,'/',test2,'-out3-',location2,'-gdd_stuff.mat'])
ll_E_qw=E_cav_qw;
ll_gdd_qw=gdd_qw;
ll_gdd_x=gdd_x;
ll_gdd_abs=gdd_abs;
ll_refSpecQW=refSpecQW;
ll_refSpecAbs=refSpecABS;
ll_w=w_gain;

load([saveKey_fold,'/',test3,'-out3-',location3,'-gdd_pulsed_passive.mat'])
load([saveKey_fold,'/',test3,'-out3-',location3,'-gdd_stuff.mat'])
t3_E_qw=E_cav_qw;
t3_gdd_qw=gdd_qw;
t3_gdd_x=gdd_x;
t3_gdd_abs=gdd_abs;
t3_refSpecQW=refSpecQW;
t3_refSpecAbs=refSpecABS;
t3_w=w_gain;
% 
% load([saveKey_fold,'/',test4,'-out3-',location4,'-gdd_pulsed_passive.mat'])
% load([saveKey_fold,'/',test4,'-out3-',location4,'-gdd_stuff.mat'])
% t4_E_qw=E_cav_qw;
% t4_gdd_qw=gdd_qw;
% t4_gdd_x=gdd_x;
% t4_gdd_abs=gdd_abs;
% t4_refSpecQW=refSpecQW;
% t4_refSpecAbs=refSpecABS;
% t4_w=w_gain;
% 
% 
% load([saveKey_fold,'/',test5,'-out3-',location5,'-gdd_pulsed_passive.mat'])
% load([saveKey_fold,'/',test5,'-out3-',location5,'-gdd_stuff.mat'])
% t5_E_qw=E_cav_qw;
% t5_gdd_qw=gdd_qw;
% t5_gdd_x=gdd_x;
% t5_gdd_abs=gdd_abs;
% t5_refSpecQW=refSpecQW;
% t5_refSpecAbs=refSpecABS;
% t5_w=w_gain;
% 
% 
% load([saveKey_fold,'/',test6,'-out3-',location6,'-gdd_pulsed_passive.mat'])
% load([saveKey_fold,'/',test6,'-out3-',location6,'-gdd_stuff.mat'])
% t6_E_qw=E_cav_qw;
% t6_gdd_qw=gdd_qw;
% t6_gdd_x=gdd_x;
% t6_gdd_abs=gdd_abs;
% t6_refSpecQW=refSpecQW;
% t6_refSpecAbs=refSpecABS;
% t6_w=w_gain;

% asd %For micro polarization comparisons
% tmp_fig=figure;
% set(tmp_fig,'Position',POS);
% semilogy(t(1:10:end)/ps,abs(real(p_qw(1:10:end,14))),'r-',...
%      t(1:10:end)/ps,abs(real(p_fm1_qw(1:10:end,14))),'k-',...
%     t(1:10:end)/ps,abs(real(p_qw(1:10:end,14))-real(p_fm1_qw(1:10:end,14))),'b-' );
% grid on;
% xlabel('t [ps]');
% ylabel('real(E) [V/m]');
% lgd=legend('Mock-Linear', 'Linear', 'Difference');
% saveas(tmp_fig,[saveKey,'Preal_log.png']);



% tmp_fig=figure;
% set(tmp_fig,'Position',POS);
% plot((1:length(ml_E_qw))*0.1*fs/ps,real(ml_E_qw),'r-',...
%      (1:length(ml_E_qw))*0.1*fs/ps,real(ll_E_qw),'k-',...
%     (1:length(ml_E_qw))*0.1*fs/ps,abs(real(ml_E_qw)-real(ll_E_qw)),'b-' );
% grid on;
% xlabel('t [ps]');
% ylabel('real(E) [V/m]');
% lgd=legend(test1, test2, 'Difference');
% saveas(tmp_fig,[saveKey,'Ereal_log.png']);
% 
% tmp_fig=figure;
% set(tmp_fig,'Position',POS);
% plot((1:length(ml_E_qw))*0.1*fs/ps,imag(ml_E_qw),'r-',...
%      (1:length(ml_E_qw))*0.1*fs/ps,imag(ll_E_qw),'k-',...
%     (1:length(ml_E_qw))*0.1*fs/ps,abs(imag(ml_E_qw)-imag(ll_E_qw)),'b-' );
% grid on;
% xlabel('t [ps]');
% ylabel('imag(E) [V/m]');
% lgd=legend(test1, test2, 'Difference');
% saveas(tmp_fig,[saveKey,'Eimag_log.png']);
% 
% tmp_fig=figure;
% set(tmp_fig,'Position',POS);
% semilogy((1:length(ml_E_qw))*0.1*fs/ps,abs(ml_E_qw),'r-',...
%      (1:length(ml_E_qw))*0.1*fs/ps,abs(ll_E_qw),'k-',...
%     (1:length(ml_E_qw))*0.1*fs/ps,abs(abs(ml_E_qw-ll_E_qw)),'b-' );
% grid on;
% xlabel('t [ps]');
% ylabel('|E| [V/m]');
% lgd=legend(test1, test2, 'Difference');
% saveas(tmp_fig,[saveKey,'E_log.png']);

delta12=ml_E_qw-ll_E_qw;
delta13=ml_E_qw-t3_E_qw;
tmp_fig=figure;
set(tmp_fig,'Position',POS);
semilogy((1:length(ml_E_qw))*0.1*fs/ps,abs(ml_E_qw),(1:length(ll_E_qw))*0.1*fs/ps,...
    abs(ll_E_qw),(1:length(t3_E_qw))*0.1*fs/ps,abs(t3_E_qw),...
    (1:length(delta12))*0.1*fs/ps,abs(delta12),(1:length(delta13))*0.1*fs/ps,abs(delta13));
grid on
xlabel('time [ps]')
ylabel('|E| (V/m)')
legend('Front-Inner-Back','Front-Back-Inner', 'Back-Front-Inner', '\Delta_{1,2}','\Delta_{1,3}');
saveas(tmp_fig,[saveKey,'field_comparison.png']);

tmp_fig=figure;
set(tmp_fig,'Position',POS);
plot(hbar*ml_w/e ,(100*ml_refSpecQW),my_lineStyle{1},...
     hbar*ll_w/e,(100*ll_refSpecQW),my_lineStyle{2},...
     hbar*t3_w/e,(100*t3_refSpecQW),my_lineStyle{3});%,...
     %hbar*t4_w/e,(100*t4_refSpecQW),my_lineStyle{4},...
     %hbar*t5_w/e,(100*t5_refSpecQW),my_lineStyle{5},...
     %hbar*t6_w/e,(100*t6_refSpecQW),my_lineStyle{6},...
     %hbar*ll_w/e,(100*ll_refSpecAbs));
     %hbar*t5_w/e,(100*t5_refSpecQW),...
%hbar*ml_w/e, abs(100*(ml_refSpecQW-ll_refSpecQW)),'b-' )
grid on
xlim([1.14,1.26])
ylim([0,4])
xlabel('Energy [eV]')
ylabel('|R(w)| [%]')
legend('Front-Inner-Back','Front-Back-Inner', 'Back-Front-Inner');
%legend('n=1.8e16','n=1.85e16','n=1.9e16','n=2.0e16','n=2.1e16','n=2.2e16','n=5.0e14');
%legend('theta=0','theta=10','theta=15','theta=20','theta=25','theta=30');
saveas(tmp_fig,[saveKey,'gain_comparison.png']);

% 
% tmp_fig=figure;
% set(tmp_fig,'Position',POS);
% plot(hbar*ml_w/e ,(100*ml_refSpecQW),'r-',...
%      hbar*ll_w/e,(100*ll_refSpecQW),'k-',...
% hbar*ml_w/e, (100*(ml_refSpecQW-ll_refSpecQW)),'b-' )
% grid on
% xlim([1.16,1.26])
% ylim([-1,3])
% xlabel('Energy [eV]')
% ylabel('|R(w)| [%]')
% lgd=legend(test1, test2, 'Difference');
% saveas(tmp_fig,[saveKey,'gain.png']);
% 
% tmp_fig=figure;
%   set(tmp_fig,'Position',POS);
% semilogy(hbar*ml_gdd_x/e ,abs(ml_gdd_qw)/(fs*fs) ,'r-',...
%      hbar*ll_gdd_x/e ,abs(ll_gdd_qw)/(fs*fs) ,'k-',...
%      hbar*ml_gdd_x/e ,abs(ml_gdd_qw-ll_gdd_qw)/(fs*fs) ,'b-')
% xlabel('Energy [eV]')
% ylabel('GDD [(fs^2)]')
% grid on
% xlim([1.16,1.26])
% lgd=legend(test1, test2, 'Difference');
% saveas(tmp_fig,[saveKey,'GDD_log.png']);

tmp_fig=figure;
  set(tmp_fig,'Position',POS);
plot(hbar*ml_gdd_x/e ,(ml_gdd_qw)/(fs*fs) ,my_lineStyle{1},...
     hbar*ll_gdd_x/e ,(ll_gdd_qw)/(fs*fs) ,my_lineStyle{2},...
     hbar*t3_gdd_x/e ,(t3_gdd_qw)/(fs*fs) ,my_lineStyle{3});%,...
     %hbar*t4_gdd_x/e ,(t4_gdd_qw)/(fs*fs) ,my_lineStyle{4},...
     %hbar*t5_gdd_x/e ,(t5_gdd_qw)/(fs*fs) ,my_lineStyle{5},...    
     %hbar*t6_gdd_x/e ,(t6_gdd_qw)/(fs*fs) ,my_lineStyle{6});
xlabel('Energy [eV]')
ylabel('GDD [fs^2]')
grid on
xlim([1.16,1.26])
ylim([-2000,2000]);
%lgd=legend(test1, test2, 'Difference');
%legend('theta=0','theta=10','theta=15','theta=20','theta=25','theta=30');
legend('Front-Inner-Back','Front-Back-Inner', 'Back-Front-Inner');
saveas(tmp_fig,[saveKey,'GDD_comparison.png']);