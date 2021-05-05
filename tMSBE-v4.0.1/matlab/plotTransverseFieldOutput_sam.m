%% Plot transver field outputs for tMSBE
%General output functions for tMSBE
clear all
close all

%IOs for each plot located here
IO_convergence=0;
IO_temp_dens=0;
IO_finalOutput=1;

%% Preliminaries
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
n=3.4453; %Background refractive index (change as needed)

set(0,'defaulttextinterpreter','tex') %Default
%set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontsize', 50)
set(0,'defaultAxesLinewidth', 5)
set(0,'defaultTextFontsize', 50)
set(0,'defaultlineMarkerSize',26)
set(0,'defaultlinelinewidth',5) %Thin lines
maps %Load maps file for color scheme
POS=[1,1,1200,800];
POS2=[1,1,500,900];
%POS_plot=[0.14,0.1100,0.53125,0.85];
POS_plot=[0.1330,0.1100,0.5930,0.8150];
my_lineStyle={'-kd','-bo','-r+', '-gx', '-c*','-ko','-k+','-kx','-k*','-b+','-bx','-b*','-g+'};

outKey = '../run/out__';
date='041620';
location='OUTPUT'; %Field location for uploading and saving
plot_num =1; %Output number +1
test='T6-VCAV';
saveKey_local='Summer2019/Modelocking/Graphics/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test];
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

round_trip_time = loadD([outKey,'round_trip_time.dat']);
transverse_grid_y = loadD([outKey,'transverse_grid_y.dat']);
transverse_grid_device_y = loadD([outKey,'transverse_grid_device_y.dat']);
dx=0.9*abs(transverse_grid_device_y(3)-transverse_grid_device_y(2));
[~, ind_device_y] = intersect(transverse_grid_y, transverse_grid_device_y);

FOCUS = 10; % Compensate for lens in output: file output is from SESAM side, which has a spot (w/focus) where w is on the GAIN chip spot

t = loadD([outKey,num2str(plot_num-1),'_t.dat']);

%Exclude edges in time
init_ind=find((t-t(1))/ps>=0,1);
final_ind=find((t-t(1))/ps>=(t(end)-t(1))/ps,1);
t=t(init_ind:final_ind);

%Loading and plotting parameters
num_peaks=1; %Number of peaks for computing peaks
num_peaks_plot=1;%num_peaks;
%tmp_width_pulse=5.0; 
PulseTLIM=[-0.2,0.2]; %Output plot temporal width (ps)
PulseTTICK=[-0.2,-0.1,0,0.1,0.2];
%tmp_width_pulse_log=1.5; %Output plot temporal width (ps) for log plot
spc_width=40; %Spacial width of plots (um)
spc_width_density=200;
PulseLogTLIM=[-5,15.0]; 
PulseLogTTICK=[-5,0,5,10,15.0]; %Log Pulse Plot temporal tick marks
PulseLogXTICK=[-spc_width,-20,0,20,spc_width]; %Log Pulse spatial tick marks
PowerTLIM=[-0.2,0.2]; %Power plot  temporal width
PowerTTICK=[-0.2,-0.1,0,0.1,0.2]; %Power plot tick marks
tmp_width_power=round_trip_time/(2*ps); %Energy integration parameter



NUM_TRANSVERSE = length(dir([outKey,num2str(plot_num-1),'_E_re_',location,'_T*.dat']));
NUM_TRANSVERSE_DEVICE = length(dir([outKey,num2str(plot_num-1),'_E_re_QW6_T*.dat']));
NUM_TRANSVERSE_DEVICES_SHORT=length(find(abs(transverse_grid_device_y)<spc_width*um+dx));
%pulse_tmp=zeros(NUM_TRANSVERSE_DEVICES_SHORT,final_ind-init_ind+1);
%out_pulse=zeros(NUM_TRANSVERSE_DEVICES_SHORT,final_ind-init_ind+1);
transverse_grid_short_y=zeros(NUM_TRANSVERSE_DEVICES_SHORT,1);
counter=0;

switch location
    case 'QW6'
    for i = 0:(NUM_TRANSVERSE_DEVICE-1)
        pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
        pulse_im = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
        pulse = (pulse_re(init_ind:final_ind) + 1i*pulse_im(init_ind:final_ind)).*exp(-1i*t*w0);
        if abs(transverse_grid_device_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_device_y(1+i);
            counter=counter+1;
        end
    end
    case 'ABS1'
    for i = 0:(NUM_TRANSVERSE_DEVICE-1)
        pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
        pulse_im = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
        pulse = (pulse_re(init_ind:final_ind) + 1i*pulse_im(init_ind:final_ind)).*exp(-1i*t*w0);
        if abs(transverse_grid_device_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_device_y(1+i);
            counter=counter+1;
        end
    end
    case 'OUTPUT'
    for i = 0:NUM_TRANSVERSE-1
        pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(i),'.dat']);
        pulse_im = loadD([outKey,num2str(plot_num-1),'_E_re_',location,'_T',num2str(i),'.dat']);
        pulse = (pulse_re(init_ind:final_ind) + 1i*pulse_im(init_ind:final_ind)).*exp(-1i*t*w0);
        if abs(transverse_grid_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
            pulse_tmp(1+counter,:)= pulse;
            out_pulse(1+counter,:) = 0.5*n*eps0*c0*abs(pulse).^2;
            transverse_grid_short_y(1+counter)=transverse_grid_y(1+i);
            counter=counter+1;
        end  
    end
    otherwise
    disp([location,' is not loadable. No data exists']);
    return;  
end
%% Find peak_num distinct peaks
t_cur=(t-t(1))/ps; %Current time refernced from t(1) in picoseconds
min_time=find(t_cur>=0,1); %First index for time>2ps
max_time=find(t_cur>=max(t_cur),1); %t_cur(end),1); %First index for max>2ps
ind_trans=ones(counter,20); %indices for peaks
peaks_trans=zeros(counter,20); %Values for peaks
points=zeros(num_peaks,1); %time values for peaksinitialized
tmp_pulse=out_pulse;%(:,min_time:max_time); %Temporary pulse variable for manipulation
% Find top peaks, away from edges for all transverse points
for ll=1:counter
    [peaks,time_ind]=findpeaks(tmp_pulse(ll,:),'SortStr','descend');
    ii=min(50,length(peaks));
    %ind_trans(ll,1:ii)=time_ind(1:ii);
    peaks_trans(ll,1:ii)=peaks(1:ii);
end

% Out of found peaks, find max separated peaks
jk=1;
while jk<=num_peaks
    [row,col]=find(tmp_pulse==max(max(peaks_trans)),1);
    if min(abs(points-t_cur(col)))>2.5
        points(jk)=t_cur(col);
        jk=jk+1;
    end
    tmp_pulse(row,col)=0;
    [row,col]=find(peaks_trans==max(max(peaks_trans)),1);
    peaks_trans(row,col)=0;
end
points=sort(points);
%% Plot final output
if IO_finalOutput==1
disp('Plotting final output figures')
tmp_fig=figure(11);
set(tmp_fig,'Name','Final Output');
set(tmp_fig,'Position',POS);
  contourf((t-t(1))/ps, transverse_grid_short_y/um,(out_pulse*cm*cm/1e6),'edgecolor','none')
  ylim([-spc_width,spc_width])
  colormap(map3);
  %set(gca,'Linewidth',4);
  xlim([min((t-t(1))/ps),max((t-t(1))/ps)]);
  colorbar
 saveas(tmp_fig,[saveKey,'finalOutput.png']);

%% Integrate power across transverse dimension
power=zeros(length(t),1);
energy=zeros(length(points),1);

for j=1:length(t)
power(j)=trapz(transverse_grid_short_y,out_pulse(:,j));
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
hold on
    if kk==1
        txt='t_0';
     elseif kk==2
         txt='t_0+T_{rt}';
    else
%        txt=['t_0+',num2str(floor(points(kk)-points(1))),'ps'];
        txt=['t_0+',num2str(kk-1),'T_{rt}'];
    end
plot(t_cur(1:10:length(t)),power(1:10:length(t))/max_power,my_lineStyle{kk},'DisplayName',txt);
pulse_ind_left=find(t_cur>-tmp_width_power,1);
if isempty(pulse_ind_left)
   pulse_ind_left=1;
end
pulse_ind_right=find(t_cur>tmp_width_power,1);
if isempty(pulse_ind_right)
   pulse_ind_right=length(t_cur); 
end
energy(kk)=trapz(t_cur(pulse_ind_left:pulse_ind_right),power(pulse_ind_left:pulse_ind_right));
 ylim([0,1]);
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
    hold off
     box on
     grid on
     legend show;
     saveas(tmp_fig,[saveKey,'QWpowerOverlay.png']);

% figure(2269)
% plot(0:length(points)-1,energy)
% xlabel('Round trip');
% ylabel('Pulse Energy(a.u.)');

%% Plots of individual output pulses 

for kk=1:length(points)
%A(kk)=subplot(1,length(points),kk);
tmp_fig=figure(1000+kk);
set(tmp_fig,'Name',['Final Output Pulse RT', num2str(kk-1)]);
set(tmp_fig,'Position',POS);
if (NUM_TRANSVERSE > 1)
    t_cur=(t-t(1))/ps-points(kk);
    if max(abs(transverse_grid_short_y))~=0
    contourf(t_cur, transverse_grid_short_y/um,(out_pulse*cm*cm/1e6),'edgecolor','none');
    ylim([-spc_width,spc_width])
    colormap(map3);
    caxis([min(min((out_pulse*cm*cm/1e6))) max(max((out_pulse*cm*cm/1e6)))])
    else
    surf((t-t(1))/ps, transverse_grid_y/um,out_pulse*cm*cm/1e6,'edgecolor','none')
    ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
    end
    xlim(PulseTLIM)
%     if kk==1
%         title('t_0','FontSize',24);
%     else
%         title(['t_0+',num2str(kk-1),'T_{RT}'],'FontSize',24);
%     end
%     if kk==1
%    ylabel('y [\mum]','Fontsize',24)
%     else
    set(gca,'YTick', [-spc_width,0,spc_width])
%     end
    zlabel('I [MW/cm^2]')
    set(gca,'XTick', PulseTTICK)
      grid on
else
    plot((t-t(1))/ps, out_pulse*cm*cm/1e6, 'b-')
    ylabel('I [MW/cm^2]');   
end
B=colorbar;
saveas(tmp_fig,[saveKey,'QWpulseRT',num2str(kk-1),'.png']);
end


tmp_fig=figure(100);
set(tmp_fig,'Name','Final Output Pulse Overlay');
hold on
set(tmp_fig,'Position',POS);
for kk=1:num_peaks_plot
    t_cur=(t-t(1))/ps-points(kk);
    pulse_ind=find(t_cur==0,1);
    if kk==1
        txt='t_0';
    elseif kk==2
        txt='t_0+T_{rt}';
    else
 %       txt=['t_0+',num2str(floor(points(kk)-points(1))),'ps'];
        txt=['t_0+',num2str(kk-1),'T_{rt}'];
    end
plot(transverse_grid_short_y/um,out_pulse(:,pulse_ind)*cm*cm/1e6,my_lineStyle{kk},'DisplayName',txt); 
end
xlim([-spc_width,spc_width]);
ylim([0,ceil(max(max(out_pulse/10))*cm*cm/1e6)*10]);
set(gca,'XTick', PulseLogXTICK);
set(gca,'YTick', [0,ceil(max(max(out_pulse/10))*cm*cm/1e6)*5,ceil(max(max(out_pulse/10))*cm*cm/1e6)*10]);
hold off
grid on
box on
legend show
saveas(tmp_fig,[saveKey,'QWpulseOverlay.png']);

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
end

%% Peak intensity and FWHM (time and transverse dim) and Pulse energy
if IO_convergence==1
disp('Plotting convergence figures')
max_I = zeros(counter, plot_num);
max_It = zeros(1,plot_num);
fwhm_t = zeros(1,plot_num);
fwhm_y = zeros(1,plot_num);
energy_evol=zeros(1,plot_num);
temperature_profile = zeros(NUM_TRANSVERSE_DEVICES_SHORT, 1);
density_profile = zeros(NUM_TRANSVERSE_DEVICES_SHORT, 1);
transverse_grid_short_y=zeros(NUM_TRANSVERSE_DEVICES_SHORT,1);
power=zeros(length(t),1);
out_pulse=0;

for j = 0:(plot_num-1)
    current_plot=j;
    t = loadD([outKey,num2str(j),'_t.dat']);
    trans_profile = zeros(NUM_TRANSVERSE,length(t));
    counter=0;
    for i = 0:(NUM_TRANSVERSE_DEVICE-1)
     if abs(transverse_grid_device_y(1+i))<spc_width*um+dx %max(abs(transverse_grid_y))/8
        pulse_re = loadD([outKey,num2str(j),'_E_re_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
        pulse_im = loadD([outKey,num2str(j),'_E_im_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
        tmp = loadD([outKey,'lattice_setup_','QW6_T',num2str(ind_device_y(1+i)),'.dat']);
        density_profile(1+counter) = tmp(1);
        temperature_profile(1+counter) = tmp(2);
        pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*w0);
        pulse_tmp(1+counter,1:length(pulse))= pulse;
        out_pulse(1+counter,1:length(pulse)) = 0.5*n*eps0*c0*abs(pulse).^2;
        out_pulse(1+counter,length(pulse)+1:end)=0;
       transverse_grid_short_y(1+counter)=transverse_grid_device_y(1+i);
        max_I(1+counter,1+j) = max(0.5*n*eps0*c0*abs(pulse).^2);
        counter=counter+1;
     end
    end
    for jj=1:length(t)
        power(jj)=trapz(transverse_grid_short_y,out_pulse(:,jj));
    end

    % Find a pulse to focus on
    if (NUM_TRANSVERSE>1)
        [~,ind] = findpeaks(out_pulse(ceil(counter/2),:),'sortstr','descend');      
        ind = sort(ind(1:2));
        fwhm_y(1+j) = findFWHM(transverse_grid_device_y, out_pulse(:,ind(1)) );
        %fwhm_y(1+j) = findSpotSize(transverse_grid_y,abs(trans_profile(:,ind(1))));
        fwhm_t(1+j) = findFWHM(t, out_pulse(ceil(counter/2),:) );
        pulse_ind_left=find((t-t(ind(1)))/ps>-tmp_width_power,1);
    if isempty(pulse_ind_left)
        pulse_ind_left=1; 
    end
    pulse_ind_right=find((t-t(ind(1)))/ps>tmp_width_power,1);
    if isempty(pulse_ind_right)
        pulse_ind_right=length(t);
    end
    energy_evol(j+1)=trapz(t(pulse_ind_left:pulse_ind_right),power(pulse_ind_left:pulse_ind_right));
   
    else
        fwhm_t(1+j) = findFWHM(t, abs(trans_profile) );
        energy_evol(j)=0;
    end
    
    max_It(1+j) = mean(t);
end

yy = max_I*cm*cm/1e6;
%yy = log(yy/max(yy(:)));

tmp_fig=figure(101);
set(tmp_fig,'Name','Convergence');
set(tmp_fig,'Position',POS);
subplot(4,1,1)
if (NUM_TRANSVERSE>1)
    surf(max_It/ns,transverse_grid_short_y/um,yy,'edgecolor','none')
    ylim([min(transverse_grid_short_y),max(transverse_grid_short_y)]/um)
    ylabel('y [\mum]','Fontsize',24)
    zlabel('I(t) [MW/cm^2]','Fontsize',24)
    view(25,22)
else
    plot(max_It/ns, yy)
    ylabel('I(t) [MW/cm^2]','Fontsize',24)
end
set(gca,'FontSize', 24)

zl = zlim;
hold on
yy2 = (density_profile-min(density_profile))/(max(density_profile)-min(density_profile));
plot3(0*ones(size(transverse_grid_short_y))/ps, transverse_grid_short_y/um, zl(1) + yy2*(zl(2)-zl(1)),'ko-')
hold off

xlabel('t [ns]')
grid on

%Plots central intensity
subplot(4,1,2)
if (NUM_TRANSVERSE>1)
    plot(max_It/ns,yy(ceil(counter/2),:),'bo-');
else
    plot(max_It/ns, yy);
end
grid on
xlabel('t [ns]','Fontsize',24)
ylabel('I(t) [MW/cm^2]','Fontsize',24)
set(gca,'FontSize', 24)

%Plots FWHM in time and space
subplot(4,1,3)
[AX,H1,H2] = plotyy(max_It/ns,fwhm_t/fs, max_It/ns,fwhm_y/um);
set(AX,{'ycolor'},{'b';'r'})
set(H1,'color','b','LineStyle','-','marker','o')
set(H2,'color','r','LineStyle','-','marker','o')
grid on
xlabel('t [ns]','Fontsize',24)

axes(AX(1))
ylabel('FWHM [fs]')
set(AX(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
set(gca,'FontSize', 24)

axes(AX(2))
ylabel('FWHM [\mum]')
set(AX(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
set(gca,'FontSize', 24) 

subplot(4,1,4)
plot(max_It/ns,energy_evol,'bo-')
grid on
xlabel('t [ns]')
ylabel('Energy [a.u.]')
set(gca,'FontSize', 24)
saveas(tmp_fig,[saveKey,'QWconvergence.png']);
end

% 
% figure
% [AX,H1,H2] = plotyy(transverse_grid_y/um, density_profile/1e16,transverse_grid_y/um, temperature_profile);
% set(AX,{'ycolor'},{'b';'r'})
% set(H1,'color','b','LineStyle','-','marker','o')
% set(H2,'color','r','LineStyle','-','marker','o')
% 
% grid on
% xlabel('y [\mum]')
% 
% axes(AX(1))
% ylim([0, 2.5])
% ylabel('Density [10^{16} m^{-2}]')
% set(AX(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
% 
% axes(AX(2))
% ylim([300, 390])
% ylabel('Temperature[K]')
% set(AX(2), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')

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
      Ne = loadD([outKey,num2str(plot_num-1),'_Nsum_e_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
      %Te = loadD([outKey,num2str(plot_num-1),'_inst_temp_e_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);

      qw_density(counter,:) = Ne(init_ind:final_ind);
      %qw_temp(1+i,:) = Te;
     
 %     Ne = loadD([outKey,num2str(plot_num-1),'_Nsum_e_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);
      %Te = loadD([outKey,num2str(plot_num-1),'_inst_temp_e_',location,'_T',num2str(ind_device_y(1+i)),'.dat']);

 %     abs_density(1+counter,:) = Ne;
      %abs_temp(1+i,:) = Te;      
      counter=counter+1;
    end
end
 

tmp_fig=figure(111);
set(tmp_fig,'Name','Density');
set(tmp_fig,'Position',POS2);
x_ind_left=find(transverse_grid_short_y/um>=-1.05*spc_width_density,1);
x_ind_right=find(transverse_grid_short_y/um>=spc_width_density,1);
contourf((t-t(1)-points(1)*ps)/round_trip_time,transverse_grid_short_y(x_ind_left:x_ind_right)/um,qw_density(x_ind_left:x_ind_right,:),'edgecolor','none')
caxis([min(min(qw_density(x_ind_left:x_ind_right,:))), max(max(qw_density(x_ind_left:x_ind_right,:)))]);
%ylim([-spc_width,spc_width]);
ylim([-spc_width_density,spc_width_density]);
xlim([-0.5*ps/round_trip_time,num_peaks_plot+0.01])%1.01*(points(num_peaks_plot)-points(1))/(round_trip_time/ps)])%num_peaks_plot+1*ps/round_trip_time])
set(gca,'YTick', [-spc_width_density,0,spc_width_density]);
ratio=(points(end)-points(1))/(round_trip_time/ps);
set(gca,'XTick',0:num_peaks-1);%ceil(ratio*10)/20,ceil(ratio*10)/10]);
set(gca,'FontSize',28);
set(gca,'Position',POS_plot);
%xlim([points(1)-2,max((t-t(1))/ps)])
%xlim([0,max((t-t(1))/ps)])
colormap('jet')
%ylabel('y [\mum]','FontSize',24);
%xlabel('t [ps]','Fontsize',24);
colorbar
saveas(tmp_fig,[saveKey,'QWdensity.png']);

% 
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
%  
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
% 
% 
% figure
% subplot(2,1,1)
% if (NUM_TRANSVERSE>1)
%     surf((t-t(1))/ps,transverse_grid_device_y/um,qw_density/1e16,'edgecolor','none')
%     ylabel('y [um]')
% else
% 	plot((t-t(1))/ps, qw_density/1e16)
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
% 	plot((t-t(1))/ps, qw_temp)
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
% 	plot((t-t(1))/ps, abs_density/1e14)
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
% 	plot((t-t(1))/ps, abs_temp)
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


