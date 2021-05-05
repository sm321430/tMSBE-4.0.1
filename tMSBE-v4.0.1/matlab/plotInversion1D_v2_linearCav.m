%% plotInversion takes output from tMSBE file
%% 03/18/19: v1 created with basic implementaiton
%% 03/20/19: Added final output plotting
%% 06/24/19: v2: Added transparent surface plot for single figure output


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

set(0,'defaulttextinterpreter','tex') %Default
%set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontsize', 50)
set(0,'defaultAxesLinewidth', 5)
set(0,'defaultTextFontsize', 50)
set(0,'defaultlineMarkerSize',26)
set(0,'defaultlinelinewidth',5) %Thin lines
set(0,'defaultfigureposition',POS)

%IOs
IO_inversion=1; %Plots occupation numbers (all orders, all peaks)
IO_video=0; %Make and save video to desktop (needs IO_inversion to be on)
IO_inv_surf=0; %Surface plot with guidelines (needs IO_inversion to be on)
IO_finalOutput=1; %Plot final output
IO_pol=1; %Plot microscopic polarizations (all orders, all peaks)

max_k=5;
XLIM=[-250,250]; %Transverse plotting width (um)
KLIM=[0,max_k];
num_peaks=1;


outKey = '../run/refSpecQW__';
location='QW6'; %Field location for uploading and saving
plot_num =1; %Output number +1
date='060220';
test='linear';
test_folder='test';
saveKey_local='Fall2019-Spring2020/VCAV/';
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
w0 = loadD('../run/refSpecQW__w0.dat');
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])
round_trip_time = loadD([outKey,'round_trip_time.dat']);

t = loadD([outKey,num2str(plot_num-1),'_t.dat']); %Load time
t_0=t(1);
k_QW6=loadD([outKey,'K_QW6_T1.dat']); %k array: Future versions will load separate k-files 
KLIM(1)=min(k_QW6);

FOCUS = 5; % Compensate for lens in output: file output is from SESAM side, which has a spot (w/focus) where w is on the GAIN chip spot
Nt=length(t); %number of time points
t=t-t(1); %Shift time to start at zero
dt=(t(2)-t(1));
%t_init=42.6*ps;     %20.3*ps; %Start time
%t_max=44.0*ps;      %21.8*ps; %Max time
t_freq=4*fs; %output slice frequency
t_fast_freq=100*fs; %Fast slice frequency
SLICE_FREQ=round(t_freq/dt); %Frequency of frame slices
FAST_SLICE_FREQ=round(t_fast_freq/dt); %Frequency of fast frame slices
Nk=length(k_QW6);

%Cut off k values above a certain point with a check for appropriate bounds
k_ind_max=find(k_QW6>=max_k,1); %Plot only k values below max_k
if k_ind_max==0
   k_ind_max=length(k_QW6);
end

%% Load field data for peak finding
 
 pulse_re_qw=loadD([outKey,num2str(plot_num-1),'_E_re_QW6_T1.dat']);
 pulse_im_qw=loadD([outKey,num2str(plot_num-1),'_E_im_QW6_T1.dat']);
 pulse_qw=(pulse_re_qw+1i*pulse_im_qw).*exp(-1i*t*w0);
 out_pulse=0.5*eps0*c0*abs(pulse_qw);
 
 if IO_pol==1
 p_re_qw=loadDA([outKey,num2str(plot_num-1),'_p_re_QW6_T1.dat'],[Nt,Nk]);
 p_im_qw=loadDA([outKey,num2str(plot_num-1),'_p_im_QW6_T1.dat'],[Nt,Nk]);
 p_qw=p_re_qw+1i*p_im_qw;


% renorm_re=loadDA([outKey,num2str(plot_num-1),'_renormalized_pfp1_re_QW6_T1.dat'],[Nt,Nk]);
% renorm_im=loadDA([outKey,num2str(plot_num-1),'_renormalized_pfm1_re_QW6_T1.dat'],[Nt,Nk]);
% renorm_pfm3=renorm_pfm3_re+1i*renorm_pfm3_im;
 end
 
 if IO_inversion==1
    %% Loading variables: Requires knowledge of lengths [Nt,Nk]
    Ne_qw=loadDA([outKey,num2str(plot_num-1),'_ne_QW6_T1.dat'],[Nt,Nk]);
    Nh_qw=loadDA([outKey,num2str(plot_num-1),'_nh_QW6_T1.dat'],[Nt,Nk]);

    Inv_qw=Ne_qw+Nh_qw-1; 
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
plot(t/ps, out_pulse*cm*cm/1e6)  
xlabel('t (ps)');
ylabel('Intensity (MW/cm^2)');
saveas(tmp_fig,[saveKey,'total_field.png']);

tmp_fig=figure('Name', 'Phases of E');
plot(t/ps,real(pulse_qw),t/ps,imag(pulse_qw))
lgd=legend('real(E)','imag(E)');
xlabel('t (ps)');
ylabel('E (a.u.)');
saveas(tmp_fig,[saveKey,'ind_field.png']);
end

if IO_pol==1
for m=1:num_peaks
    t_init=points(m)-0.15;
    t_max=points(m)+0.25;
    t_mid_ind=find(t_cur>=points(m)-0.01,1);
    t_min_ind=find(t_cur>=t_init,1);
    t_max_ind=find(t_cur>=t_max,1);
    t_diff=t_max_ind-t_min_ind;
      
    tmp_fig=figure('Name','P_qw_k');
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),abs(p_qw(t_min_ind:t_max_ind,1:k_ind_max)));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'p_qw_',num2str(m),'.png']);
    
end   
    
end

if IO_inversion==1
%Inv_abs=Ne_abs-Nh_abs-1;

    for m=1:num_peaks
    t_init=points(m)-0.2;
    t_max=points(m)+0.2;
    t_mid_ind=find(t_cur>=points(m)-0.01,1);
    t_min_ind=find(t_cur>=t_init,1);
    t_max_ind=find(t_cur>=t_max,1);
    t_diff=t_max_ind-t_min_ind;
      
    tmp_fig=figure('Name','Ne+Nh-1')
    set(tmp_fig,'Position',POS);
    contourf(k_QW6(1:k_ind_max),t_cur(t_min_ind:t_max_ind),Inv_qw(t_min_ind:t_max_ind,1:k_ind_max));
    xlabel('k [1/a_0]');
    ylabel('t [ps]');
    colorbar
    saveas(tmp_fig,[saveKey,'Inversion_',num2str(m),'.png']);
    end
   

    %Set static axis meshes
    [X_line,K_line]=meshgrid(x_line,k_QW6(1:k_ind_max));
    [X,K]=meshgrid(transverse_grid_device_y(DEVICE_EXC+1:end-DEVICE_EXC)/um,k_QW6(1:k_ind_max));
    
    %Plot transparent surfaces
    for j=1:TOTAL_SLICES_SURF
        tmp=squeeze(Inv_qw(SURF_SLICE(j),:,:));
        surf(X,K,tmp,'FaceColor',surf_colormap(j,:), 'FaceAlpha',0.1, 'EdgeColor','none');
        set(gca,'Ydir','reverse');
        set(gca,'Xdir','reverse');
        if j==1
            hold on
        end
    end
    
    %Plot guidelines
    for j=1:TOTAL_LINES_SLICES_SURF
        tmp=squeeze(Inv_qw(LINES_SLICE(j),:,:));
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

if IO_video==1
    t_init=points(1)-250*fs;
    t_end=points(4)+390*fs;
    t_fast_lhs=250*fs;
    t_fast_rhs=400*fs;
    t_min_ind=find(t>=t_init,1);
    t_max_ind=find(t>=t_end,1);
    myVideo = VideoWriter([saveKey,'INVvideo.avi']);
    myVideo.FrameRate = 30;
  %  open(myVideo);
    fig=figure('Name', 'Inv video');
    slice=t_min_ind;
    j=1;
    
    while t(slice)<t_end
        tmp=squeeze((Inv_qw(slice,:,:)));
        plot(transpose(tmp),'edgecolor','none');
        ylim(XLIM);
        set(gca,'YTICK',[XLIM(1),0,XLIM(2)]);
        xlim(KLIM);
        set(gca,'FontSize',50)
        zlim(ZLIM);
        title(' ')
        %set(gca,'Xdir','reverse');
        %xlabel('k [1/a_0]');
        %ylabel(' x [\mu m] ');
        %zlabel('n_e+n_h-1');
        colormap(map3)
        caxis(color_axis);
        colorbar;
        
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
    else
        slice=slice+FAST_SLICE_FREQ;
        save_slice=[0;0;0];
    end
         title(['QW Inversion: t=',num2str(floor(t_0/ns)),'ns+',num2str(floor((t(slice)-t(1))/ps*100)/100),'ps  ']);
        
        drawnow
        
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
 %   writeVideo(myVideo, M);
 %   close(myVideo);
end

% 4.3 7.3
    %11.5 7.3
    %4.3 12.8
% 11.5 12.8

%1726,1241

