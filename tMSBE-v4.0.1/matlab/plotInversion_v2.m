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
IO_inversion=1; %Loads occupation numbers
IO_video=0; %Make and save video to desktop (needs IO_inversion to be on)
IO_inv_surf=0; %Surface plot with guidelines (needs IO_inversion to be on)
IO_finalOutput=0; %Plot final output

DEVICE_EXC=22; %Excludes side devices
max_k=4;
XLIM=[-250,250]; %Transverse plotting width (um)
KLIM=[0,max_k];
ZLIM=[-0.2,0.45]; %Height
POS=[267,138,873,657];
spc_width=max(abs(XLIM)); %(um)
plot_num = 77; %Output number +1
num_peaks=5;
%videoFile='../../tMSBE_June24.avi';

TOTAL_SLICES_SURF=2; %Number of surface slices
TOTAL_LINES_SURF=5; %Number of lines on a slice, divides domain equally)
TOTAL_LINES_SLICES_SURF=3; %Number of slices of lines
SLICE_FREQ_SURF=710; %Surface frequency plot (currently does nothing)
SLICE_FREQ_LINE=355; %Line frequency plot (currently does nothing)

outKey = '../run/refSpecQW__';
location='CAVP'; %Field location for uploading and saving
plot_num =4; %Output number +1
date='050520';
test='1DLinear-noCoulomb';
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

%Setting colormaps for surface plots
%Surfaces in transparent surface and guideline plot
surf_colormap=[0.0 0.0 0.0
               1.0 0.0 0.0
               0.0 0.0 0.0
               0.0 0.0 1.0]; 


%Lines in transparent surface and guideline plot
line_colormap=[0.0 0.0 0.0
               0.5 0.0 0.0
               0.8 0.0 0.0
               0.0 0.0 0.0
               0.0 0.0 0.5
               0.0 0.0 0.8]; 
%line_colormap=surf_colormap;

%Output key
outKey = '../run/out__';
w0 = loadD([outKey,'w0.dat']);
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])
round_trip_time = loadD([outKey,'round_trip_time.dat']);
transverse_grid_y = loadD([outKey,'transverse_grid_y.dat']);
transverse_grid_device_y = loadD([outKey,'transverse_grid_device_y.dat']);
[~, ind_device_y] = intersect(transverse_grid_y, transverse_grid_device_y);
NUM_TRANSVERSE = length(dir([outKey,num2str(plot_num-1),'_E_re_QW6_T*.dat']))
NUM_TRANSVERSE_DEVICE = length(dir([outKey,num2str(plot_num-1),'_E_re_QW6_T*.dat']))

transverse_grid_y = loadD([outKey,'transverse_grid_y.dat']);
transverse_grid_device_y = loadD([outKey,'transverse_grid_device_y.dat']);
[~, ind_device_y] = intersect(transverse_grid_y, transverse_grid_device_y);
t = loadD([outKey,num2str(plot_num-1),'_t.dat']); %Load time
t_0=t(1);
k_QW6=loadD([outKey,'K_QW6_T',num2str(ind_device_y(1)),'.dat']); %k array: Future versions will load separate k-files 
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
counter=0;
transverse_grid_short_y=0;
for i = 0:(NUM_TRANSVERSE_DEVICE-1)
    pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
    pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_QW6_T',num2str(ind_device_y(1+i)),'.dat']);
    pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*w0);
if abs(transverse_grid_device_y(1+i))<spc_width*um %max(abs(transverse_grid_y))/8
        pulse_tmp(1+counter,:)= pulse;
        out_pulse(1+counter,:) = 0.5*eps0*c0*abs(pulse).^2;
        transverse_grid_short_y(1+counter)=transverse_grid_device_y(1+i);
        counter=counter+1;
end
end

%% Find peak_num distinct peaks
min_time=find(t>3*ps,1); %First index for time>2ps
max_time=find(t>t(end)-3*ps,1); %First index for max>2ps
ind_trans=zeros(counter,20); %indices for peaks
points=-round_trip_time/ps*ones(num_peaks,1); %time values for peaksinitialized
tmp_pulse=out_pulse(:,min_time:max_time); %Temporary pulse variable for manipulation

% Find top peaks, away from edges for all transverse points
for ll=1:counter
    [~,tmp_ind]=findpeaks(tmp_pulse(ll,:),'SortStr','descend');
    ind_trans(ll,1:20)=tmp_ind(1:20);
end

% Out of found peaks, find largest peaks
jk=1;
while jk<=num_peaks
 [row,col]=find(tmp_pulse==max(max(tmp_pulse)));
 tmp_pulse(row,col)=0;
 if min(abs(points-t(col+min_time)))>round_trip_time/1.5
    points(jk)=t(col+min_time);
    jk=jk+1;
 end
end
points=sort(points); %Vector of peak points

%% Plot final output
if IO_finalOutput==1
tmp_fig=figure(111);
    if max(abs(transverse_grid_short_y))~=0
    surf(t/ps, transverse_grid_short_y/um,out_pulse*cm*cm/1e6,'edgecolor','none')
    ylim([min( transverse_grid_short_y),max(transverse_grid_short_y)]/um)
    else
     surf(t/ps, transverse_grid_y/um,out_pulse*cm*cm/1e6,'edgecolor','none')
    ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
    end
    zlim([0,max(max(out_pulse))*cm*cm/1e6]);
  %  ylabel('y [\mum]')
  %  xlabel('t (ps)')
 %   zlabel('I [MW/cm^2]')
    view(77,72.4)
       saveas(tmp_fig,[saveKey,'-QWfield3D.png']);
end


if IO_inversion==1
%% Loading variables: Requires knowledge of lengths [Nt,Nk]
Ne_qw=zeros(Nt,k_ind_max,NUM_TRANSVERSE_DEVICE-2*DEVICE_EXC);
Nh_qw=zeros(Nt,k_ind_max,NUM_TRANSVERSE_DEVICE-2*DEVICE_EXC);
Ne_abs=zeros(Nt,k_ind_max,NUM_TRANSVERSE_DEVICE-2*DEVICE_EXC);
Nh_abs=zeros(Nt,k_ind_max,NUM_TRANSVERSE_DEVICE-2*DEVICE_EXC);

for i = DEVICE_EXC:(NUM_TRANSVERSE_DEVICE-1-DEVICE_EXC)
    %Load each occupation number independently from array
    tmp=loadDA([outKey,num2str(plot_num-1),'_ne_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
    Ne_qw(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);

      tmp=loadDA([outKey,num2str(plot_num-1),'_nh_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
    Nh_qw(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);
    
      tmp=loadDA([outKey,num2str(plot_num-1),'_ne_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
    Ne_abs(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);
    
      tmp=loadDA([outKey,num2str(plot_num-1),'_nh_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
    Nh_abs(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);
end

Inv_qw=Ne_qw+Nh_qw-1;
Inv_abs=Ne_abs-Nh_abs-1;

if IO_inv_surf==1 %transparent surface plot with lines
    for m=1:num_peaks
    tmp_fig=figure(100+m);
    set(tmp_fig,'Position',POS);
    t_init=points(m)-200*fs;
    t_max=points(m)+50*fs;
    t_mid_ind=find(t>=points(m)-10*fs,1);
    t_min_ind=find(t>=t_init,1);
    t_max_ind=find(t>=t_max,1)-1;
    t_diff=t_max_ind-t_min_ind;
    
%     t_init=points(3)-200*fs;
%     t_max=points(3)+50*fs;
%     t_min_ind2=find(t>=t_init,1);
%     t_max_ind2=find(t>=t_max,1)-1;
    
    TRANS_FREQ_SURF=floor((length(transverse_grid_device_y)-2*DEVICE_EXC)/TOTAL_LINES_SURF);
    TRANS_FREQ_INIT=1+floor(TRANS_FREQ_SURF/2); 
    
    %Set particular slice indices
    SURF_SLICE=[t_min_ind,t_max_ind];
    LINES_SLICE=[t_min_ind, t_mid_ind, t_max_ind];
    
    %Set guideline points
    x_line=zeros(TOTAL_LINES_SURF,1);
    line_ind=zeros(TOTAL_LINES_SURF,1);
    line_mesh=zeros(k_ind_max,TOTAL_LINES_SURF);
   
    for j=1:TOTAL_LINES_SURF
        x_line(j)=transverse_grid_device_y(DEVICE_EXC+TRANS_FREQ_SURF*(j-1)+TRANS_FREQ_INIT)/um;
        line_ind(j)=TRANS_FREQ_SURF*(j-1)+TRANS_FREQ_INIT;
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
end

if IO_video==1
    t_init=points(1)-250*fs;
    t_end=points(4)+390*fs;
    t_fast_lhs=250*fs;
    t_fast_rhs=400*fs;
    t_min_ind=find(t>=t_init,1);
    t_max_ind=find(t>=t_end,1);
    [K,X]=meshgrid(k_QW6(1:k_ind_max),transverse_grid_device_y(DEVICE_EXC+1:end-DEVICE_EXC)/um);
    color_axis=[min(min(min(Inv_qw(t_min_ind:t_max_ind,:,:)))) max(max(max(Inv_qw(t_min_ind:t_max_ind,:,:))))];
    myVideo = VideoWriter([saveKey,'INVvideo.avi']);
    myVideo.FrameRate = 60;
  %  open(myVideo);
    fig=figure(100001);
    slice=t_min_ind;
    j=1;
    save_slice=[0;0;0]; %Flag for whether we have saved a figure from this peak already
    
    while t(slice)<t_end
        tmp=squeeze((Inv_qw(slice,:,:)));
        contourf(K,X,transpose(tmp),'edgecolor','none');
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
end

% 4.3 7.3
    %11.5 7.3
    %4.3 12.8
% 11.5 12.8

%1726,1241

