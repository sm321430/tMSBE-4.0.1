%% plotInversion takes output from tMSBE file
%% 03/18/19: v1 created with basic implementaiton
%% 03/20/19: Added final output plotting
%% 06/24/19: v2: Added transparent surface plot for single figure output
%% 07/30/20: VCAV: Modified for twoArmDevice output


clear all
close all

setupConstants
setupPlot

%IOs
IO_finalOutput=0; %Plot final output
IO_inversion=1; %Loads occupation numbers
IO_video=1; %Make and save video to desktop (needs IO_inversion to be on)
IO_inv_surf=0; %Surface plot with guidelines (needs IO_inversion to be on)
IO_higherOrderSBE=1; %Use higher order occupation number terms as well (neds IO_inversion to be on)
IO_abs=1; %Load absorber as well
IO_save=1; %Save plots/videos

DEVICE_EXC=45; %Excludes side devices
max_k=3.5;
XLIM=[-35,35]; %Transverse plotting width (um)
KLIM=[0,max_k];
ZLIM=[-0.2,0.45]; %Height
spc_width=max(abs(XLIM)); %(um)
plot_num = 9; %Output number +1
num_peaks=1;
%videoFile='../../tMSBE_June24.avi';

TOTAL_SLICES_SURF=2; %Number of surface slices
TOTAL_LINES_SURF=5; %Number of lines on a slice, divides domain equally)
TOTAL_LINES_SLICES_SURF=3; %Number of slices of lines
SLICE_FREQ_SURF=710; %Surface frequency plot (currently does nothing)
SLICE_FREQ_LINE=355; %Line frequency plot (currently does nothing)

outKey = '/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV122-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5-seedPulse-2e5-highRes-QWoutput/run/out__';
date='042821';
test='tMSBE-RCAV122-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5-seedPulse-2e5-highRes-QWoutput-narrow';
test_folder='test';
saveKey_local='Fall2020-Summer2021/RingCAV/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey = [saveKey_fold,'/',test,'-out',num2str(plot_num-1),'-QW6-'];
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

%% Load basic cavity parameters 
w0 = loadD([outKey,'w0.dat']);
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])
round_trip_time = loadD([outKey,'round_trip_time.dat']);
transverse_grid_y = loadD([outKey,'transverse_grid_y.dat']);
transverse_grid_device_y = loadD([outKey,'transverse_grid_device_y.dat']);
[~, ind_device_y] = intersect(transverse_grid_y, transverse_grid_device_y);
NUM_TRANSVERSE = length(dir([outKey,num2str(plot_num-1),'_E_re_ABS1_T*.dat']));
NUM_TRANSVERSE_DEVICE = length(dir([outKey,num2str(plot_num-1),'_E_re_ABS1_T*.dat']));

t = loadD([outKey,num2str(plot_num-1),'_t.dat']); %Load time
t_0=t(1);
k_QW6=loadD([outKey,'K_QW6_T',num2str(ind_device_y(1)),'.dat']); %k array: Future versions will load separate k-files 
KLIM(1)=0.85;%min(k_QW6);
FOCUS = 5; % Compensate for lens in output: file output is from SESAM side, which has a spot (w/focus) where w is on the GAIN chip spot
Nt=length(t); %number of time points
t=t-t(1); %Shift time to start at zero
dt=(t(2)-t(1));
%t_init=42.6*ps;     %20.3*ps; %Start time
%t_max=44.0*ps;      %21.8*ps; %Max time
t_freq=2.0*fs; %output slice frequency
t_fast_freq=10.0*fs; %Fast slice frequency
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
    i
    pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
    pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
    pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*w0);
if abs(transverse_grid_device_y(1+i))<spc_width*um %max(abs(transverse_grid_y))/8
        pulse_tmp(1+counter,:)= pulse;
        out_pulse(1+counter,:) = 0.5*eps0*c0*abs(pulse).^2;
        transverse_grid_short_y(1+counter)=transverse_grid_device_y(1+i);
        counter=counter+1;
end
end

%% Find peak_num distinct peaks
min_time=find(t>14.0*ps,1); %First index for time>2ps
max_time=find(t>15.5*ps,1); %First index for max>2ps
ind_trans=zeros(counter,10); %indices for peaks
points=-round_trip_time/ps*ones(num_peaks,1); %time values for peaksinitialized
tmp_pulse=out_pulse(:,min_time:max_time); %Temporary pulse variable for manipulation

% Find top peaks, away from edges for all transverse points
for ll=1:counter
    [~,tmp_ind]=findpeaks(tmp_pulse(ll,:),'SortStr','descend');
    if length(tmp_ind)<10
       tmp_ind(length(tmp_ind)+1:10)=length(tmp_pulse(ll,:)); 
    end
    ind_trans(ll,1:10)=tmp_ind(1:10);
end

% Out of found peaks, find largest peaks
jk=1;
while jk<=num_peaks
 [row,col]=find(tmp_pulse==max(max(tmp_pulse)));
 tmp_pulse(row,col)=0;
 if min(abs(points-t(col+min_time-1)))>2.0e-12%round_trip_time/1.5
    points(jk)=t(col+min_time-1);
    jk=jk+1;
 end
end
points=sort(points); %Vector of peak points


%% Plot final output
if IO_finalOutput==1
    t_cur=t-points(1);
    tmp_fig=figure(111);
    if max(abs(transverse_grid_short_y))~=0
        contourf(t_cur(min_time:max_time)/ps, transverse_grid_short_y/um,out_pulse(:,min_time:max_time)*cm*cm/1e6,'edgecolor','none')
        ylim([min( transverse_grid_short_y),max(transverse_grid_short_y)]/um)
    else
        contourf(t_cur(min_time:max_time)/ps, transverse_grid_y/um,out_pulse(:,min_time:max_time)*cm*cm/1e6,'edgecolor','none')
        ylim([min(transverse_grid_y),max(transverse_grid_y)]/um)
    end
    xlim([min(t_cur(min_time:max_time)),max(t_cur(min_time:max_time))]/ps);
    zlim([0,max(max(out_pulse))*cm*cm/1e6]);
    colormap(map3);
    ylabel('y [\mum]')
    xlabel('t [ps]')
 %   zlabel('I [MW/cm^2]')
  %  view(77,72.4)
    if IO_save==1
       saveas(tmp_fig,[saveKey,'-QWfield3D.png']);
    end
end


if IO_inversion==1
    %% Loading variables: Requires knowledge of lengths [Nt,Nk]
    Ne_qw=zeros(Nt,k_ind_max,NUM_TRANSVERSE_DEVICE-2*DEVICE_EXC);
    Nh_qw=zeros(Nt,k_ind_max,NUM_TRANSVERSE_DEVICE-2*DEVICE_EXC);
    Ne_p2_qw=zeros(Nt,k_ind_max,NUM_TRANSVERSE_DEVICE-2*DEVICE_EXC);
    Nh_p2_qw=zeros(Nt,k_ind_max,NUM_TRANSVERSE_DEVICE-2*DEVICE_EXC);
    Ne_abs=zeros(Nt,k_ind_max,NUM_TRANSVERSE_DEVICE-2*DEVICE_EXC);
    Nh_abs=zeros(Nt,k_ind_max,NUM_TRANSVERSE_DEVICE-2*DEVICE_EXC);

    for i = DEVICE_EXC:(NUM_TRANSVERSE_DEVICE-1-DEVICE_EXC)
        i
        
%         %Load each occupation number independently from array
%         tmp=loadDA([outKey,num2str(plot_num-1),'_ne_00_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
%         Ne_qw(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);
% 
%           tmp=loadDA([outKey,num2str(plot_num-1),'_nh_00_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
%         Nh_qw(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);
%         
%         if IO_higherOrderSBE==1
%             tmp=loadDA([outKey,num2str(plot_num-1),'_ne_p2_re_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
%             %tmp2=loadDA([outKey,num2str(plot_num-1),'_ne_p2_im_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
%             Ne_p2_qw(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);%+1i*tmp2(:,1:k_ind_max);
% 
%             tmp=loadDA([outKey,num2str(plot_num-1),'_nh_p2_re_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
%             %tmp2=loadDA([outKey,num2str(plot_num-1),'_nh_p2_im_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
%             Nh_p2_qw(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);%+1i*tmp2(:,1:k_ind_max);
%         end
        
        if IO_abs==1
              tmp=loadDA([outKey,num2str(plot_num-1),'_ne_00_ABS1_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
            Ne_abs(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);

              tmp=loadDA([outKey,num2str(plot_num-1),'_nh_00_ABS1_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
            Nh_abs(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);
            
            if IO_higherOrderSBE==1
            tmp=loadDA([outKey,num2str(plot_num-1),'_ne_p2_re_ABS1_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
            %tmp2=loadDA([outKey,num2str(plot_num-1),'_ne_p2_im_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
            Ne_p2_abs(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);%+1i*tmp2(:,1:k_ind_max);

            tmp=loadDA([outKey,num2str(plot_num-1),'_nh_p2_re_ABS1_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
            %tmp2=loadDA([outKey,num2str(plot_num-1),'_nh_p2_im_QW6_T',num2str(ind_device_y(1+i)),'.dat'],[Nt,Nk]);
            Nh_p2_abs(:,:,i+1-DEVICE_EXC) = tmp(:,1:k_ind_max);%+1i*tmp2(:,1:k_ind_max);
            end
        end
    end

    %Inv_qw=Ne_qw+Nh_qw+2.0*real(Ne_p2_qw+Nh_p2_qw)-1;
    Inv_abs=Ne_abs+Nh_abs+2.0*real(Ne_p2_abs+Nh_p2_abs)-1;
    Inv_qw=Inv_abs;
    
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
            %    'VerticalAlignment','mid','HorizontalAlignment','left');
            %tmp_z=zlabel('n_e+n_h-1','FontName','Arial','FontSize',30,...
            %    'rotation',90,'VerticalAlignment','baseline','HorizontalAlignment','Center');
            %tmp_y=ylabel('    k(a^{-1}_0)','FontName','Arial','FontSize',28,'rotation',0,...
            %   'HorizontalAlignment','left','VerticalAlignment','mid');
            az=165.7000;
            el=30.8000;
            view(az,el);
            %tmp=title([num2str(m-1),'RT'],'FontSize',20);
            %tmp_pos=get(tmp,'Position');
            %set(tmp,'Position',[tmp_pos(1)+25,tmp_pos(2)-2.5,tmp_pos(3)])
            saveas(tmp_fig,[saveKey,'Inversion-',num2str(m-1),'RT.png']);
       end
    end

    if IO_video==1
         for i = 0:(NUM_TRANSVERSE_DEVICE-1)
                i
                pulse_re = loadD([outKey,num2str(plot_num-1),'_E_re_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
                pulse_im = loadD([outKey,num2str(plot_num-1),'_E_im_ABS1_T',num2str(ind_device_y(1+i)),'.dat']);
                pulse = (pulse_re + 1i*pulse_im).*exp(-1i*t*w0);
                out_pulse(1+i,:) = 0.5*eps0*c0*abs(pulse).^2;
        end
        t_init=points(1)-500*fs;
        t_end=points(num_peaks)+500*fs;
        t_fast_lhs=100*fs;
        t_fast_rhs=150*fs;
        t_min_ind=find(t>=t_init,1);
        t_max_ind=find(t>=t_end,1);
        [K,X]=meshgrid(k_QW6(1:k_ind_max),transverse_grid_device_y(DEVICE_EXC+1:end-DEVICE_EXC)/um);
        color_axis=[-0.3,0.25];%[min(min(min(Inv_qw(t_min_ind:t_max_ind,:,:)))) max(max(max(Inv_qw(t_min_ind:t_max_ind,:,:))))];
        slice=t_min_ind;
        j=1;
        save_slice=[0;0;0]; %Flag for whether we have saved a figure from this peak already
        inversion_slices=zeros([num_peaks*2,size(squeeze((Inv_qw(1,:,:))))]);
        save_counter=1;
        while t(slice)<t_end
            fig=figure(11);
            tmp=squeeze((Inv_qw(slice,:,:)));
            contourf(K,X,transpose(tmp),30,'edgecolor','none');
            ylim(XLIM);
            set(gca,'YTICK',[XLIM(1),0,XLIM(2)]);
            xlim(KLIM);
            set(gca,'FontSize',50)
            zlim(ZLIM);
            title(' ')
            ylabel('y [\mu m]')
            xlabel('k [1/a_0]')
            colormap(map3)
            caxis(color_axis);
            colorbar;
            text(1.0,30,['t=',num2str(floor((t(slice)-points(1))/ps*100)/100),'ps  ']);
            
            axes1=axes('position', [0.69 0.705 0.155 0.2]);
            fig2=plot(transverse_grid_device_y(DEVICE_EXC+1:end-DEVICE_EXC)/um,out_pulse(DEVICE_EXC+1:end-DEVICE_EXC,slice));
           set(gca,'XDir','reverse');
            set(gca,'YTickLabel',[]);
           set(gca,'XTickLabel',[-35,0,35]);
           ylim([0,max(max(out_pulse))]);
            xlim(XLIM);
            camroll(-90)
            
            
            %set(gca,'Xdir','reverse');
            %xlabel('k [1/a_0]');
            %ylabel(' x [\mu m] ');
            %zlabel('n_e+n_h-1');
            
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
            inversion_slices(save_counter,:,:)=squeeze((Inv_qw(slice,:,:)));
            save_counter=save_counter+1;
            end
            if t(slice)>=points(ind) && save_slice(2)==0 %Peak
            saveas(gcf,[saveKey,'InversionPeak-',num2str(ind-1),'RT.png']);
            save_slice(2)=1;
            end
            if t(slice)>=points(ind)+100*fs && save_slice(3)==0 %100fs after peak
            saveas(gcf,[saveKey,'InversionPostPeak-',num2str(ind-1),'RT.png']);
            save_slice(3)=1;
            inversion_slices(save_counter,:,:)=squeeze((Inv_qw(slice,:,:)));
            save_counter=save_counter+1;
            end
            slice=slice+SLICE_FREQ;
        else
            slice=slice+FAST_SLICE_FREQ;
            save_slice=[0;0;0];
        end
        %title(['QW Inversion: t=',num2str(floor(t_0/ns)),'ns+',num2str(floor((t(slice)-t(1))/ps*100)/100),'ps  ']);
        drawnow

        %Video creation/saving tools
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        ax.Units = 'normalized';
        M(j)=getframe(fig);
        %Update counter
        j=j+1;
        close all;
        end
        
        %movie(fig,M);
        
        if IO_save==1
            myVideo = VideoWriter([saveKey,'INVvideo.avi']);
            myVideo.FrameRate = 5;
            open(myVideo);
            writeVideo(myVideo, M);
            close(myVideo);
        end
        
        for j=1:num_peaks
          tmp=figure;
          hole_depth=transpose(squeeze(inversion_slices(2*j,:,:))-squeeze(inversion_slices(2*j-1,:,:)));
          contourf(K,X,hole_depth,'edgecolor','none'); 
          ylim(XLIM);
          set(gca,'YTICK',[XLIM(1),0,XLIM(2)]);
          xlim(KLIM);
          xlabel('k [1/a_0]');
          ylabel('n_k^e+n_k^h-1\mid_{t_0}^{t_1}');
          caxis([-0.175,0.05]);
          colormap(map3);
          colorbar;
          if IO_save==1
            saveas(tmp,[saveKey,'inversionDifference',num2str(j-1),'RT.png']);
          end
          total_hole_depth(j)=trapz(trapz(hole_depth));
          hole_depth_x(j,1:k_ind_max)=trapz(transverse_grid_device_y(DEVICE_EXC+1:end-DEVICE_EXC)/um,hole_depth);
          %diff_integrated_k(j,1:length(transverse_grid_device_y(DEVICE_EXC+1:end-DEVICE_EXC)))=trapz(transpose(diff));
        end
        tmp=figure;
        hold on
        load hole_depth_x_VCAV86;
        plot(k_QW6(1:k_ind_max),hole_depth_x_VCAV86,'DisplayName','6.7mm, 1:1');
        for j=1:num_peaks
        plot(k_QW6(1:k_ind_max),hole_depth_x(j,:),'DisplayName',['6.7mm, 7:1, t_', num2str(j)])
        end
        xlim(KLIM);
        grid on
        xlabel('k [1/a_0]');
        ylabel('Kinetic Hole Depth [-]');
        if num_peaks>1
           legd=legend('show');
           set(legd,'Location','Best');
        end
        if IO_save==1
            saveas(tmp,[saveKey,'integratedInversion_x.png']);
        end
        
%         tmp=figure;
%         hold on
%         for j=1:num_peaks
%         plot(transverse_grid_device_y(DEVICE_EXC+1:end-DEVICE_EXC)/um,diff_integrated_k(j,:)/max(max(diff_integrated_k)),'DisplayName',['Pass ', num2str(j)])
%         end
%         xlim(XLIM);
%         xlabel('x [1/\mu m]');
%         ylabel('\int_{k}N(t_0)-N(t_1)dk [a.u.]');
%         legend show
%         if IO_save==1
%             saveas(tmp,[saveKey,'integratedInversion_k.png']);
%         end
        
    end
end

