% Coding up basic diffusion
%df/dt=-Dd^2f/dt^2
close all
clear all

%% IO flags
IO_sourceContour=0; %Source term contour
IO_densityPlot=0; %Plot overlay of output
IO_densityContour=1; %Contour of output
IO_densityVideo=1; %Video of output
IO_break=0; %Break after source plot
IO_save=0; %Does nothing currently

%% Constants and Setup
setupPlot
setupConstants
ns=1.0e-9;
ps=1.0e-12;
fs=1.0e-15;
um=1.0e-6;

%% Space intialization
Nx=288;
XLIM=2.5*336*um;
x=linspace(-XLIM,XLIM,Nx);
dx=x(2)-x(1);

%% Time intiialization
dt=10.0*fs;
dt_skip=10;
t=0:dt:5*22.1*ps;
Nt=length(t);

%% Physical Constants
D=0.01;%0.01;
tau=100.0*ps;%1.0*ns; %Recombination rate
theta=1.5*pi/180; %Angle for interference fringes
lambda=1030*nm; %Wavelength
N_0=2.0e16; %Peak pump density
DiffLength=sqrt(D*tau); %Diffraction length
f_0=exp(-(x./(336*um)).^12);%-0.05.*exp(-(x./(201.6*um)).^12).*(1+0.5*(cos(2*pi.*theta*x/lambda))); %Initial condition
f_0=f_0*N_0; %Initial conditions


%g=N_0*exp(-(x./(336*um)).^12)'; %Pump/Source conditions (Single Gaussian)

%Pulse hitting gain chip on successive passses
pulse_time= @(time) exp(-(mod(time,22.1*ps)-11.05*ps).^2/(0.2*ps).^2);
pulse_space= @(space) exp(-(space./(201.6*um)).^2).*(1+0.5*(cos(2*pi.*theta*space/lambda)));
decayed_inf= @(time) 5.0*exp(-time/(max(t)));
g = @(time) f_0';
pulse= @(time) 7.0E25*decayed_inf(time)*pulse_time(time)*pulse_space(x)';

source=zeros(Nt,Nx);
for j=1:Nt
   source(j,:)=g(t(j));
end

if IO_sourceContour==1
    tmp_fig=figure(101);
    contourf(t/ps,x/um,source');
    ylim([min(x)/2,max(x)/2]/um);
    xlabel('t [ps]');
    ylabel('x [um]');
    colorbar;
end

if IO_break==1
   ABRUPTEND 
end


%% Zeroing Initialization
f_now=zeros(Nx,Nt);
k1=zeros(Nx,1);
k2=k1;
k3=k1;
k4=k1;
f_now(:,1)=f_0;

%% Stencil
M = diag((-2)*ones(1,Nx)) + diag(ones(1,Nx-1),1) + diag(ones(1,Nx-1),-1);

%% PDE
dfdt=@(D,y,dx,tau,g, pulse) (D/dx^2).*M*y-(y-g)/tau-pulse;


%% RK4
for j=2:length(t)
        k1=dfdt(D,f_now(:,j-1),dx,tau,g(t(j-1)),pulse(t(j-1)));
        k2=dfdt(D,f_now(:,j-1)+k1*dt/2,dx,tau,g(t(j-1)+dt/2),pulse(t(j-1)+dt/2));
        k3=dfdt(D,f_now(:,j-1)+k2*dt/2,dx,tau,g(t(j-1)+dt/2),pulse(t(j-1)+dt/2));
        k4=dfdt(D,f_now(:,j-1)+k3*dt,dx,tau,g(t(j)),pulse(t(j)));
        f_now(:,j)=f_now(:,j-1)+dt*(k1+2*k2+2*k3+k4);
        j
end

%% VISUALIZATION
if IO_densityContour==1
    tmp_fig=figure(101);
    set(tmp_fig,'Name','Density Contour');
    contourf(t(1:dt_skip:end)/ps,x/um,f_now(:,1:dt_skip:end),'Edgecolor','none')
    ylim([min(x)/2,max(x)/2]/um)
    xlabel('t [ps]');
    ylabel('x [\mu m]');
    colorbar;
end

if IO_densityPlot==1
    tmp_fig=figure(1010);
    set(tmp_fig,'Name','Density Plot');
    hold on
    counter=0;
    for j=1:max(1,floor(Nt/6)):Nt
        counter=counter+1;
       plot(x/um,f_now(:,j));
       lgd_entry{counter}=['time=',num2str(t(j)/ps),'ps'];
    end
    lgd=legend(lgd_entry);
    xlim([min(x)/2,max(x)/2]/um)
    ylabel(' Density [a.u.]');
    xlabel('x [\mu m]');
    hold off
end

if IO_densityVideo==1
    tmp_fig=figure(1010);
    set(tmp_fig,'Name','Density Plot');
    
    counter=0;
    for j=1:100:Nt
        plot(x/um,f_now(:,j));
        xlim([min(x)/2,max(x)/2]/um)
        ylabel(' Density [a.u.]');
        xlabel('x [\mu m]');
        title(['time=',num2str(t(j)/ps),'ps']);
        drawnow
        counter=counter+1;
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        ax.Units = 'normalized';
        F(counter)=getframe(gcf);
    end
    %movie(tmp_fig,F);
end