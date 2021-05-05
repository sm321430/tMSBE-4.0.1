function plotStructure()
%% Plot structure from structure file
close all
clear all

global nm;
global um;

global w0;
global c0;
global QW_THICKNESS;

um = 1.0e-6;
nm = 1.0e-9;

c0 = 2.99792458E+08;
e = 1.60217657E-19;
hbar = 1.054571726E-34;
QW_THICKNESS = 8.35*nm;

W = what;
addpath(W.path)
addpath([W.path,'/','matlab2tikz/src'])
addpath([W.path,'/','matlab2tikz/tools'])


outKey ='/Volumes/SAMbackup/tMSBE-VCAV-data-2020/tMSBE-v3.7-RCAV3-1D-n2p5-colThresh-2em2-theta4-6400lam-spontEmis-wExpSBE/run/';
date='113020';
test='tMSBE-v3.7-RCAV3-1D-n2p5-colThresh-2em2-theta4-6400lam-spontEmis-wExpSBE-V';
runKey = [outKey,'refSpecQW__'];
%outKey = '../run/refSpecQW__';
location='CAVOC'; %Field location for uploading and saving
point='T0';
plot_num =1; %Output number +1
test_folder='test';
saveKey_local='Fall2020-Summer2021/RingCav/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey = [saveKey_fold,'/',test,'-out',num2str(plot_num-1),'-',location,'-'];
user_entry = input(['saveKey=',saveKey,'? (y to continue): '], 's');
if user_entry~='y'
    gafuggle
end
if exist(saveKey_fold,'dir')~=7
    mkdir(saveKey_fold);
end


w0 = loadD([runKey,'w0.dat']);
disp(['Load: w0 = ',num2str(w0*hbar/e,'%.3f'),' [eV]'])
disp(['lambda   = ',num2str(2*pi*c0/w0/nm,'%.2f'),' [nm]'])

tmp_fig=figure(1);
[ni, z0, z1,qw_ind] = loadStruct([runKey,'system_structure.dat']);
[ht,XLIM] = plotStruct(ni,z0,z1,qw_ind,w0);

[XTOT, YTOT] = plotSingleFreqAll(w0,ni,z0,z1,3.63,1.0,'L',qw_ind);

for j=1:length(ni)
   disp([num2str(ni(j),'%3f')]);
end
hold on
plot(XTOT,abs(YTOT).^2,'b-')
hold off
title('Active region')
setXlim(XLIM,0,qw_ind);
ylim([0,4])
ylabel('n')

saveas(tmp_fig,[saveKey,'VCSELStruct.png']);
%cleanfigure();
%matlab2tikz('filename','plot_design_vcsel.tikz','height', '\figureheight', 'width', '\figurewidth');

figure(2)
[ni, z0, z1,qw_ind] = loadStruct([outKey,'refSpecABS__system_structure.dat']);
[ht,XLIM] = plotStruct(ni,z0,z1,qw_ind,w0);
[XTOT, YTOT] = plotSingleFreqAll(w0,ni,z0,z1,1.0,3.63,'R',qw_ind);
hold on
plot(XTOT,abs(YTOT).^2,'b-')
hold off
title('Arm1: SESAM')
setXlim(XLIM,z1(end),qw_ind);
ylim([0,4])
ylabel('n')
saveas(tmp_fig,[saveKey,'SESAMyStruct.png']);
% oldDr = cd('../pics');
% saveas(gcf,'plot_design_arm1_sesam','eps2c')
% cleanfigure();
% matlab2tikz('filename','plot_design_arm1_sesam.tikz','height', '\figureheight', 'width', '\figurewidth');
% cd(oldDr);


figure(3)
[ni, z0, z1,qw_ind] = loadStruct('outV2__system_structure.dat');
ht = plotStruct(ni,z0,z1,qw_ind,w0);
title('Arm2: Output')
setXlim([],z1(end),qw_ind);
ylim([0,4])
ylabel('n')

oldDr = cd('../pics');
saveas(gcf,'plot_design_arm2_output','eps2c')
cleanfigure();
matlab2tikz('filename','plot_design_arm2_output.tikz','height', '\figureheight', 'width', '\figurewidth');
cd(oldDr)

asd


ht = figure(1);
plot(XTOT/um,abs(YTOT)/normLize,'b-')
%     XTOT/um,atan2(imag(YTOT),real(YTOT))/pi,'r--')
hold off
xlabel('x [um]')
uistack(ht, 'down')

%xlim([0,0.6])      % DBR: RHS
%xlim([0,2])         % VECSEL: ALL
%xlim([5.5,8])       % VECSEL: GAIN
%xlim([3217.3, 3217.8])  % ABS: ALL

%ylim([0, 4])
%ylim([0,1e-6])

h4 = figure(2);
copyobj(get(ht,'children'),h4);
xlim((z1(AIR_IND(1))+1*um*[-1,1])/um)
title('PROBE')


W = what;
addpath([W.path,'/','matlab2tikz/src'])
addpath([W.path,'/','matlab2tikz/tools'])

h2 = figure(3);
copyobj(get(ht,'children'),h2);
xlim(LEFT_DEVICE_PLOT_RANGE/um)
title('Active Mirror')
ylim([0,4])
ylabel('n')
xlim([5.779,7.6])

dirName = '../pics/';

if ~exist(dirName,'dir')
    mkdir(dirName)
end

oldDirPic = cd(dirName);
cleanfigure();
saveas(gca,'designActiveMirror','eps2c')
matlab2tikz('filename','designActiveMirror.tikz','height', '\figureheight', 'width', '\figurewidth');
cd(oldDirPic)

h3 = figure(4);
copyobj(get(ht,'children'),h3);
xlim(RIGHT_DEVICE_PLOT_RANGE/um)
title('SESAM')
ylim([0,4])
xlim([3212.7, 3213.16])
ylabel('n')

oldDirPic = cd(dirName);
cleanfigure();
saveas(gca,'designSESAM','eps2c')
matlab2tikz('filename','designSESAM.tikz','height', '\figureheight', 'width', '\figurewidth');
cd(oldDirPic);


cd(oldDir)

end

function setXlim(XLIM,offset,qw_ind)

global um;

h = findobj(gca,'Type','line');
old = get(h,'xdata');
new = cellfun(@(x)(x-offset)/um,old,'UniformOutput',false);
set(h,{'xdata'}, new);


    function x = tmpFun(x)
        
        x(1) = (x(1)-offset)/um;
        x(3) = x(3)/um;
    end

if (numel(qw_ind)>0)
    
    h = findobj(gca,'Type','rectangle');
    old = get(h,'Position');
    if (iscell(old))
        new = cellfun(@tmpFun,old,'UniformOutput',false);
        set(h,{'Position'}, new);
    else
        new = old; new(1) = (new(1) - offset)/um; new(3) = new(3)/um;
        set(h,'Position', new);
    end
    
end

if (numel(XLIM)>0)
    xlim((XLIM-offset)/um)
end

xlabel('x [um]')

end

function v = loadD(name)
% Read a single double or a list of doubles into v
% v is a row vector

fid = fopen(name,'rb');
v = fread(fid,'double');
fclose(fid);

end


function [ht,XLIM] = plotStruct(ni,z0,z1,qw_ind,w0)

global c0;

% Replace all QW indices with the correct one
ni(qw_ind) = real(getQWind(2*pi*c0./w0));

%ht = figure;
ht = 0;

hold on

%% Plot all edges

for i = 1:numel(z0)
   
    
    x0 = z0(i);
    x1 = z1(i);
    h = ni(i);
    
    plot([x0,x0],[0,h],'k','LineWidth',1);
    plot([x0,x1],[h,h],'k','LineWidth',1);
    plot([x1,x1],[0,h],'k','LineWidth',1);
    %uistack(h2, 'bottom')
    
end


%% Plot QW's

for i = 1:numel(qw_ind)
        
   x0 = z0(qw_ind(i));
   x1 = z1(qw_ind(i));
   qw_n = ni(qw_ind(i));
   heigth = qw_n;
   
   % Barrier around single QW
   H = rectangle('Position',[x0,0,x1-x0,heigth],'FaceColor',0.9*[1,1,1],'EdgeColor','none');
   uistack(H,'bottom');
    
end


hold off

%xlabel('x [um]')
%title('Active region')

%{
% Plot only QWs with it s own gain region
if (numel(qw_ind)>0)
   xmin = z0(qw_ind(1)  -1);
   xmax = z1(qw_ind(end)+1);
   xlim([xmin, xmax]/um)
end
%}

% Plot QW region pluss an extra bit
if (numel(qw_ind)>0)
    
    ind_air = find(ni==1);
    lambda = 2*pi*c0/w0;
    if (z0(1) > 0)
        
        xmin = z1(ind_air) - lambda/4;
        xmax = z1(end-1);
        
        if (qw_ind(end)+2 > length(z1))
            xmax = z1(end);
        else
            xmax = (z0(qw_ind(end)  + 2)+z1(qw_ind(end)  + 2))/2;
        end
    else
        
        xmax = z0(ind_air) + lambda/4;
        
        if (qw_ind(1)-2 > length(z0))
            xmin = z0(1);
        else
            xmin = (z0(qw_ind(1)  -2)+z1(qw_ind(1)  - 2))/2;
        end
    end
    
    %xlim([xmin, xmax])
    
else
    ind_air = find(ni==1);
    if (numel(ind_air)>=numel(ni)/2)
        xmin = z0(1);
        xmax = z1(end-numel(ind_air)+1);
    else
        
        lambda = 2*pi*c0/w0;
        if (z0(1) > 0)
            xmin = z1(ind_air) - lambda/4;
            xmax = z1(end);
        else
            xmin = z0(1);
            xmax = z0(ind_air) + lambda/4;
        end
    end
    
    %xlim([xmin, xmax])
end

XLIM = [xmin,xmax];

end


function [n_s, z0_s, z1_s, qw_ind] = loadStruct(fileName)
%% Plot structure in background
% Input:
% runKey - What to load from
% Output:
% n_s  - Refractive index in layer s 
% z0_s - Position of Left edge of layer s
% z1_s - Position of Right edge of layer s
% qw_pos - Position of QW's

str = load(fileName);
[n,m] = size(str);

% Import refractive indices and positions
z0_s = [];
z1_s = [];
n_s = [];
qw_z0_ind = [];
cnt = 1;
cnt2 = 1;
for i = 1:n
    
        
    if (str(i,1)==2)    % QW
        qw_z0_ind(cnt2) = cnt;
        cnt2 = cnt2 + 1;
        
    end
    
    if (str(i,1) == 1)  % CAVITY
        
        z0_s(cnt) = str(i,2);
        z1_s(cnt) = str(i,3);
        n_s(cnt)  = str(i,4);
        cnt = cnt + 1;
    end
end

%Find and remove air layers, replace them with 2um of air
ind = find(n_s==1);
%if (numel(ind) <= length(n_s)/2)
    z0_s(ind) = [];
    z1_s(ind) = [];
    n_s(ind) = [];

    for i = 1:length(qw_z0_ind)

        % Find how many layers before qw_z0_ind(i) is deleted
        num = sum(ind<qw_z0_ind(i));
        qw_z0_ind(i) = qw_z0_ind(i) - num;

    end


    AIR_WIDTH = 1e-6;
    if (z0_s(1)>AIR_WIDTH)

        z1_s = [z0_s(1)          , z1_s];
        z0_s = [z0_s(1)-AIR_WIDTH, z0_s];
        n_s  = [1, n_s];

        % All layers are infront of the first one
        qw_z0_ind = qw_z0_ind + 1;

    else
        z0_s = [z0_s, z1_s(end)];
        z1_s = [z1_s, z0_s(end)+AIR_WIDTH];
        n_s  = [n_s,1];
    end
%end

width = z1_s - z0_s;


%===========================
%       Active QW's
%===========================
% Modify the index layers that have QW's in them

global N_QW_w;
global N_QW;
global nm;
global QW_THICKNESS;

%%SAM-hardcode 128 for number of transverse points
NUM_TRANS=128;
for i = 1:numel(qw_z0_ind)/NUM_TRANS
   
   ind = qw_z0_ind(i);
   
   % Modyfy lengths
   z1_s(ind-1) = z1_s(ind-1) - QW_THICKNESS/2;
   z0_s(ind)   = z0_s(ind)   + QW_THICKNESS/2;
   
   layer_z0 = z1_s(ind-1);
   layer_z1 = z0_s(ind);
   
   % Add in new layer
   n_s  = [ n_s(1:ind-1),-1, n_s(ind:end)];
   z0_s = [z0_s(1:ind-1),layer_z0,z0_s(ind:end)];
   z1_s = [z1_s(1:ind-1),layer_z1,z1_s(ind:end)];
   
    
   % Add one to the index becuase the layer is extended
   qw_z0_ind(i:end) = qw_z0_ind(i:end)+1;
    
end

width = z1_s - z0_s;
terrible_ind = find(width <= 0);


if numel(terrible_ind)>0
    
    disp('----------------------')
    disp(' QWs are overlapping..')
    disp(' Ignoring for now')
    disp('----------------------')
end



qw_ind = find(n_s==-1);

% Load active QW indices
oldDir = cd('../matlab/index');
[N_QW_w , N_QW ] = loadIndexTable('singleQWindex.txt');
%[N_ABS_w, N_ABS] = loadIndexTable('singleABSindex.dat');
cd(oldDir);


end

function [w,n] = loadIndexTable(name)


%dir
%oldDir = cd('../pics');

data = load(name);
w = data(:,1);
n = data(:,2) + 1i*data(:,3);

%cd(oldDir);

end

function [nl,lambda] = getQWind(lambda)
% Return complex index of refraction

global N_QW;
global N_QW_w;
global c0;

l = 2*pi*c0./N_QW_w;

if ((max(lambda) > max(l))||(min(lambda) < min(l)))
    disp('getQWind(): Requesting lambda outside of region')
    disp(['Asking for l = ',num2str(lambda/1.0e-9),' [nm]'])
    disp(['Max Lambda = ',num2str(max(l)/1.0e-9),' [nm]'])
    disp(['Min Lambda = ',num2str(min(l)/1.0e-9),' [nm]'])
    asd
end

nl = interp1(l, N_QW,lambda);

end


function [XTOT,YTOT] = plotSingleFreqAll(w0,ni,z0,z1,na,nb,start,qw_ind)

global c0;

% Replace all QW indices with the correct one
ni(qw_ind) = real(getQWind(2*pi*c0./w0));

width = z1-z0;

coff = zeros(2,length(z0)+1);

if (start=='L')
    
    %% Frequency
    % First layer, on right side of layer
    rho  = (ni(1)-na)/(na + ni(1));
    t    =    2*ni(1)/(na + ni(1));
    coff(1,1) = rho/t;
    coff(2,1) = 1/t;
    
    %Fill in coff
    for i = 2:length(z0)
        n_left  = ni(i-1);
        n_right = ni(i);
        w_left  = width(i-1);
        coff(:,i) =  setForwardCoff(w0,n_left,n_right,w_left,coff(:,i-1));
    end
    
    %First layer
    rho = (nb-ni(end))/(nb + ni(end));
    t = 2*nb/(nb + ni(end));
    coff(:,end) = ([1  , rho;
        rho,   1]/t)*coff(:,end-1);
    
    
elseif (start == 'R')
    
    %Last layer
    rho = (ni(end)-na)/(na + ni(end));
    t = 2*ni(end)/(na + ni(end));
    coff(1,end) = 1/t;
    coff(2,end) = rho/t;
    
    %Fill in coff
    for i = length(z0):-1:2
        n_left  = ni(i-1);
        n_right = ni(i);
        w_right  = width(i);
        coff(:,i) =  setBackwardsCoff(w0,n_left,n_right,w_right,coff(:,i+1));
    end
    
    rho  = (nb-ni(1))/(nb + ni(1));
    t    =    2*nb/(nb + ni(1));
    coff(:,1) = ([1  , rho;
                  rho,   1]/t)*coff(:,2);
    
else
    
    disp('No option for plotting the standing wave this way...')
    asd
end


%nrm = norm(coff(:,qw_ind(1)));
%ind_air = find(ni==1);
%nrm = sqrt(2)*norm(coff(:,ind_air(1)));
%coff = coff/nrm;


XTOT = [];
YTOT = [];
normLize = 0;
for i = 1:length(z0)

    
    if (start=='L')
        ind_air = find(ni==1);
        nrm = sqrt(2)*norm(coff(:,ind_air(1)));
        coff = coff/nrm;
        
        Npoints = numPlotPoints(z1(i)-z0(i),w0,ni(i));
        x = linspace(z0(i), z1(i),Npoints);
        gn = w0*ni(i)/c0;
        
        y = coff(1,i)*exp(-1i*gn*(x-z0(i))) + coff(2,i)*exp(1i*gn*(x-z0(i)));
    elseif (start=='R')
        
        ind_air = find(ni==1);
        nrm = sqrt(2)*norm(coff(:,ind_air(1)+1));
        coff = coff/nrm;
        
        Npoints = numPlotPoints(z1(i)-z0(i),w0,ni(i));
        x = linspace(z0(i), z1(i),Npoints);
        gn = w0*ni(i)/c0;
        y = coff(1,i+1)*exp(-1i*gn*(x-z1(i))) + coff(2,i+1)*exp(1i*gn*(x-z1(i)));
    else
        disp('No option for plotting the standing wave this way...')
        asd
    end
    
    XTOT = [XTOT, x];
    YTOT = [YTOT, y];
    
end

YTOT = YTOT;
XTOT = XTOT;

end


function N = numPlotPoints(L,w0,n)

global c0;

numLambda = L/(2*pi*c0/(w0*n));

if (numLambda <= 0.25)
    N = 25;
elseif (numLambda < 1)
    N = 100;
else
    N = 100*floor(numLambda);
end


end


function coff1 = setForwardCoff(wF,n_left,n_right,width_left,coff0)
%% Set Coefficients for stationary wave analysis

% Removed global for constant
c0 = 2.99792458E+08;

rho = (n_right-n_left)/(n_right + n_left);
t = 2*n_right/(n_right + n_left);
arg_ikl = 1i*width_left*(wF/c0)*n_left;

%{
ema = exp(-arg_ikl);
epa = exp( arg_ikl);

% Transfer matrix through a linear medium moving from left side to right
% side. The other direction is the inverse

T = [ema, 0;
     0,  epa];
% Matching matrix for boundary condition, from left to right
M = [1, rho;
     rho, 1]/t;

% Match layer then transfer
coff1 = M*T*coff0;
%}

ema = exp(-arg_ikl)*coff0(1);
epa = exp( arg_ikl)*coff0(2);

coff1 = coff0;
coff1(1) = (    ema + rho*epa)/t;
coff1(2) = (rho*ema +     epa)/t;

end

function coff1 = setBackwardsCoff(wF,n_left,n_right,width_right,coff0)
%% Set Coefficients for stationary wave analysis

% Removed global for constant
c0 = 2.99792458E+08;

rho = (n_left-n_right)/(n_right + n_left);
t = 2*n_left/(n_right + n_left);
arg_ikl = 1i*width_right*(wF/c0)*n_right;

%{
ema = exp(-arg_ikl);
epa = exp( arg_ikl);

% Transfer matrix through a linear medium moving from right side to left
% side. The other direction is the inverse

T = [epa, 0;
     0,  ema];
% Matching matrix for boundary condition, from right to left
M = [1, rho;
     rho, 1]/t;

% Transfer then match layer
coff1 = M*T*coff0;
%}

epa = exp( arg_ikl)*coff0(1);
ema = exp(-arg_ikl)*coff0(2);

coff1 = coff0;
coff1(1) = (    epa + rho*ema)/t;
coff1(2) = (rho*epa +     ema)/t;

end




