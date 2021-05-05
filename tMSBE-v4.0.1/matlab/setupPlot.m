maps %Load maps file for color scheme
POS=[1,1,1200,800];
POS2=[1,1,500,900];
POS_pump=[1,1,1300,900];
my_lineStyle={'-k',':k','-b',':b','-r',':r', '-g', ':g', '-c', ':c', '-m', ':m', '-y',':y',':r', ':g', ':c',':m','-kx','-bx','-rx', '-gx', '-cx','-mx','-k+','-kx','-k*','-b+','-bx','-b*','-g+'};
my_plot_color={[0         0    1.0000]
    [1.0000         0         0]
    [    0    1.0000         0]
    [     0         0    0.1724]
    [1.0000    0.1034    0.7241]
    [1.0000    0.8276         0]
    [     0    0.3448         0]
    [0.5172    0.5172    1.0000]};
my_plot_style={':','-'};


%% Fat lines and text
% set(0,'defaulttextinterpreter','tex') %Default
% %set(0,'defaulttextinterpreter','latex')
% set(0,'defaultAxesFontName', 'Arial')
% set(0,'defaultAxesFontsize', 60)
% set(0,'defaultAxesLinewidth', 8)
% set(0,'defaultTextFontsize', 60)
% set(0,'defaultlineMarkerSize',50)
% set(0,'defaultlinelinewidth',10) %Thin lines
% set(0,'defaultfigureposition',POS)

%% Smaller text
% set(0,'defaulttextinterpreter','tex') %Default
% %set(0,'defaulttextinterpreter','latex')
% set(0,'defaultAxesFontName', 'Arial')
% set(0,'defaultAxesFontsize', 30)
% set(0,'defaultAxesLinewidth', 4)
% set(0,'defaultTextFontsize', 30)
% set(0,'defaultlineMarkerSize',20)
% set(0,'defaultlinelinewidth',7) %Thin lines
% set(0,'defaultfigureposition',POS)

%% Original 
set(0,'defaulttextinterpreter','tex') %Default
%set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontsize', 35) 
set(0,'defaultAxesLinewidth', 4)
set(0,'defaultTextFontsize', 45) 
set(0,'defaultlineMarkerSize',30)
set(0,'defaultlinelinewidth',7) %Thin lines
set(0,'defaultfigureposition',POS)