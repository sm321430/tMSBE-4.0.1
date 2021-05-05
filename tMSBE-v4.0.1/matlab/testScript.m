%VARIOUS TESTING FUNCTIONS
% setupConstants
% 
% F=10;
% w0=(20.16*10/F)*um/sqrt(2*log(2));
% lambda=1030*nm;
% zR=pi*w0^2/lambda;
% z=sqrt(F^2-1)*zR;
% format long
% f=@(tau,mag) sqrt(-2*log(mag)*(tau^2)); %Compute time Gaussian decay to mag with width tau

clear all
setupConstants
shortCav_normalLength=3200.5*1030*nm;
longCav_normalLength=6400.5*1030*nm;
theta4_totalLength=5*336*um*tand(4);
shortCav_shortArm_normalLength=((3200/8)+0.25)*1030*nm;
longCav_shortArm_normalLength=((6400/8)+0.25)*1030*nm;
shortCav_longArm_normalLength=((3200*7/8)+0.25)*1030*nm;
longCav_longArm_normalLength=((6400*7/8)+0.25)*1030*nm;

shortCav_totalLength=shortCav_normalLength+theta4_totalLength;
longCav_totalLength=longCav_normalLength+theta4_totalLength;
shortCav_shortArm_averageLength=shortCav_shortArm_normalLength+theta4_totalLength/2;
longCav_shortArm_averageLength=longCav_shortArm_normalLength+theta4_totalLength/2;
shortCav_longArm_averageLength=shortCav_longArm_normalLength+theta4_totalLength/2;
longCav_longArm_averageLength=longCav_longArm_normalLength+theta4_totalLength/2;

shortCav_totalTransit=shortCav_totalLength*2/c0;
longCav_totalTransit=longCav_totalLength*2/c0;
shortCav_shortArm_averageTransit=shortCav_shortArm_averageLength*2/c0;
longCav_shortArm_averageTransit=longCav_shortArm_averageLength*2/c0;
shortCav_longArm_averageTransit=shortCav_longArm_averageLength*2/c0;
longCav_longArm_averageTransit=longCav_longArm_averageLength*2/c0;

shortCav_totalTransit_summed=shortCav_longArm_averageTransit+shortCav_shortArm_averageTransit;
longCav_totalTransit_summed=longCav_longArm_averageTransit+longCav_shortArm_averageTransit;
shortCav_totalTransit_nominal=2*3.4e-3/c0;
longCav_totalTransit_nominal=2*6.7e-3/c0;