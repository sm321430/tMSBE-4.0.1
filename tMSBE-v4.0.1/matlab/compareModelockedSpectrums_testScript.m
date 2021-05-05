%%S.A.McLaren
%%Loads set of modelocked simulation data that has been preprocessed by plotTransverseFieldOutput_v2
%%Plots an overlay of computed spectrums
%Last updated 12/08/20
close all
clear all
% 
% key(1)={'tMSBE-v3.7-VCAV64-equalArms-n2p0-fullThresh-2em2-Ny288-theta2-out14-OUTPUT-pulseSpectrum1.mat'};
% key(2)={'tMSBE-v3.7-VCAV63-equalArms-n2p0-fullThresh-2em2-Ny288-theta4-out13-OUTPUT-pulseSpectrum1.mat'};
% key(3)={'tMSBE-v3.7-VCAV68-equalArms-n2p0-fullThresh-2em2-Ny288-theta8-out11-OUTPUT-pulseSpectrum1.mat'};
% key(4)={'tMSBE-v3.7-VCAV75-equalArms-n2p0-fullThresh-2em2-Ny288-theta12-out12-OUTPUT-pulseSpectrum1.mat'};
% key(5)={'tMSBE-v3.7-VCAV69-equalArms-n2p0-fullThresh-2em2-Ny288-theta16-out11-OUTPUT-pulseSpectrum1.mat'};
% key(6)={'tMSBE-v3.7-VCAV76-equalArms-n2p0-fullThresh-2em2-Ny288-theta20-out15-OUTPUT-pulseSpectrum1.mat'};
% 
% savedData(1).legend={'\theta=2^\circ'};
% savedData(2).legend={'\theta=4^\circ'};
% savedData(3).legend={'\theta=8^\circ'};
% savedData(4).legend={'\theta=12^\circ'};
% savedData(5).legend={'\theta=16^\circ'};
% savedData(6).legend={'\theta=20^\circ'};

% key(1)={'tMSBE-v3.7-VCAV63-equalArms-n2p0-fullThresh-2em2-Ny288-theta4-out13-OUTPUT-pulseSpectrum1.mat'};
% key(2)={'tMSBE-v3.7-VCAV67-Arms7to1-n2p0-fullThresh-2em2-Ny288-theta4-out12-OUTPUT-pulseSpectrum1.mat'};
% key(3)={'tMSBE-v3.7-VCAV86-equalArms-6400lambda-n2p0-fullThresh-2em2-Ny288-theta4-out11-OUTPUT-pulseSpectrum1.mat'};
% key(4)={'tMSBE-v3.7-VCAV85-Arms7to1-6400lambda-n2p0-fullThresh-2em2-Ny288-theta4-out12-OUTPUT-pulseSpectrum1.mat'};
% % 
% 
% savedData(1).legend={'3.4mm, 1:1'};
% savedData(2).legend={'3.4mm, 7:1'};
% savedData(3).legend={'6.7mm, 1:1'};
% savedData(4).legend={'6.7mm, 7:1'};

% key(1)={'tMSBE-RCAV5-1D-n2p5-colThresh-2em2-theta2-6400lam-spontEmis-noExpSBE-restarted-out47-OUTPUT-pulseSpectrum1.mat'};
% key(2)={'tMSBE-RCAV5-1D-n2p5-colThresh-2em2-theta2-6400lam-spontEmis-noExpSBE-restarted-out47-OUTPUTBACK-pulseSpectrum1.mat'};
% key(3)={'tMSBE-RCAV4-1D-n2p5-colThresh-2em2-theta2-6400lam-spontEmis-wExpSBE-restarted-out45-OUTPUT-pulseSpectrum1.mat'};
% key(4)={'tMSBE-RCAV4-1D-n2p5-colThresh-2em2-theta2-6400lam-spontEmis-wExpSBE-restarted-out45-OUTPUTBACK-pulseSpectrum1.mat'};
% 
% savedData(1).legend={'CW 1st order'};
% savedData(2).legend={'CCW 1st order'};
% savedData(3).legend={'CW 3rd order'};
% savedData(4).legend={'CCW 3rd order'};

key(1)={'tMSBE-RCAV10-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus10-out29-OUTPUT-pulseSpectrum1.mat'};
key(2)={'tMSBE-RCAV10-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus10-out29-OUTPUTBACK-pulseSpectrum1.mat'};
savedData(1).legend={'CW'};
savedData(2).legend={'CCW'};

test_folder='test';
date='121820';
saveKey_local='Fall2020-Summer2021/RingCAV/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey=[saveKey_fold,'/RCAV4-RCAV5-'];

for l=1:length(key)
    savedData(l).key = strcat(saveKey_fold,'/',key(l));
end

compareModelockedSpectrums(savedData,saveKey)