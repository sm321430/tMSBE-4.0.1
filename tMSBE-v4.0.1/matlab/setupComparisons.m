IO_focusTest=0; %Turn on focus based legend extraction from file names

%simNums=[20,21,25:32,38:43,51:83];
%simNums=[45:46,90,47:50];
%simNums=[84,85,91,86:89];

% simNums=[21,25,20,10,27];
% simType={'R','R','R','R','R'};
% CAVsim(1).legend={'0.824mm, CW'};
% CAVsim(2).legend={'1.318mm, CW'};
% CAVsim(3).legend={'1.545mm, CW'};
% CAVsim(4).legend={'1.648mm, CW'};
% CAVsim(5).legend={'1.751mm, CW'};

%simNums=[45,1,2,4];
%simType={'R','M','M','M'};
%CAVsim(1).legend={'OGS ring cavity, CW'};
%CAVsim(2).legend={'OGSG ring cavity, CW'};
%CAVsim(3).legend={'OSGS ring cavity, CW'};
%CAVsim(4).legend={'OGSGS ring cavity, CW'};
% 
% simNums=[6,3,5,11,10,12];
% simType={'M','M','M','M','M', 'M'};
% CAVsim(1).legend={'V-cavity'};
% CAVsim(2).legend={'W-cavity'};
% CAVsim(3).legend={'VW-cavity'};
% CAVsim(4).legend={'Balanced W-cavity'};
% CAVsim(5).legend={'Balanced VW-cavity'};
% CAVsim(6).legend={'Single SESAM W-cavity'};

% simNums=111:114;
% simType={'R','R','R','R'};
% CAVsim(1).legend={'6.7mm, CW'};
% CAVsim(2).legend={'13.4mm, CW'};
% CAVsim(3).legend={'33.5mm, CW'};
% CAVsim(4).legend={'67mm, CW'};

% simNums=[113,115:117];
% simType={'R','R','R','R'};
% CAVsim(1).legend={'2.7e16m^{-2}, CW'};
% CAVsim(2).legend={'3.1e16m^{-2}, CW'};
% CAVsim(3).legend={'3.5e16m^{-2}, CW'};
% CAVsim(4).legend={'3.9e16m^{-2}, CW'};

% simNums=[96,102,104];
% simType={'R','R','R','R','R'};
% CAVsim(1).legend={'2.5e16m^{-2}, CW'};
% CAVsim(2).legend={'2.7e16m^{-2}, CW'};
% CAVsim(1).legend={'2.9e16m^{-2}, CW'};
% CAVsim(1).legend={'3.1e16m^{-2}, CW'};
% CAVsim(3).legend={'3.3e16m^{-2}, CW'};
% CAVsim(2).legend={'3.5e16m^{-2}, CW'};
 %CAVsim(3).legend={'3.7e16m^{-2}, CW'};
% CAVsim(3).legend={'3.9e16m^{-2}, CW'};
% CAVsim(5).legend={'4.1e16m^{-2}, CW'};
 
% simNums=[04,09];
% simType={'M','M'};
% CAVsim(1).legend={'2.7e16m^{-2}, CW'};
% CAVsim(2).legend={'3.1e16m^{-2}, CW'};

% simNums=[01,08,13:15];
% simType={'M','M','M','M','M'};
% CAVsim(1).legend={'1.95e16m^{-2}, CW'};
% CAVsim(2).legend={'2.2e16m^{-2}, CW'};
% CAVsim(3).legend={'2.5e16m^{-2}, CW'};
% CAVsim(4).legend={'2.9e16m^{-2}, CW'};
% CAVsim(5).legend={'3.3e16m^{-2}, CW'};

% simNums=[08,16:18];
% simType={'M','M','M','M','M'};
% CAVsim(1).legend={'2.2e16m^{-2}, CW'};
% CAVsim(2).legend={'2.5e16m^{-2}, CW'};
% CAVsim(3).legend={'2.9e16m^{-2}, CW'};
% CAVsim(4).legend={'3.3e16m^{-2}, CW'};
%CAVsim(5).legend={'3.7e16m^{-2}, CW'};

% simNums=[120,107,102,106,119];
% simType={'R','R','R','R','R','R','R','R','R',};
% CAVsim(1).legend={'1.854mm, F=10, CW'};
% CAVsim(2).legend={'1.803mm, F=10, CW'};
% CAVsim(3).legend={'1.751mm, F=10, CW'};
% CAVsim(4).legend={'1.700mm, F=10, CW'};
% CAVsim(5).legend={'1.545mm, F=10, CW'};

% simNums=[20,25];
% simType={'M','M','M'};
% CAVsim(1).legend={'2.2e16m^{-2}, 2.2e16m^{-2}, CW'};
% CAVsim(2).legend={'2.2e16m^{-2}, 2.5e16m^{-2}, CW'};
% CAVsim(2).legend={'1.9e16m^{-2}, 2.2e16m^{-2}, CW'};

% simNums=[13,21];
% simType={'M','M','M','M'};
% CAVsim(1).legend={'F=10, CW'};
% CAVsim(2).legend={'F=4, CW'};
% CAVsim(3).legend={'F=2, CW'};
% CAVsim(4).legend={'F=1, CW'};

% simNums=26:28;
% simType={'M','M','M','M'};
% %CAVsim(1).legend={'F=10, CW'};
% CAVsim(1).legend={'F=4, CW'};
% CAVsim(2).legend={'F=2, CW'};
% CAVsim(3).legend={'F=1, CW'};

% simNums=29:31;
% simType={'M','M','M'};
% CAVsim(1).legend={'309\mu m, CW'};
% CAVsim(2).legend={'618\mu m, CW'};
% CAVsim(3).legend={'927\mu m, CW'};

% simNums=[34,35,39,40];
% simType={'M','M','M','M','M'};
% CAVsim(1).legend={'1.9e16m^{-2}, 1.9e16m^{-2}, CW'};
% %CAVsim(2).legend={'1.9e16m^{-2}, 2.2e16m^{-2}, CW'};
% %CAVsim(1).legend={'1.9e16m^{-2}, 2.35e16m^{-2}, CW'};
% %CAVsim(4).legend={'1.75e16m^{-2}, 2.2e16m^{-2}, CW'};
% %CAVsim(2).legend={'1.75e16m^{-2}, 2.35e16m^{-2}, CW'};
% %CAVsim(3).legend={'1.6e16m^{-2}, 2.35e16m^{-2}, CW'};
% CAVsim(4).legend={'1.9e16m^{-2}, 1.8e16m^{-2}, CW'};
% 
% simNums=[31,41,42];
% simType={'M','M','M','M','M'};
% CAVsim(1).legend={'1.9e16m^{-2}, 1.9e16m^{-2}, CW'};
% CAVsim(2).legend={'2.1e16m^{-2}, 1.7e16m^{-2}, CW'};
% CAVsim(3).legend={'2.3e16m^{-2}, 1.5e16m^{-2}, CW'};

% simNums=[28,37,38];
% simType={'M','M','M','M','M'};
% CAVsim(1).legend={'OC=1.0%, CW'};
% CAVsim(2).legend={'OC=0.75%, CW'};
% CAVsim(3).legend={'OC=0.5%, CW'};

% 
%  simNums=[5,5];
%  simType={'D','D'};
%  CAVsim(1).legend={'ordinary'};
%  CAVsim(2).legend={'extraordinary'};
%  CAVsim(1).location='OUTPUT';
%  CAVsim(2).location='OUTPUTBACK';

CAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-dualProp-data-2021/tMSBE-v4.2-KerrLensSML1_untested_DualProp_sepChips/run/out__'};
CAVsim(2).outKey=CAVsim(1).outKey;
%CAVsim(3).outKey=CAVsim(1).outKey;
%CAVsim(1).location='CAVL';
CAVsim(1).location='OUTPUT';
CAVsim(2).location='OUTPUTBACK';
%CAVsim(1).legend=CAVsim(1).location;
CAVsim(1).legend={'ordinary'};
CAVsim(2).legend={'extraordinary'};
%CAVsim(1).legend=CAVsim(1).location;
%CAVsim(2).legend=CAVsim(2).location;
simNums=[];

for j=1:length(simNums)
    CAVsim(j).location='OUTPUT';
    %populate list of simulation outKeys
    if simType{j}=='R'
        mainDir ='/Volumes/SAMbackup/tMSBE-RCAV-data-2021/';
        full=strcat(mainDir,'tMSBE-v3.8-RCAV',num2str(simNums(j)),'*');
    elseif simType{j}=='M'
        if simNums(j)<10
            mainDir ='/Volumes/SAMbackup/tMSBE-ECAV-data-2021/';
            full=strcat(mainDir,'tMSBE-v3.8-ECAV-0',num2str(simNums(j)),'*');
        else
            mainDir ='/Volumes/SAMbackup/tMSBE-ECAV-data-2021/';
            full=strcat(mainDir,'tMSBE-v3.8-ECAV-',num2str(simNums(j)),'*');    
        end
    elseif simType{j}=='D'
        mainDir ='/Volumes/SAMbackup/tMSBE-dualProp-data-2021/';
        full=strcat(mainDir,'tMSBE-v4.0.1-DualPropCavity',num2str(simNums(j)),'*');
    else
       disp('Error: simType not recognized. BREAKING');
       asd
    end
    fstruct = dir(full);
    if length(fstruct)>1
       disp(strcat('Multiple RCAV',num2str(simNums(j)),' directories.'));
       simName=fstruct(1).name;
    else
        simName=fstruct.name;
    end
    CAVsim(j).outKey=strcat(mainDir,simName,'/run/out__');
    
    if IO_focusTest==1
        %populate phaseSpace coordinates
        dashes=strfind(simName,strcat('-'));
        focus_startInd=strfind(simName,'focus')+5;
        focus_stopInd=dashes(find(dashes>focus_startInd,1));
        if isempty(focus_stopInd)
            focus_stopInd=length(simName);
        end
        CAVsim(j).phaseSpacePos(1)=str2double(simName(focus_startInd:focus_stopInd));

        OC_startInd=strfind(simName,'theta2-')+7;
        OC_stopInd=dashes(find(dashes>OC_startInd,1))-1;
        OC_length=str2double(simName(OC_startInd:OC_stopInd));
        SESAM_stopInd=dashes(find(dashes>OC_stopInd+1,1))-1;
        SESAM_length=str2double(simName(OC_stopInd+2:SESAM_stopInd));
        CAVsim(j).phaseSpacePos(2)=OC_length/6400;%/SESAM_length;
        CAVsim(j).legend={strcat('F=',simName(focus_startInd:focus_stopInd),', CW')};
    end
    j
end

CAVsim(1).cavType='RCAV'; %Includes both directions automatically
%CAVsim(1).phaseSpaceLabel={'Focus [-]','OC Arm Length [L_{tot}]'};
num_sims=length(CAVsim);

% asd
% CAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV28-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% %CAVsim(1).phaseSpacePos=[10,1.5656]; %position in phase space (focus, length in mm)
% CAVsim(1).phaseSpacePos=[10,0.2625];
% 
% CAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV30-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus8/run/out__'};
% %CAVsim(2).phaseSpacePos=[8,1.5656];
% CAVsim(2).phaseSpacePos=[8,0.2625];
% 
% CAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV29-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus5/run/out__'};
% %CAVsim(3).phaseSpacePos=[5,1.5656];
% CAVsim(3).phaseSpacePos=[5,0.2625];
% 
% CAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV31-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus4/run/out__'};
% %CAVsim(4).phaseSpacePos=[4,1.5656];
% CAVsim(4).phaseSpacePos=[4,0.2625];
% 
% CAVsim(5).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV32-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus2/run/out__'};
% %CAVsim(5).phaseSpacePos=[2,1.5656];
% CAVsim(5).phaseSpacePos=[2,0.525];
% 
% CAVsim(6).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV27-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% %CAVsim(6).phaseSpacePos=[10,1.545];
% CAVsim(6).phaseSpacePos=[10,0.2656];
% 
% CAVsim(7).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV40-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus8/run/out__'};
% %CAVsim(7).phaseSpacePos=[8,1.545];
% CAVsim(7).phaseSpacePos=[8,0.2656];
% 
% CAVsim(8).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV41-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5/run/out__'};
% %CAVsim(8).phaseSpacePos=[5,1.545];
% CAVsim(8).phaseSpacePos=[5,0.2656];
% 
% CAVsim(9).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV42-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus4/run/out__'};
% %CAVsim(9).phaseSpacePos=[4,1.545];
% CAVsim(9).phaseSpacePos=[4,0.2656];
% 
% CAVsim(10).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV43-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus2/run/out__'};
% %CAVsim(10).phaseSpacePos=[2,1.545];
% CAVsim(10).phaseSpacePos=[2,0.5312];
% 
% CAVsim(11).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV39-1D-n2p5-noThresh-theta2-1600-1600-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% %CAVsim(11).phaseSpacePos=[10,1.648];
% CAVsim(11).phaseSpacePos=[10,0.125];
% 
% CAVsim(12).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV20-1D-n2p5-noThresh-theta2-1500-1700-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% %CAVsim(12).phaseSpacePos=[10,1.751];
% CAVsim(12).phaseSpacePos=[10,0.2344];
% 
% CAVsim(13).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV25-1D-n2p5-noThresh-theta2-1280-1920-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% %CAVsim(13).phaseSpacePos=[10,1.9776];
% CAVsim(13).phaseSpacePos=[10,0.2];
% 
% CAVsim(14).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV26-1D-n2p5-noThresh-theta2-1920-1280-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% %CAVsim(14).phaseSpacePos=[10,1.3184];
% CAVsim(14).phaseSpacePos=[10,0.3];
% 
% CAVsim(15).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV21-1D-n2p5-noThresh-theta2-800-2400-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% %CAVsim(15).phaseSpacePos=[10,2.472];
% CAVsim(15).phaseSpacePos=[10,0.25];
% 
% CAVsim(16).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV38-1D-n2p5-noThresh-theta2-1650-1550-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% %CAVsim(16).phaseSpacePos=[10,1.5965];
% CAVsim(16).phaseSpacePos=[10,0.2578];


% % CAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-CAV-data-2020/tMSBE-v3.7-RCAV4-1d-n2p5-colThresh-2em2-theta2-6400lam-spontEmis-wExpSBE/run/out__'};
% % CAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-CAV-data-2020/tMSBE-v3.7-RCAV4-1d-n2p5-colThresh-2em2-theta2-6400lam-spontEmis-wExpSBE/run/out__'};
% % CAVsim(1).legend={'Clockwise'};
% % CAVsim(2).legend={'Ctr. Clockwise'};
% % CAVsim(1).location='OUTPUT';
% % CAVsim(2).location='OUTPUTBACK';

% CAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.7-RCAV5-1D-n2p5-colThresh-2em2-theta2-6400lam-spontEmis-noExpSBE-restarted/run/out__'};
% CAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.7-RCAV5-1D-n2p5-colThresh-2em2-theta2-6400lam-spontEmis-noExpSBE-restarted/run/out__'};
% CAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.7-RCAV4-1D-n2p5-colThresh-2em2-theta2-6400lam-spontEmis-wExpSBE-restarted/run/out__'};
% CAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.7-RCAV4-1D-n2p5-colThresh-2em2-theta2-6400lam-spontEmis-wExpSBE-restarted/run/out__'};
% CAVsim(1).legend={'CW 1st order'};
% CAVsim(2).legend={'CCW 1st order'};
% CAVsim(3).legend={'CW 3rd order'};
% CAVsim(4).legend={'CCW 3rd order'};
% CAVsim(1).location='OUTPUT';
% CAVsim(2).location='OUTPUTBACK';
% CAVsim(3).location='OUTPUT';
% CAVsim(4).location='OUTPUTBACK';

% CAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV10-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV10-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV12-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus8/run/out__'};
% CAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV12-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus8/run/out__'};
% CAVsim(5).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV11-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus5/run/out__'};
% CAVsim(6).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV11-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus5/run/out__'};
% CAVsim(7).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV23-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus4/run/out__'};
% CAVsim(8).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV23-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus4/run/out__'};
% CAVsim(9).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV24-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus2/run/out__'};
% CAVsim(10).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV24-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus2/run/out__'};
% CAVsim(11).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV13-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus1/run/out__'};
% CAVsim(12).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV13-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus1/run/out__'};
% 
% 
% CAVsim(1).legend={'F=10, CW'};
% CAVsim(2).legend={'F=10, CCW'};
% CAVsim(3).legend={'F=8, CW'};
% CAVsim(4).legend={'F=8, CCW'};
% CAVsim(5).legend={'F=5, CW'};
% CAVsim(6).legend={'F=5, CCW'};
% CAVsim(7).legend={'F=4, CW'};
% CAVsim(8).legend={'F=4, CCW'};
% CAVsim(9).legend={'F=2, CW'};
% CAVsim(10).legend={'F=2, CCW'};
% CAVsim(11).legend={'F=1, CW'};
% CAVsim(12).legend={'F=1, CCW'};
% 
% CAVsim(1).location='OUTPUT';
% CAVsim(2).location='OUTPUTBACK';
% CAVsim(3).location='OUTPUT';
% CAVsim(4).location='OUTPUTBACK';
% CAVsim(5).location='OUTPUT';
% CAVsim(6).location='OUTPUTBACK';
% CAVsim(7).location='OUTPUT';
% CAVsim(8).location='OUTPUTBACK';
% CAVsim(9).location='OUTPUT';
% CAVsim(10).location='OUTPUTBACK';
% CAVsim(11).location='OUTPUT';
% CAVsim(12).location='OUTPUTBACK';

% CAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV10-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV10-1D-n2p5-noThresh-theta2-6400lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(1).legend={'CW, control'};
% CAVsim(2).legend={'CCW, control'};
% CAVsim(1).location='OUTPUT';
% CAVsim(2).location='OUTPUTBACK';
% 
% CAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV17-1D-2np5-noThresh-theta2-6400lam-noSpont-wExpSBE-focus10/run/out__'};
% CAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV17-1D-2np5-noThresh-theta2-6400lam-noSpont-wExpSBE-focus10/run/out__'};
% CAVsim(3).legend={'CW, noSpont'};
% CAVsim(4).legend={'CCW, noSpont'};
% CAVsim(3).location='OUTPUT';
% CAVsim(4).location='OUTPUTBACK';
% 
% 
% CAVsim(5).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV19-1D-n2p5-colThresh-theta2-6400lamFixed-SpontEmis-wExpSbe-focus10/run/out__'};
% CAVsim(6).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV19-1D-n2p5-colThresh-theta2-6400lamFixed-SpontEmis-wExpSbe-focus10/run/out__'};
% CAVsim(5).legend={'CW, colThresh'};
% CAVsim(6).legend={'CCW, colThresh'};
% CAVsim(5).location='OUTPUT';
% CAVsim(6).location='OUTPUTBACK';

% CAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV16-debugger/run/out__'};
% CAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV16-debugger/run/out__'};
% CAVsim(1).legend={'CW'};
% CAVsim(2).legend={'CCW'};
% CAVsim(1).location='OUTPUT';
% CAVsim(2).location='OUTPUTBACK';

% CAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV39-1D-n2p5-noThresh-theta2-1600-1600-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV39-1D-n2p5-noThresh-theta2-1600-1600-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV20-1D-n2p5-noThresh-theta2-1500-1700-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV20-1D-n2p5-noThresh-theta2-1500-1700-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(5).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV27-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(6).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV27-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(7).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV25-1D-n2p5-noThresh-theta2-1280-1920-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(8).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV25-1D-n2p5-noThresh-theta2-1280-1920-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(9).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV26-1D-n2p5-noThresh-theta2-1920-1280-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(10).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV26-1D-n2p5-noThresh-theta2-1920-1280-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(11).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV21-1D-n2p5-noThresh-theta2-800-2400-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(12).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV21-1D-n2p5-noThresh-theta2-800-2400-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(13).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV28-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(14).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV28-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(15).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV38-1D-n2p5-noThresh-theta2-1650-1550-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(16).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV38-1D-n2p5-noThresh-theta2-1650-1550-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(1).legend={'CW, 1:1'};
% CAVsim(2).legend={'CCW, 1:1'};
% CAVsim(3).legend={'CW, 1:1.13'};
% CAVsim(4).legend={'CCW, 1:1.13'};
% CAVsim(5).legend={'CW, 1.13:1'};
% CAVsim(6).legend={'CCW, 1.13:1'};
% CAVsim(7).legend={'CW, 2:3'};
% CAVsim(8).legend={'CCW, 2:3'};
% CAVsim(9).legend={'CW, 3:2'};
% CAVsim(10).legend={'CCW, 3:2'};
% CAVsim(11).legend={'CW, 1:3'};
% CAVsim(12).legend={'CCW, 1:3'};
% CAVsim(13).legend={'CW, 1.11:1'};
% CAVsim(14).legend={'CCW, 1.11:1'};
% CAVsim(15).legend={'CW, 1.06:1'};
% CAVsim(16).legend={'CCW, 1.06:1'};

% CAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV27-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV27-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV40-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus8/run/out__'};
% CAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV40-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus8/run/out__'};
% CAVsim(5).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV41-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5/run/out__'};
% CAVsim(6).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV41-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5/run/out__'};
% CAVsim(7).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV42-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus4/run/out__'};
% CAVsim(8).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV42-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus4/run/out__'};
% CAVsim(9).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV43-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus2/run/out__'};
% CAVsim(10).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV43-1D-n2p5-noThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus2/run/out__'};
% 
% CAVsim(1).legend={'CW, F=10'};
% CAVsim(2).legend={'CCW, F=10'};
% CAVsim(3).legend={'CW, F=8'};
% CAVsim(4).legend={'CCW, F=8'};
% CAVsim(5).legend={'CW, F=5'};
% CAVsim(6).legend={'CCW, F=5'};
% CAVsim(7).legend={'CW, F=4'};
% CAVsim(8).legend={'CCW, F=4'};
% CAVsim(9).legend={'CW, F=2'};
% CAVsim(10).legend={'CCW, F=2'};
% 
% for j=1:length(CAVsim)/2
%     CAVsim(2*j-1).location='OUTPUT';
%     CAVsim(2*j).location='OUTPUTBACK';
% end
% 
% CAVsim(1).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV28-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(2).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV28-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus10/run/out__'};
% CAVsim(3).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV30-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus8/run/out__'};
% CAVsim(4).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV30-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus8/run/out__'};
% CAVsim(5).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV29-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus5/run/out__'};
% CAVsim(6).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV29-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus5/run/out__'};
% CAVsim(7).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV31-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus4/run/out__'};
% CAVsim(8).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV31-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus4/run/out__'};
% % CAVsim(9).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV32-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus2/run/out__'};
% % CAVsim(10).outKey={'/Volumes/SAMbackup/tMSBE-RCAV-data-2021/tMSBE-v3.8-RCAV32-1D-n2p5-noThresh-theta2-1680-1520-3200lam-spontEmis-wExpSBE-focus2/run/out__'};
% CAVsim(1).legend={'CW, F=10'};
% CAVsim(2).legend={'CCW, F=10'};
% CAVsim(3).legend={'CW, F=8'};
% CAVsim(4).legend={'CCW, F=8'};
% CAVsim(5).legend={'CW, F=5'};
% CAVsim(6).legend={'CCW, F=5'};
% CAVsim(7).legend={'CW, F=4'};
% CAVsim(8).legend={'CCW, F=4'};
% % CAVsim(9).legend={'CW, F=2'};
% % CAVsim(10).legend={'CCW, F=2'};
% for j=1:length(CAVsim)/2
%     CAVsim(2*j-1).location='OUTPUT';
%     CAVsim(2*j).location='OUTPUTBACK';
% end
% 
% num_sims=length(CAVsim);
% 
% for j = 1 : num_sims
%    if(mod(j,2)==0)
%        CAVsim(j).location='OUTPUTBACK';
%    else
%        CAVsim(j).location='OUTPUT';
%    end
% end