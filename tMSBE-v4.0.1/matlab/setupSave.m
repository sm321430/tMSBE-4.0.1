%% Setup load and save parameters
user_entry = input('saveKey Updated??? (y to continue): ', 's');
if user_entry~='y'
    gafuggle
end

location='OUTPUT-OUTPUTBACK'; %Field location for uploading and saving
date='050321';
test='tMSBE-dualProp-sepChips1';
%test='tMSBE-ECAV-OGSG-Focus1-OCcomp-ECAV28-ECAV37-38';
%test='tMSBE-dualProp-fullCavity2-n3p0-comp';
%test='tMSBE-DualPropCavity5-n3p0-BRClength-10lambda';
%test='tMSBE-v4.0.1-DualPropCavity5-n3p0-BRClength-10lambda';
%test='tMSBE-v4.0.1-DualPropCavity8-n3p0-BRClength-1000lambda';
%test='tMSBE-v4.0.1-DualPropCavity3-cavityModesBalance';
%test='tMSBE-dualProp-fullCavity4-twoToOneCavity-comp';
%test='tMSBE-ECAV-biDirControl-chipPlacement-ECAV31-ECAV41-42';
%test='tMSBE-ECAV-biDirControl-900lambda-ECAV34-35-ECAV39';
%test='tMSBE-GSGO-highDensFocusTest-n2p5-ECAV26-28';
%test='tMSBE-ECAV-biDirControl-ECAV20-ECAV25';
%test='tMSBE-GSGO-highDensFocusTest-n2p5-ECAV13-ECAV21';
%test='tMSBE-RCAV-n3p5-lengthsComp-RCAV102-RCAV106-109-RCAV119-120-reordered';
%test='tMSBE-OGSG-unequalPumpComp-ECAV08-ECAV16-18';
%test='tMSBE-PentaCAV-ECAV04-ECAV09';
%test='tMSBE-OGSG-ECAV01-ECAV08-ECAV13-15';
%test='tMSBE-RCAV-2Dcomp-RCAV99-RCAV122';
%test='tMSBE-RCAV10-RCAV20-21-RCAV25-28';
%test='tMSBE-lengthsComp-RCAV10-RCAV20-21-RCAV25-RCAV27';
%test='tMSBE-RCAV-PhaseSpace-20_21_25-32_38-43_51-83';
%test='tMSBE-RCAV-1D-n2p7-RCAV45-50-RCAV90';
%test='tMSBE-RCAV-1D-n2p9-RCAV84-89-RCAV91';
%test='tMSBE-RCAV-1D-n3p1-RCAV96-98-RCAV100';
%test='tMSBE-RCAV45-ECAV01-02-ECAV04';
%test='tMSBE-1D-densityComp-Focus-10-3p1-3p5-3p9';
%test='tMSBE-extendedVCAV-ECAV03-ECAV05-ECAV06-ECAV10-12';
%test='tMSBE-RCAV-longCav1-RCAV111-114';
%test='tMSBE-RCAV-longCav2-33p5mm-densityComp-RCAV113-RCAV115-117';
%test='tMSBE-RCAV96-98-1D-n3p1';
%test='tMSBE-RCAV37-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5';
%test='tMSBE-RCAV95-2D-n2p5-colThresh-theta2-1700-1500-3200lam-spontEmis-wExpSBE-focus5-seedPulse-2e5';
test_folder='test';
saveKey_local='Fall2020-Summer2021/DualProp/';
saveKey_fold=['../../Research/Notes/',saveKey_local,date,'/',test_folder];
saveKey = [saveKey_fold,'/',test,'-out','-',location,'-'];

if exist(saveKey_fold,'dir')~=7
    mkdir(saveKey_fold);
end