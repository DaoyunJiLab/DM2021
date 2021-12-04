function [eeg, eegdata] = DataManager_EEGComputeInit_Callback(eeg, eegdata, eegind, vv)
%compute initial variables of a eegdatabase
if (~isempty(eegind))
   disp('------>Compute intitial variables of the eeg database');
   
   %disp('-------> compute spectral properties');
   %[eeg, eegdata] = DataManager_FindSpectralProp(eeg, eegdata, eegind, vv);
  
   %disp('-------> compute ripple properties');
   %[eeg, eegdata] = DataManager_FindeegRippleProp(eeg, eegdata, eegind, vv);
   
   %disp('-------> compute theta properties');
   %[eeg, eegdata] = DataManager_FindeegThetaProp(eeg, eegdata, eegind, vv);

end 


