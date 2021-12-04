function [behav, bhdata] = DataManager_BehavComputeInit_Callback(behav, bhdata, dateind, vv)
%compute initial variables of a bhdatabase
if (~isempty(dateind))
   disp('----->Compute intitial variables of the behaviroal database');
   
   disp('--------> compute position data properties');
   [behav, bhdata] = DataManager_FindBehavPositionProp(behav, bhdata, dateind, vv);
  
   %disp('-----> compute event properties');
   %[behav, bhdata] = DataManager_FindBehavEventProp(behav, bhdata, dateind, vv);
   
%get behavioral/fields information -- not done yest!
%%%%%%%%%disp('-----> get behavioral properties for all spikes');
%behav = bhdataManager_FindBehaviorInfo(behav);
        
%%get spatial information -- not done yet!
%behav = Manage_FindSpatialInfo(behav);
  
%get interneuron/pyramidal cell classification -- NOT DONE yet!, don't do it at the beginning: it requires user interference 
% disp('-----> classify neuron type');
% behav = Manage_FindNeuronType(behav);

%get phase information - NOT done yet!!! 
%disp('-----> get phase properties for all spikes');
%behav = Manage_FindPhaseInfo(behav); There Maybe a new program doing this - chech TrajSwitch - phase programs
end 


