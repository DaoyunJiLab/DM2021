function DataManager_Save_Callback
hf = gcbf; fname = get(hf, 'Name');
mmind = strfind(fname, '__'); extfile = fname(mmind+2:numel(fname)); 
[pp, nn, ee] = fileparts(extfile);

if (strcmp(ee, '.behavdb'))
   behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
   
   save(extfile, 'behav', 'bhdata');
   %save(extfile, 'behav', 'bhdata', '-mat', '-v7.3');
elseif (strcmp(ee, '.eegdb'))
   eeg = getappdata(hf, 'eeg'); eegdata = getappdata(hf, 'eegdata');
   save(extfile, 'eeg', 'eegdata');
   %save(extfile, 'eeg', 'eegdata', '-mat', '-v7.3');
else %(strcmp(ee, '.spikedb'))
   pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
   s = whos('data'); Gb1 = s.bytes/(1024^3);
   s = whos('pinfo'); Gb2 = s.bytes/(1024^3);
   if Gb1 + Gb2 < 2
      save(extfile, 'pinfo', 'data'); %disp('Small one');
   else
      save(extfile, 'pinfo', 'data', '-mat', '-v7.3'); %disp('Big one');
   end
end
