function DataManager_OpenEEGDatabase_Callback
hmain = gcbf;
cpathname = fullfile(cd, '*.eegdb'); 
[fname, pname] = uigetfile(cpathname, 'Select an eeg database file to open:');
if (fname ~= 0)
filename = strcat(pname, fname);
disp('-----> read database');
S = load(filename, '-mat'); %load the file
eeg = S.eeg; eegdata = S.eegdata; %load the plot structure
if (~isempty(eeg.general.eegfile))
   setappdata(hmain, 'eeg', eeg); setappdata(hmain, 'eegdata', eegdata); S = [];
   DataManager_PlotEEGDatabase(hmain, eeg, eegdata);
   set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', filename));
else
    disp('-------------> no eeg files in the data base!');
end
end
