function DataManager_Open_Callback
hmain = gcbf;
cpathname = fullfile(cd, '*.db'); 
[fname, pname] = uigetfile(cpathname, 'Select a database file to open:');
filename = strcat(pname, fname);
disp('-----> read database');
S = load(filename, '-mat'); %load the file
pinfo = S.pinfo; %load the plot structure
setappdata(hmain, 'pinfo', pinfo); S = [];
if (pinfo.databasetype == 'spike')
    DataManager_PlotSpikeDatabase(hmain, pinfo);
elseif (pinfo.databasetype == 'eeg')
    DataManager_PlotEEGDatabase(hmain, pinfo);
end