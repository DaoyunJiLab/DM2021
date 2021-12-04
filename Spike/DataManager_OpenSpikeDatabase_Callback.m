function DataManager_OpenSpikeDatabase_Callback
hmain = gcbf;
cpathname = fullfile(cd, '*.*db'); 
[fname, pname] = uigetfile(cpathname, 'Select a spike database file to open:');
if (fname ~= 0)
filename = strcat(pname, fname);
disp('-----> read database');
S = load(filename, '-mat'); %load the file
pinfo = S.pinfo; data = S.data; %load the plot structure
if (~isempty(pinfo.general.parmfile))
   setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); S = [];
   DataManager_PlotSpikeDatabase(hmain, pinfo, data);
   set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', filename));
else
    disp('-------------> no cells in the data base file!');
end
end
