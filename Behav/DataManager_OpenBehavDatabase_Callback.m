function DataManager_OpenBehavDatabase_Callback
hmain = gcbf;
cpathname = fullfile(cd, '*.behavdb'); 
[fname, pname] = uigetfile(cpathname, 'Select a behavior database file to open:');
if (fname ~= 0)
filename = strcat(pname, fname);
disp('-----> read database');
S = load(filename, '-mat'); %load the file
behav = S.behav; bhdata = S.bhdata; %load the plot structure
if (~isempty(behav.general.datedir))
   setappdata(hmain, 'behav', behav); setappdata(hmain, 'bhdata', bhdata); S = [];
   DataManager_PlotBehavDatabase(hmain, behav, bhdata);
   set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', filename));
else
    disp('-------------> no cells in the data base file!');
end
end
