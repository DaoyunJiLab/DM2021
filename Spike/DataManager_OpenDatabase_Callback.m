function DataManager_OpenDatabase_Callback
hmain = gcbf;
cpathname = fullfile(cd, '*.*db'); 
[fname, pname] = uigetfile(cpathname, 'Select a database file to open:');
if (fname ~= 0)
   filename = strcat(pname, fname);
   [~, ~, ee] = fileparts(filename);
   switch ee
       case '.spikedb'
             disp('-----> reading spike database');
             S = load(filename, '-mat'); %load the file
             pinfo = S.pinfo; data = S.data; %load the plot structure
             if (~isempty(pinfo.general.parmfile))
                 setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); S = [];
                 DataManager_PlotSpikeDatabase(hmain, pinfo, data);
                 set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', filename));
             else
                 disp('-------------> no cells in the data base file!');
             end
       case '.eegdb'
             disp('-----> read EEG database');
             S = load(filename, '-mat'); %load the file
             eeg = S.eeg; eegdata = S.eegdata; %load the plot structure
             if (~isempty(eeg.general.eegfile))
                setappdata(hmain, 'eeg', eeg); setappdata(hmain, 'eegdata', eegdata); S = [];
                DataManager_PlotSpikeDatabase(hmain, eeg, eegdata, 'EEG File', 'eegID', '.eegdb');
                set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', filename));
             else
                disp('-------------> no eeg files in the data base!');
             end
       case '.behavdb'
             disp('-----> read behavioral database');
             S = load(filename, '-mat'); %load the file
             behav = S.behav; bhdata = S.bhdata; %load the plot structure
             if (~isempty(behav.general.datedir))
                setappdata(hmain, 'behav', behav); setappdata(hmain, 'bhdata', bhdata); S = [];
                DataManager_PlotSpikeDatabase(hmain, behav, bhdata, 'Session', 'sessID', '.behavdb');
                set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', filename));
             else
                disp('-------------> no behavioral sesesions in the data base file!');
             end
       case '.crrdb'
           disp('-----> read cross correlation database');
             S = load(filename, '-mat'); %load the file
             pinfo = S.pinfo; data = S.data; %load the plot structure
             if (~isempty(pinfo.general.clname))
                setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); S = [];
                DataManager_PlotSpikeDatabase(hmain, pinfo, data, 'CellPair', 'clname', '.crrdb'); 
                set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', filename));
             else
                disp('-------------> no cell pairs in the data base file!');
             end
       case '.seqdb'
           disp('-----> read sequence database');
             S = load(filename, '-mat'); %load the file
             pinfo = S.pinfo; data = S.data; %load the plot structure
             if (~isempty(pinfo.general.finaldir))
                setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); S = [];
                if isfield(pinfo.parm, 'seqtype') && contains(pinfo.parm.seqtype{1}, 'evtitemized')
                   DataManager_PlotSpikeDatabase(hmain, pinfo, data, 'Event', 'eventname', '.seqdb');
                elseif isfield(pinfo.general, 'tmpID')
                   DataManager_PlotSpikeDatabase(hmain, pinfo, data, 'Template', 'tmpID', '.seqdb');
                end
                set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', filename));
             else
                disp('-------------> no cell pairs in the data base file!');
             end
   end
end
