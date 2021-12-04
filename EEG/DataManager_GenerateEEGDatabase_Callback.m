function DataManager_GenerateEEGDatabase_Callback
%%orgainze EEG data from an animal's all recording dates
%%calculate features associate for each spike and set to corresponding fields
%%standard features are generated here
%%resulted structure is saved as a -mat type file with extension .spikedb
%%new features can be added later by openning the file

%get the animal name (animal root directory)
%animalpath = uigetdir(cd); %get just one animal dir
datepath = []; ok = 1;  curdir = cd; %%%current directory
aaa = questdlg('Import a final-directory list?');
if (strcmp(aaa, 'Yes'))
   [fname, pname] = uigetfile(fullfile(cd, '*.lst'), 'Choose a list:');
   if fname == 0
       ok = 0; 
   else
       datepath = textread(fullfile(pname, fname), '%s', 'delimiter', ' ', 'commentstyle', 'matlab');
       %%%%check datepath
       iii = []; nw = 0;
       for (i = 1:numel(datepath))
         if (exist(datepath{i}, 'dir') == 7) && (~isempty(strfind(datepath{i}, '\final')))
            iii = union(iii, i); 
         elseif ~isempty(datepath{i})
            nw = nw + 1;
            disp(['---> warning: not a \final dir or not exist: ', datepath{i}]); 
         end
       end
      if (nw > 0)
        aaa = questdlg(['number of dirs that will not be processed: ', num2str(nw), '. Continue?']);
        if strcmp(aaa, 'Yes')
           datepath = datepath(iii);
        else
           datepath = [];
        end
      else
        datepath = datepath(iii);
      end
   end
elseif (strcmp(aaa, 'No'))
   [datepath, ok] = DirDlg(curdir);
   cd(curdir);
end
datepath = unique(datepath);
if ok
    copyover = questdlg('Copy over raw EEG data?'); eegdir = [];
    if (strcmp(copyover, 'Cancel')) 
        ok = 0; 
    elseif strcmp(copyover, 'Yes')
        eegdir = uigetdir(curdir, 'Choose a directory for copying eeg data:');
        if (eegdir ==0) ok = 0; end
    end
end
if ok && (numel(datepath)>0)
    disp('Create new EEG database');
    %disp(['---> ', animalpath]);
    disp(['---> number of date final directories: ', num2str(numel(datepath))]);
    [fname, pname] = uiputfile(fullfile(cd, '*.eegdb'), 'Write the new EEG database to:');
    if (numel(fname)>1)
        eeg = []; eeg.work = struct([]); 
       writefilename = fullfile(pname, fname);
       disp('-----> get EEG data on all dates');
       [eeg, eegdata] = DataManager_FindAllEEG(datepath,eeg, copyover, eegdir);
       disp('-----> Set EEG parameters');
       [eeg, eegdata] = DataManager_FindEEGParm(eeg, eegdata);
       %%%compute initial variables
       eegind = 1:numel(eeg.general.eegfile); vv = 0;
       [eeg,eegdata] = DataManager_EEGComputeInit_Callback(eeg, eegdata, eegind, vv);
        
       %%%%set group list
       if (~isempty(eeg.general.eegfile))
          eegdata.grouplist.groupname{1} = 'List0'; eegdata.grouplist.groupindex{1} = 1:numel(eeg.general.eegfile);
          eegdata.grouplist.grouptype{1} = 'Manual'; eegdata.grouplist.groupcrit{1} = [];
          %%%Save the result
          save(writefilename, 'eeg', 'eegdata');
          %%plot reult    
          %hmain = gcbf; DataManager_PlotEEGDatabase(hmain, eeg, eegdata); 
          hmain = gcbf; DataManager_PlotSpikeDatabase(hmain, eeg, eegdata, 'EEG File', 'eegID', '.eegdb'); 
          set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
          setappdata(hmain, 'eeg', eeg); setappdata(hmain, 'eegdata', eegdata); eeg = []; eegdata = [];
       else
          disp(['-------------> no eeg files found!']);
       end
    end
end
disp('**********************');

