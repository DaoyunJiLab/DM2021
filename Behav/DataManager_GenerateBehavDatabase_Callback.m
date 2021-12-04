function DataManager_GenerateBehavDatabase_Callback
%%orgainze an animal's all recording date and spikes
%%calculate features associate for each spike and set to corresponding fields
%%standard features are generated here
%%resulted structure is saved as a -mat type file with extension .spikedb
%%new features can be added later by openning the file

%get the animal name (animal root directory)
%animalpath = uigetdir(cd); %get just one animal dir
datepath = []; parentfile = [];
aaa = questdlg('Import a final-directory list?');
if (strcmp(aaa, 'Yes'))
   [fname, pname] = uigetfile(fullfile(cd, '*.lst'), 'Choose a list:');
   if (fname ~= 0)
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
      parentfile = fullfile(pname, fname);
   end
elseif (strcmp(aaa, 'No'))
   [datepath, ok] = DirDlg(cd);
end
datepath = unique(datepath);
if (numel(datepath)>0)
    disp('Create new behav database');
    %disp(['---> ', animalpath]);
    disp(['---> number of date final directories: ', num2str(numel(datepath))]);
    [fname, pname] = uiputfile(fullfile(cd, '*.behavdb'), 'Write the new behav database to:');
    if (numel(fname)>1)
        behav = []; behav.work = struct([]); 
       writefilename = fullfile(pname, fname);
       disp('-----> get behavioral data on all dates');
       [behav, bhdata] = DataManager_FindAllBehav(datepath, behav);
       disp('-----> Set behavioral parameters');
       [behav, bhdata] = DataManager_FindBehavParm(behav, bhdata);
       %%%compute initial variables
       dateind = 1:numel(behav.general.datedir); vv = 0; %%%dateind here is actually a session ind: multiple sessions on a given day
       [behav,bhdata] = DataManager_BehavComputeInit_Callback(behav, bhdata, dateind, vv);
        
       %%%%set group list
       if (~isempty(behav.general.datedir))
          bhdata.grouplist.groupname{1} = 'List0'; bhdata.grouplist.groupindex{1} = 1:numel(behav.general.datedir);
          bhdata.grouplist.grouptype{1} = 'Manual'; bhdata.grouplist.groupcrit{1} = []; bhdata.parentfile = parentfile;
          %%%Save the result
          save(writefilename, 'behav', 'bhdata');
          %%plot reult    
          hmain = gcbf; DataManager_PlotSpikeDatabase(hmain, behav, bhdata, 'Sessions', 'sessID', '.behavdb'); 
          set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
          setappdata(hmain, 'behav', behav); setappdata(hmain, 'bhdata', bhdata); behav = []; bhdata = [];
       else
          disp(['-------------> no cells found!']);
       end
    end
end
disp('**********************');

