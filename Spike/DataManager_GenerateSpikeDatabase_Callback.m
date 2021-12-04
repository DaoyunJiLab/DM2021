function DataManager_GenerateSpikeDatabase_Callback
%%orgainze an animal's all recording date and spikes
%%calculate features associate for each spike and set to corresponding fields
%%standard features are generated here
%%resulted structure is saved as a -mat type file with extension .spikedb
%%new features can be added later by openning the file

%Select multiple date directories to process
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
%%%%%animalpath = uigetdir(cd); %get just one animal dir
if (numel(datepath)>0)
    disp('Create new spike database');
    disp(['---> number of date final directories: ', num2str(numel(datepath))]);
    [fname, pname] = uiputfile(fullfile(cd, '*.spikedb'), 'Write the new spike database to:');
    if (numel(fname)>1)
       pinfo = []; pinfo.work = struct([]); 
       writefilename = fullfile(pname, fname);
       disp('-----> get all spikes');
         
       [pinfo, data] = DataManager_FindAllSpike(datepath, pinfo);
       %get parameters for all the calculations
     
       disp('-----> Set all parameters for all spikes');
       [pinfo, data] = DataManager_FindParm(pinfo, data);
       
       %%%compute initial variables
       vv = 0; cellind = 1:numel(pinfo.general.parmfile); 
       [pinfo,data] = DataManager_ComputeInit_Callback(pinfo, data, cellind, vv);
        
       %%%%set group list
       if (~isempty(pinfo.general.parmfile))
          data.grouplist.groupname{1} = 'List0'; data.grouplist.groupindex{1} = 1:numel(pinfo.general.parmfile);
          data.grouplist.grouptype{1} = 'Manual'; data.grouplist.groupcrit{1} = []; data.grouplist.groupparents{1} = [];
          data.parentfile = parentfile;
          %%%Save the result
          save(writefilename, 'pinfo', 'data');
          %%plot reult    
          hmain = gcbf; DataManager_PlotSpikeDatabase(hmain, pinfo, data); 
          set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
          setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); pinfo = []; data = [];
       else
          disp(['-------------> no cells found!']);
       end
    end
end
disp('**********************');

