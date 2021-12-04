function DataManager_SaveStripped_Callback
hf = gcbf; fname = get(hf, 'Name');
mmind = strfind(fname, '__'); extfile = fname(mmind+2:numel(fname)); 
[pp, nn, ee] = fileparts(extfile);
[MCroot, MCname, DAname, DEname, ANname, MAname] = CurrentVersion;
if (strcmp(ee, '.spikedb'))
   pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); 
   grouplist = data.grouplist; data = []; data.grouplist = grouplist;
   cpathname = fullfile(cd, '*.spikedb'); 
   [fname, pname] = uiputfile(cpathname, 'Spike database file to save:');
   if (fname ~= 0)
       filename = fullfile(pname, fname);
       save(filename, 'pinfo', 'data');
       %set(hf, 'Name', strcat(MAname, '__', filename));
   end
elseif (strcmp(ee, '.behavdb'))
   behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
   grouplist = data.grouplist; bhdata = []; bhdata.grouplist = grouplist;
   cpathname = fullfile(cd, '*.behavdb'); 
   [fname, pname] = uiputfile(cpathname, 'Behav database file to save:');
   if (fname ~= 0)
       filename = fullfile(pname, fname);
       save(filename, 'behav', 'bhdata');
       %set(hf, 'Name', strcat(MAname, '__', filename));
   end
elseif (strcmp(ee, '.eegdb'))
   eeg = getappdata(hf, 'eeg'); eegdata = getappdata(hf, 'eegdata');
   grouplist = eegdata.grouplist; eegdata = []; eegdata.grouplist = grouplist;
   cpathname = fullfile(cd, '*.eegdb'); 
   [fname, pname] = uiputfile(cpathname, 'EEG database file to save:');
   if (fname ~= 0)
       filename = fullfile(pname, fname);
       save(filename, 'eeg', 'eegdata');
       %set(hf, 'Name', strcat(MAname, '__', filename));
   end 
elseif (strcmp(ee, '.crrdb'))
   pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
   grouplist = data.grouplist; data = []; data.grouplist = grouplist;
   cpathname = fullfile(cd, '*.crrdb'); 
   [fname, pname] = uiputfile(cpathname, 'Crr database file to save:');
   if (fname ~= 0)
       filename = fullfile(pname, fname);
       save(filename, 'pinfo', 'data');
       %set(hf, 'Name', strcat(MAname, '__', filename));
   end
elseif (strcmp(ee, '.seqdb'))
   pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
   grouplist = data.grouplist; data = []; data.grouplist = grouplist;
   cpathname = fullfile(cd, '*.seqdb'); 
   [fname, pname] = uiputfile(cpathname, 'Seq database file to save:');
   if (fname ~= 0)
       filename = fullfile(pname, fname);
       save(filename, 'pinfo', 'data');
       %set(hf, 'Name', strcat(MAname, '__', filename));
   end
end

