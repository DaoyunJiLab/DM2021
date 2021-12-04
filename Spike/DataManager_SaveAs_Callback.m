function DataManager_SaveAs_Callback
hf = gcbf; fname = get(hf, 'Name');
mmind = strfind(fname, '__'); extfile = fname(mmind+2:numel(fname)); 
[pp, nn, ee] = fileparts(extfile);
[MCroot, MCname, DAname, DEname, ANname, MAname] = CurrentVersion;
if (strcmp(ee, '.spikedb'))
   pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
   cpathname = fullfile(cd, '*.spikedb'); 
   [fname, pname] = uiputfile(cpathname, 'Spike database file to save:');
   if (fname ~= 0)
       filename = fullfile(pname, fname);
       s = whos('data'); Gb1 = s.bytes/(1024^3);
       s = whos('pinfo'); Gb2 = s.bytes/(1024^3);
       if Gb1 + Gb2 < 2
          save(filename, 'pinfo', 'data'); %disp('Small one');
       else
          save(filename, 'pinfo', 'data', '-mat', '-v7.3'); %disp('Big one');
       end
       set(hf, 'Name', strcat(MAname, '__', filename));
   end
elseif (strcmp(ee, '.behavdb'))
   behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
   cpathname = fullfile(cd, '*.behavdb'); 
   [fname, pname] = uiputfile(cpathname, 'Behav database file to save:');
   if (fname ~= 0)
       filename = fullfile(pname, fname);
       save(filename, 'behav', 'bhdata'); %, '-mat', '-v7.3');
       set(hf, 'Name', strcat(MAname, '__', filename));
   end
elseif (strcmp(ee, '.eegdb'))
   eeg = getappdata(hf, 'eeg'); eegdata = getappdata(hf, 'eegdata');
   cpathname = fullfile(cd, '*.eegdb'); 
   [fname, pname] = uiputfile(cpathname, 'EEG database file to save:');
   if (fname ~= 0)
       filename = fullfile(pname, fname);
       save(filename, 'eeg', 'eegdata'); %, '-mat', '-v7.3');
       set(hf, 'Name', strcat(MAname, '__', filename));
   end 
elseif (strcmp(ee, '.crrdb'))
   pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
   cpathname = fullfile(cd, '*.crrdb'); 
   [fname, pname] = uiputfile(cpathname, 'Crr database file to save:');
   if (fname ~= 0)
       filename = fullfile(pname, fname);
       s = whos('data'); Gb1 = s.bytes/(1024^3);
       s = whos('pinfo'); Gb2 = s.bytes/(1024^3);
       if Gb1 + Gb2 < 2
          save(filename, 'pinfo', 'data'); %disp('Small one');
       else
          save(filename, 'pinfo', 'data', '-mat', '-v7.3'); %disp('Big one');
       end
       set(hf, 'Name', strcat(MAname, '__', filename));
   end
elseif (strcmp(ee, '.seqdb'))
   pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
   cpathname = fullfile(cd, '*.seqdb'); 
   [fname, pname] = uiputfile(cpathname, 'Seq database file to save:');
   if (fname ~= 0)
       filename = fullfile(pname, fname);
       s = whos('data'); Gb1 = s.bytes/(1024^3);
       s = whos('pinfo'); Gb2 = s.bytes/(1024^3);
       if Gb1 + Gb2 < 2
          save(filename, 'pinfo', 'data'); %disp('Small one');
       else
          save(filename, 'pinfo', 'data', '-mat', '-v7.3'); %disp('Big one');
       end
       set(hf, 'Name', strcat(MAname, '__', filename));
   end
end

