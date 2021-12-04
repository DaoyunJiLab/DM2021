function DataManager_SaveSelect_Callback
%%Save selected spikes into a new spike database
%%%%%%% Attention %%%%%%%%%%%%%
%%%%%Only cell arrays and non-empty numeric matrix assignments are allowed as entries such as:
%%%%%  pinfo.general.celltype{i}  or   pinfo.general.fieldnumber(i)

hf = gcbf; fname = get(hf, 'Name'); hgroup = getappdata(hf, 'hgroup');

groupselection = getappdata(hgroup, 'selection'); cellind = []; 
grpind = find(groupselection == 1); 

if (~isempty(grpind))
    mmind = strfind(fname, '__'); extfile = fname(mmind+2:numel(fname)); 
   [pp, nn, ee] = fileparts(extfile);
   if (strcmp(ee, '.behavdb'))
      pinfo = getappdata(hf, 'behav'); data = getappdata(hf, 'bhdata');
   elseif (strcmp(ee, '.eegdb'))
      pinfo = getappdata(hf, 'eeg'); data = getappdata(hf, 'eegdata');
   else %(strcmp(ee, '.spikedb'))
      pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
   end
   for (kk = 1:numel(grpind)) 
       if (~strncmpi(data.grouplist.groupname{grpind(kk)}, 'List0', 5))
          cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); 
       end
   end
   if (~isempty(cellind))
       disp('Save selected items into a new database');
       disp('-----> sorting data structure');
       pinfo = sortspike(pinfo, cellind);
       data = sortspike(data, cellind);
       if (~isempty(pinfo)) && (~isempty(data)) 
          disp('-----> writing selected group data to a file');
          [fname, pname] = uiputfile(fullfile(cd, ee), 'Write the new database to:');
          writefilename = fullfile(pname, fname);
          if (strcmp(ee, '.behavdb'))
             behav = pinfo; bhdata = data; save(writefilename, 'behav', 'bhdata'); %, '-mat', '-v7.3');
          elseif (strcmp(ee, '.eegdb'))
             eeg = pinfo; eegdata = data; save(writefilename, 'eeg', 'eegdata'); %, '-mat', '-v7.3');
          else %(strcmp(ee, '.spikedb'))
             save(writefilename, 'pinfo', 'data'); %, '-mat', '-v7.3');
          end
       else
          disp('-----> either selected contents or data structure is empty. Nothing to save.'); 
       end
   else
       disp('-----> no items in the selected groups or only List0 selected. Nothing to save.');
   end
else
   disp('-----> no groups selected. Nothing to save.');
end

disp('**********************************');

function pinfoout = sortspike(pinfo, selectindex)
%keep info if selected
%%%%%get all the fields in pinfo
pinfoout = []; fieldlist = fieldnames(pinfo); nfield = numel(fieldlist);
for (iii = 1:nfield)
    if (strncmpi(fieldlist{iii}, 'grouplist', 6))
        pinfoout.grouplist.groupname = {'List0'}; pinfoout.grouplist.groupindex = {[1:numel(selectindex)]};
        pinfoout.grouplist.grouptype = {'Manual'}; pinfoout.grouplist.groupcrit = {[]};
    else
        if ~strcmp(fieldlist{iii}, 'parentfile')
           subfieldlist = fieldnames(pinfo.(fieldlist{iii})); %%all the subfields
           if isempty(subfieldlist)
              work = struct([]); pinfoout.(fieldlist{iii}) = work;
           else
              for (jjj = 1:numel(subfieldlist))
               if isempty(pinfo.(fieldlist{iii}).(subfieldlist{jjj}))
                  pinfoout.(fieldlist{iii}).(subfieldlist{jjj}) = [];
               else
                  pinfoout.(fieldlist{iii}).(subfieldlist{jjj}) = pinfo.(fieldlist{iii}).(subfieldlist{jjj})(selectindex);
               end
              end
           end
        end
    end
end

