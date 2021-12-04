function DataManager_SpikeSelect_Callback(varargin)
%%%%select spikes from spike database according to input criteria
%%%%output selected spikes to a group
%%%%criteria must contain full field name, relational operater, and values in sequence with spaces in between:
%%%%    behavior.fieldlength >= 10
%%%%    general.recarea = CA1
%%%%
hf = gcbf; dbtype = getappdata(hf, 'dbtype'); celllistvar = getappdata(hf, 'celllistvar');
%%%%%%%added options for all types of databases
pinfostr = 'pinfo'; datastr = 'data'; catok = 'general'; varok = 'parmfile';
switch dbtype
    case '.eegdb'
        pinfostr = 'eeg'; datastr = 'eegdata'; catok = 'general'; varok = 'eegfile'; %%also determine a variable (one in pinfo.general) to identify cells
    case '.behavdb'
        pinfostr = 'behav'; datastr = 'bhdata'; catok = 'general'; varok = 'sessID'; 
    case '.seq'
        catok = 'seq'; varok = 'seqID';  
    case '.seqdb'
        catok = 'general'; varok = 'tmpID';
end
pinfo = getappdata(hf, pinfostr); data = getappdata(hf, datastr); 
hspike = getappdata(hf, 'hspike'); hgroup= getappdata(hf, 'hgroup'); grouplistpos= getappdata(hf, 'grouplistpos');
spikeselection = getappdata(hspike, 'selection'); groupselection = getappdata(hgroup, 'selection');
nspike = numel(spikeselection); %number of spike (name)s
htext = getappdata(hspike, 'htext'); %every line (spike name) is a text object
ntext = numel(htext); %number of text lines displayed: not all spikes being displayed

%%%%need to modify the '.seqdb' varok
if strcmp(dbtype, '.seqdb')
    if ~isfield(pinfo.general, varok)
        varok = 'eventname';
    end
    if isfield(pinfo.parm, 'seqtype') && contains(pinfo.parm.seqtype{1}, 'evtitemized')
        varok = 'eventname';
    end
end
if (~isempty(varargin))
    tagmark = varargin{1};
else
    tagmark = get(gcbo, 'Tag');
end
if (strcmp(tagmark, 'spikeselectall')) %if just select all spikes in the table
   spikeselection = ones(1, nspike);
elseif (strcmp(tagmark, 'spikeselectallnot')) %if deselect all spikes in the table
   spikeselection = zeros(1, nspike);
elseif (strcmp(tagmark, 'showcrit')) %to show the group criteria 
   tit{1} = 'GroupName'; tit{2} = 'GroupType'; 
   val{1} = data.grouplist.groupname; val{2} = data.grouplist.grouptype; ngroup = numel(data.grouplist.groupname);
   maxcrit = 0;
   for (i = 1:ngroup)
       crit{i} = data.grouplist.groupcrit{i}; 
       if (maxcrit < numel(crit{i})) maxcrit = numel(crit{i}); end
   end
   for (j = 1:maxcrit)
       tit{2+j} = strcat('Crit', num2str(j));
       for (i = 1:ngroup)
           if (j <= numel(crit{i}))
               val{2+j}{i} = crit{i}{j};
           else
               val{2+j}{i} = [];
           end
       end
   end
   val = rearrangereal(val);
   hm = figure('Name', 'GroupDefinition', 'Units', 'inches', 'Position', [1 1 8 6], 'NumberTitle', 'off');
   pos = [0.2 0.2 7.6 5.6]; TextDisplayer_multiple(hm, pos, val, tit);
elseif (strcmp(tagmark, 'savegroup')) %to export/save group list
   grouplist = data.grouplist; 
   for ii = 1:numel(grouplist.grouptype) %%%% reset all to 'manual', auto is impossible to recreate
       grouplist.grouptype{ii} = 'Manual'; grouplist.groupcrit{ii} = [];
   end
   [nn, pp] = uiputfile(fullfile(cd, '*.glist'), 'File to save');
   save(fullfile(pp,nn), 'grouplist', '-mat');
elseif (strcmp(tagmark, 'importgroup')) %to import w grouplist and work it out
   [nn,pp] = uigetfile(fullfile(cd, '*.glist'), 'File to import groups');
   S = load(fullfile(pp,nn), '-mat'); grouplist = S.grouplist; 
   [newgrouplist, ok] = importall(pinfo, data, catok, varok, grouplist); 
   data.grouplist = newgrouplist;
   if ok
      %%%%%update the group list field
      delete(hgroup);
      hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
      setappdata(hf, 'hgroup', hgroup); setappdata(hf, datastr, data);
      %now re-set groupselectioin/spikeselection
      groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
      groupsel = find(groupselection == 1); spikeselectindex = [];
      for (i = 1:numel(groupsel))
         spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
      end
      spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
   end
elseif (strcmp(tagmark, 'outputevents')) %output time-related variables as event files to the .final\events\ directory
   hfield = getappdata(hf, 'hfield');
   spikeind = find(spikeselection == 1); %disp('What is going on?');
   outputevents(pinfo, data, spikeind, hfield);
elseif (strcmp(tagmark, 'deletegroup')) %to delete a group
   groupind = find(groupselection == 1);
   if (~isempty(groupind))
       dd = questdlg('Selected groups will be deleted. Are you sure?');
       if (strcmp(dd, 'Yes'))
           allind = 1:numel(data.grouplist.groupname); restind = setdiff(allind, groupind);
           data.grouplist.groupname = data.grouplist.groupname(restind); 
           data.grouplist.groupindex = data.grouplist.groupindex(restind);
           data.grouplist.grouptype = data.grouplist.grouptype(restind); 
           data.grouplist.groupcrit = data.grouplist.groupcrit(restind);
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, datastr, data);
           %now re-set groupselectioin/spikeselection
             groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
             groupsel = find(groupselection == 1); spikeselectindex = [];
             for (i = 1:numel(groupsel))
                 spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
             end
             spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
       else
           disp('-------> action canceled.');
       end
   else
        disp('-------> no group selected.');
   end
elseif (strcmp(tagmark, 'combinegroup')) %to combine selected group into another group
   groupsel = find(groupselection == 1); ind = []; ntype = 'Combine(';
   for (i = 1:numel(groupsel))
       ind = union(ind, data.grouplist.groupindex{groupsel(i)}); 
       if (i == numel(groupsel))
           endstr = ')';
       else
           endstr = '_';
       end
       ntype = strcat(ntype, data.grouplist.groupname{groupsel(i)}, endstr); 
   end
   if (isempty(ind))
       disp('-------> no spikes after the combination');
   else
       input = inputdlg({'Group name'}, 'Enter group name', 1);
       if (~isempty(input)) && (~isempty(input{1}))
           ngroup = numel(data.grouplist.groupname);
           data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = ind;
           data.grouplist.grouptype{ngroup+1} = ntype; 
           data.grouplist.groupcrit{ngroup+1} = [];
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, datastr, data);
           %now re-set groupselectioin/spikeselection
             groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
             groupsel = find(groupselection == 1); spikeselectindex = [];
             for (i = 1:numel(groupsel))
                 spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
             end
             spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
       else
           disp('-------> no group name entered.');
       end
   end
elseif (strcmp(tagmark, 'downsamplegroups')) %downsample groups to a certain number
   groupsel = find(groupselection == 1); ntype = 'Downsample('; ok = 1;
   for (i = 1:numel(groupsel))
       ind = data.grouplist.groupindex{groupsel(i)}; gname = data.grouplist.groupname{groupsel(i)};
       input = inputdlg({'Number of samples'; 'Group name'}, 'Enter sample number and group name', 2, {num2str(numel(ind)); gname});
       if (~isempty(input)) && (~isempty(input{2})) && (~isempty(input{1})) 
           nsample = str2num(input{1});
           if nsample < numel(ind)
              ntype = strcat(ntype, gname, ')');
              iii = randperm(numel(ind)); iii = iii(1:nsample);
              indnow = ind(iii);
              ngroup = numel(data.grouplist.groupname);
              data.grouplist.groupname{ngroup+1} = input{2}; data.grouplist.groupindex{ngroup+1} = indnow;
              data.grouplist.grouptype{ngroup+1} = ntype; 
              data.grouplist.groupcrit{ngroup+1} = [];
           else
              disp('-----------> warning: invalid number entered.'); ok = 0;
           end
       else
           disp('-------> no group name or invalid number entered.'); break; ok = 0;
       end
   end
   if ok
      %%%%%update the group list field
      delete(hgroup);
      hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
      setappdata(hf, 'hgroup', hgroup); setappdata(hf, datastr, data);
      %now re-set groupselectioin/spikeselection
      groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
      groupsel = find(groupselection == 1); spikeselectindex = [];
      for (i = 1:numel(groupsel))
           spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
      end
      spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
   end
elseif (strcmp(tagmark, 'subtractgroup')) %to subtract a small one from a large group
   groupsel = find(groupselection == 1); ind = [];
   if (numel(groupsel)<2)
       disp('-------> less or more than 2 groups selected');
   else
       ind1 = data.grouplist.groupindex{groupsel(1)}; ind2 = data.grouplist.groupindex{groupsel(2)};
       nam1 = data.grouplist.groupname{groupsel(1)}; nam2 = data.grouplist.groupname{groupsel(2)};
       if (numel(ind1)>numel(ind2))
           ind = setdiff(ind1, ind2);
       else
           ind = setdiff(ind2, ind1);
       end
       if (isempty(ind))
           disp('-------> no spikes after the subtraction');
       else
           input = inputdlg({'Group name'}, 'Enter group name', 1);
           if (~isempty(input)) && (~isempty(input{1}))
           ngroup = numel(data.grouplist.groupname);
           data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = ind;
           ntype = strcat('Subtract(', nam1, '_', nam2, ')');
           data.grouplist.grouptype{ngroup+1} = ntype; 
           data.grouplist.groupcrit{ngroup+1} = [];
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf,  datastr, data);
           %now re-set groupselectioin/spikeselection
             groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
             groupsel = find(groupselection == 1); spikeselectindex = [];
             for (i = 1:numel(groupsel))
                 spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
             end
             spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
           else
           disp('-------> no group name entered.');
           end
       end       
   end
elseif (strcmp(tagmark, 'intersectgroup')) %to choose the intersection of 2 or more groups
   groupsel = find(groupselection == 1); ind = []; ntype = 'Intersect(';
   if (numel(groupsel)<2)
       disp('-------> less than 2 groups selected');
   else
       ind = data.grouplist.groupindex{groupsel(1)}; ntype = strcat(ntype, data.grouplist.groupname{groupsel(1)}, '_');
       for (ti = 2:numel(groupsel))
           indnow = data.grouplist.groupindex{groupsel(ti)}; ind = intersect(ind, indnow);
           if (ti == numel(groupsel))
               endstr = ')';
           else
               endstr = '_';
           end
           ntype = strcat(ntype, data.grouplist.groupname{groupsel(ti)}, endstr); 
       end
       if (isempty(ind))
           disp('-------> no spikes after the reverse intersection');
       else
           input = inputdlg({'Group name'}, 'Enter group name', 1);
           if (~isempty(input)) && (~isempty(input{1}))
           ngroup = numel(data.grouplist.groupname);
           data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = ind;
           data.grouplist.grouptype{ngroup+1} = ntype; 
           data.grouplist.groupcrit{ngroup+1} = [];
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, datastr, data);
           %now re-set groupselectioin/spikeselection
             groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
             groupsel = find(groupselection == 1); spikeselectindex = [];
             for (i = 1:numel(groupsel))
                 spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
             end
             spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
           else
           disp('-------> no group name entered.');
           end
       end       
   end
elseif (strcmp(tagmark, 'reverseintgroup')) %to subtract the largest group from the intersection of 2 or more groups
   groupsel = find(groupselection == 1); ind = []; ntype = 'ReverseInt(';
   if (numel(groupsel)<2)
       disp('-------> less than 2 groups selected');
   else
       indint = data.grouplist.groupindex{groupsel(1)}; indnum(1) = numel(data.grouplist.groupindex{groupsel(1)}); 
       ntype = strcat(ntype, data.grouplist.groupname{groupsel(1)}, '_');
       for (ti = 2:numel(groupsel))
           indnow = data.grouplist.groupindex{groupsel(ti)}; indnum(ti) = numel(indnow);
           indint = intersect(indint, indnow); 
           if (ti == numel(groupsel))
               endstr = ')';
           else
               endstr = '_';
           end
           ntype = strcat(ntype, data.grouplist.groupname{groupsel(ti)}, endstr);
       end
       [~, iij] = max(indnum); indall = data.grouplist.groupindex{groupsel(iij)};
       ind = setdiff(indall, indint);
       if (isempty(ind))
           disp('-------> no spikes after the reverse intersection');
       else
           input = inputdlg({'Group name'}, 'Enter group name', 1);
           if (~isempty(input)) && (~isempty(input{1}))
           ngroup = numel(data.grouplist.groupname);
           data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = ind;
           data.grouplist.grouptype{ngroup+1} = ntype; 
           data.grouplist.groupcrit{ngroup+1} = [];
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, datastr, data);
           %now re-set groupselectioin/spikeselection
             groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
             groupsel = find(groupselection == 1); spikeselectindex = [];
             for (i = 1:numel(groupsel))
                 spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
             end
             spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
           else
           disp('-------> no group name entered.');
           end
       end       
   end
elseif (strcmp(tagmark, 'manuselect')) %if set the selected spikes into a group
   spikeselectind = find(spikeselection == 1);
   if (~isempty(spikeselectind))
       input = inputdlg({'Group name'}, 'Enter group name', 1);
       if (~isempty(input)) && (~isempty(input{1}))
           ngroup = numel(data.grouplist.groupname);
           data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = spikeselectind;
           data.grouplist.grouptype{ngroup+1} = 'Manual'; data.grouplist.groupcrit{ngroup+1} = [];
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, datastr, data);
           %now re-set groupselectioin/spikeselection
             groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
             groupsel = find(groupselection == 1); spikeselectindex = [];
             for (i = 1:numel(groupsel))
                 spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
             end
             spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
       else
           disp('-------> no group name entered.');
       end
   else
       disp('-------> no spikes selected')
   end
elseif (strcmp(tagmark, 'redoone'))
   groupsel = find(groupselection == 1);
   if (numel(groupsel)~=1)
       disp('-------> None or more than 1 group selected; aborted');
   else
       [data, ok] = redosingle(pinfo, data, groupsel, spikeselection, catok, varok); 
       if ok
          %%%%%update the group list field
          delete(hgroup);
          hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
          setappdata(hf, 'hgroup', hgroup); setappdata(hf, datastr, data);
          %now re-set groupselectioin/spikeselection
          groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
          groupsel = find(groupselection == 1); spikeselectindex = [];
          for (i = 1:numel(groupsel))
              spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
          end
          spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
       end
   end
elseif (strcmp(tagmark, 'redoall'))
   [data, ok] = redoall(pinfo, data, catok, varok); 
   if ok
      %%%%%update the group list field
      delete(hgroup);
      hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
      setappdata(hf, 'hgroup', hgroup); setappdata(hf, datastr, data);
      %now re-set groupselectioin/spikeselection
      groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
      groupsel = find(groupselection == 1); spikeselectindex = [];
      for (i = 1:numel(groupsel))
         spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
      end
      spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
   end
elseif (strcmp(tagmark, 'autoselect')) %if select spikes by criterion
   ncrt = 1;
   for (kkk = 1:ncrt)
       s1 = {'Criterion #1:'; 'Criterion #2:'; 'Criterion #3:'; 'Criterion #4:'; 'Criterion #5:';}; %select all spikes that satisfy all the criteria
       input = inputdlg(s1, 'Criterion input', 5);
       for (i = 1:numel(input))
           if (~isempty(input{i})) 
              crt{kkk}{i} = input{i};
           end
       end
   end 
       %%%autoselect within the selected group
       pgroup = []; %%%parent groups
       groupsel = find(groupselection == 1); cellind = [];
       for (tti = 1:numel(groupsel))
           cellind = union(cellind, data.grouplist.groupindex{groupsel(tti)});
           pgroup = strcat(pgroup, '_', data.grouplist.groupname{groupsel(tti)});
       end
       parmfilenow = cell(1,ncrt);
       for (kkk = 1:ncrt)
           parmfilenow{kkk} = DataManager_SearchSpikeDatabase(pinfo, crt{kkk}, catok, varok, cellind); %disp(parmfilenow{kkk}{1})
       end
       parmfile = findcommoncell(parmfilenow);
       disp(['-----> ', num2str(numel(parmfile)), ' cells match the criteria']);
       spikesel = 0*spikeselection;
       if (~isempty(parmfile))
          %mark spikes in pinfo that match the spikelist: only need to match one of the filenames: s1file, runfile or s2file
          for (i = 1:nspike)
              if (~isempty(find(strcmp(parmfile, pinfo.(catok).(varok){i}))))
                 spikesel(i) = 1;
              end
          end
       end
       spikeindnow = find(spikesel == 1);
       if (~isempty(spikeindnow))
          input = inputdlg({'Group name'}, 'Enter group name', 1);
          if (~isempty(input)) && (~isempty(input{1}))
             ngroup = numel(data.grouplist.groupname);
             data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = spikeindnow;
             data.grouplist.grouptype{ngroup+1} = strcat('Auto(', pgroup, ')'); data.grouplist.groupcrit{ngroup+1} = crt;
             %%%%%update the group list field
             delete(hgroup);
             hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
             setappdata(hf, 'hgroup', hgroup); setappdata(hf, datastr, data);
             %now re-set groupselectioin/spikeselection
             groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
             groupsel = find(groupselection == 1); spikeselectindex = [];
             for (i = 1:numel(groupsel))
                 spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
             end
             spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
          else
             disp('-------> no group name entered.');
          end
       else
          disp('-------> no spikes selected')
       end
end
for (i = 1:ntext)
   tagtag = get(htext(i), 'Tag');
   [str, rem] = strtok(tagtag, '_');
   linenum = str2num(str); %current line number that selected
   if (spikeselection(linenum) == 0)
       set(htext(i), 'Color', [0 0 0]);
   else
       set(htext(i), 'Color', [1 0 0]);
   end
end
setappdata(hspike, 'selection', spikeselection); setappdata(hspike, 'htext', htext);
disp('*************');

function [newgrouplist, ok] = importall(pinfo, data, catok, varok, grouplist)
ok = 1;
ngroup = numel(grouplist.groupname);
newgrouplist.groupname{1} = 'List0'; newgrouplist.groupindex{1} = 1:numel(pinfo.(catok).(varok));
newgrouplist.grouptype{1} = 'Manual'; newgrouplist.groupcrit{1} = []; newgrouplist.groupparents{1} = [];
if ~strcmp(grouplist.groupname{1}, 'List0')
    disp(['---------> warning: imported first group is not List0']); startg = 1;
else
    startg = 2;
end
for (i = startg:ngroup)
    [newgrouplist, ok] = importsingle_auto(pinfo, data, i, catok, varok, newgrouplist, grouplist);
end

function [newgrouplist, ok] = importsingle_auto(pinfo, data, groupsel, catok, varok, newgrouplist, grouplist)
ok = 1;
ntype = grouplist.grouptype{groupsel}; 
oldcrit = grouplist.groupcrit{groupsel};
if strncmpi(ntype, 'Manual', 6)%%%%%now changed to appending new groups for manual type
   %ii = find(strcmp(data.grouplist.groupname, grouplist.groupname{groupsel}));
   %if isempty(ii)
   %   disp(['-------------> skip manual group: ', grouplist.groupname{groupsel}]);
   %else
      disp(['-------------> keep manual group: ', grouplist.groupname{groupsel}]);
      groupnow = 1+numel(newgrouplist.groupname);
      newgrouplist.groupindex{groupnow} = grouplist.groupindex{groupsel}; 
      newgrouplist.groupname{groupnow} = grouplist.groupname{groupsel}; newgrouplist.groupcrit{groupnow} = grouplist.groupcrit{groupsel};
      newgrouplist.grouptype{groupnow} = grouplist.grouptype{groupsel}; newgrouplist.groupcrit{groupnow} = grouplist.groupcrit{groupsel}; 
      if isfield(grouplist, 'groupparents') && (groupsel <= numel(grouplist.groupparents))
         newgrouplist.groupparents{groupnow} = grouplist.groupparents{groupsel};
      else
         newgrouplist.groupparents{groupnow} = [];
      end
   %end
elseif ~isempty(strfind(ntype, 'Auto'))
           %%%%%%%%%working here: figure out pgroupind, oldcrit, groupind
           ii = strfind(ntype, '_'); ng = numel(ii); pgroupind = zeros(1, ng);
           if (ng > 0)
              for (tt = 1:ng-1)
               pgname = ntype(ii(tt)+1:ii(tt+1)-1); iii = find( strcmp(newgrouplist.groupname, pgname) ); pgroupind(tt) = iii(1);
              end
              pgname = ntype(ii(ng)+1:numel(ntype)-1); iii = find( strcmp(newgrouplist.groupname, pgname) ); pgroupind(ng) = iii(1);
           end
           if ~isempty(find( (pgroupind >= groupsel)|(pgroupind==0) ))
               ok = 1; disp(['-------------> parent group(s) not found or behind current group: ', grouplist.groupname{groupsel}, '; aborted']);
           else
               groupnow = 1+numel(newgrouplist.groupname);
               groupind = autoselectnow_auto_import(pinfo, newgrouplist, pgroupind, oldcrit, catok, varok);
               newgrouplist.groupindex{groupnow} = groupind; 
               newgrouplist.groupname{groupnow} = grouplist.groupname{groupsel}; newgrouplist.groupcrit{groupnow} = grouplist.groupcrit{groupsel};
               newgrouplist.grouptype{groupnow} = grouplist.grouptype{groupsel}; newgrouplist.groupcrit{groupnow} = grouplist.groupcrit{groupsel}; 
               if isfield(grouplist, 'groupparents') && (groupsel <= numel(grouplist.groupparents))
                  newgrouplist.groupparents{groupnow} = grouplist.groupparents{groupsel};
               else
                  newgrouplist.groupparents{groupnow} = [];
               end
           end
else %%%if subtract, intersect, ect.
    ii = strfind(ntype, '('); 
    if (~isempty(ii)) 
        keyword = ntype(1:ii(1)-1); otherword = ntype(ii(1)+1:numel(ntype)-1);
        ii = strfind(otherword, '_'); ng = numel(ii); pgind = zeros(1, ng+1); 
        if (ng > 0)
            pgname = otherword(1:ii(1)-1); iii = find( strcmp(newgrouplist.groupname, pgname) ); pgind(1) = iii(1);
            for (ti = 1:ng-1)
                pgname = otherword(ii(ti)+1:ii(ti+1)-1); iii = find( strcmp(newgrouplist.groupname, pgname) ); pgind(ti+1) = iii(1);
            end
            pgname = otherword(ii(ng)+1:numel(otherword)); iii = find( strcmp(newgrouplist.groupname, pgname) ); pgind(ng+1) = iii(1);
            if ~isempty(find( (pgind >= groupsel)|(pgind==0) ))
               ok = 1; disp(['-------------> parent group(s) not found or behind current group: ', grouplist.groupname{groupsel}, '; aborted']);
            else
               groupnow = 1+numel(newgrouplist.groupname);
               groupind = keywordgroup_import(newgrouplist, keyword, pgind);
               newgrouplist.groupindex{groupnow} = groupind; 
               newgrouplist.groupname{groupnow} = grouplist.groupname{groupsel}; newgrouplist.groupcrit{groupnow} = grouplist.groupcrit{groupsel};
               newgrouplist.grouptype{groupnow} = grouplist.grouptype{groupsel}; newgrouplist.groupcrit{groupnow} = grouplist.groupcrit{groupsel}; 
               if isfield(grouplist, 'groupparents') && (groupsel <= numel(grouplist.groupparents))
                  newgrouplist.groupparents{groupnow} = grouplist.groupparents{groupsel};
               else
                  newgrouplist.groupparents{groupnow} = [];
               end
            end
        else
            ok = 0; disp(['-------------> can not determine parent name(s): ', grouplist.groupname{groupsel}, '; aborted']);
        end
    else
        ok = 0; disp(['-------------> can not determine group type: ', grouplist.groupname{groupsel}, '; aborted']);
    end
end

function [data, ok] = redoall(pinfo, data, catok, varok)
ok = 1;
ngroup = numel(data.grouplist.groupname);
if ~strcmp(data.grouplist.groupname{1}, 'List0')
    disp(['---------> warning: first group is not List0']);
end
for (i = 2:ngroup)
    [data, ok] = redosingle_auto(pinfo, data, i, catok, varok);
end

function [data,ok] = redosingle_auto(pinfo, data, groupsel, catok, varok)
ok = 1;
ntype = data.grouplist.grouptype{groupsel}; 
oldcrit = data.grouplist.groupcrit{groupsel};
if strncmpi(ntype, 'Manual', 6)
   disp(['-------------> skipping manual group: ', data.grouplist.groupname{groupsel}]);
elseif ~isempty(strfind(ntype, 'Auto'))
           %%%%%%%%%working here: figure out pgroupind, oldcrit, groupind
           ii = strfind(ntype, '_'); ng = numel(ii); pgroupind = zeros(1, ng);
           if (ng > 0)
              for (tt = 1:ng-1)
               pgname = ntype(ii(tt)+1:ii(tt+1)-1); iii = find( strcmp(data.grouplist.groupname, pgname) ); pgroupind(tt) = iii(1);
              end
              pgname = ntype(ii(ng)+1:numel(ntype)-1); iii = find( strcmp(data.grouplist.groupname, pgname) ); pgroupind(ng) = iii(1);
           end
           if ~isempty(find( (pgroupind >= groupsel)|(pgroupind==0) ))
               ok = 1; disp(['-------------> parent group(s) not found or behind current group: ', data.grouplist.groupname{groupsel}, '; aborted']);
           else
               data = autoselectnow_auto(pinfo, data, pgroupind, oldcrit, groupsel, catok, varok);
           end
else %%%if subtract, intersect, ect.
    ii = strfind(ntype, '('); 
    if (~isempty(ii)) 
        keyword = ntype(1:ii(1)-1); otherword = ntype(ii(1)+1:numel(ntype)-1);
        ii = strfind(otherword, '_'); ng = numel(ii); pgind = zeros(1, ng+1); 
        if (ng > 0)
            pgname = otherword(1:ii(1)-1); iii = find( strcmp(data.grouplist.groupname, pgname) ); pgind(1) = iii(1);
            for (ti = 1:ng-1)
                pgname = otherword(ii(ti)+1:ii(ti+1)-1); iii = find( strcmp(data.grouplist.groupname, pgname) ); pgind(ti+1) = iii(1);
            end
            pgname = otherword(ii(ng)+1:numel(otherword)); iii = find( strcmp(data.grouplist.groupname, pgname) ); pgind(ng+1) = iii(1);
            if ~isempty(find( (pgind >= groupsel)|(pgind==0) ))
               ok = 1; disp(['-------------> parent group(s) not found or behind current group: ', data.grouplist.groupname{groupsel}, '; aborted']);
            else
               data = keywordgroup(data, keyword, pgind, groupsel);
            end
        else
            ok = 0; disp(['-------------> can not determine parent name(s): ', data.grouplist.groupname{groupsel}, '; aborted']);
        end
    else
        ok = 0; disp(['-------------> can not determine group type: ', data.grouplist.groupname{groupsel}, '; aborted']);
    end
end

function [data, ok] = redosingle(pinfo, data, groupsel, spikeselection, catok, varok)
ntype = data.grouplist.grouptype{groupsel}; ok = 1;
oldcrit = data.grouplist.groupcrit{groupsel};
if strncmpi(ntype, 'Manual', 6)
   spikeselectind = find(spikeselection == 1);
   if (~isempty(spikeselectind))
              data.grouplist.groupindex{groupsel} = spikeselectind;
   else
              disp('-------> no spikes selected')
   end
elseif ~isempty(strfind(ntype, 'Auto'))
           %%%%%%%%%working here: figure out pgroupind, oldcrit, groupind
           ii = strfind(ntype, '_'); ng = numel(ii); pgroupind = zeros(1, ng);
           if (ng > 0)
              for (tt = 1:ng-1)
               pgname = ntype(ii(tt)+1:ii(tt+1)-1); iii = find( strcmp(data.grouplist.groupname, pgname) ); pgroupind(tt) = iii(1);
              end
              pgname = ntype(ii(ng)+1:numel(ntype)-1); iii = find( strcmp(data.grouplist.groupname, pgname) ); pgroupind(ng) = iii(1);
              data = autoselectnow(pinfo, data, pgroupind, oldcrit, groupsel, catok, varok);
           end
else %%%if subtract, intersect, ect.
    ii = strfind(ntype, '('); 
    if (~isempty(ii)) 
        keyword = ntype(1:ii(1)-1); otherword = ntype(ii(1)+1:numel(ntype)-1);
        ii = strfind(otherword, '_'); ng = numel(ii); pgind = zeros(1, ng+1); 
        if (ng > 0)
            pgname = otherword(1:ii(1)-1); iii = find( strcmp(data.grouplist.groupname, pgname) ); pgind(1) = iii(1);
            for (ti = 1:ng-1)
                pgname = otherword(ii(ti)+1:ii(ti+1)-1); iii = find( strcmp(data.grouplist.groupname, pgname) ); pgind(ti+1) = iii(1);
            end
            pgname = otherword(ii(ng)+1:numel(otherword)); iii = find( strcmp(data.grouplist.groupname, pgname) ); pgind(ng+1) = iii(1);
            data = keywordgroup(data, keyword, pgind, groupsel);
        end
    end
end

function ind = keywordgroup_import(grouplist, keyword, pgind)
tagmark = strcat(lower(keyword), 'group');
if (strcmp(tagmark, 'combinegroup')) %to combine selected group into another group
   ind = [];
   for (ti = 1:numel(pgind))
        ind = union(ind, grouplist.groupindex{pgind(ti)}); 
   end
elseif (strcmp(tagmark, 'subtractgroup')) %to subtract a small one from a large group
   ind1 = grouplist.groupindex{pgind(1)}; ind2 = grouplist.groupindex{pgind(2)};
   if (numel(ind1)>numel(ind2))
       ind = setdiff(ind1, ind2);
   else
       ind = setdiff(ind2, ind1);
   end
elseif (strcmp(tagmark, 'intersectgroup')) %to choose the intersection of 2 or more groups
   ind = grouplist.groupindex{pgind(1)}; 
   for (ti = 2:numel(pgind))
        indnow = grouplist.groupindex{pgind(ti)}; ind = intersect(ind, indnow);
   end
elseif (strcmp(tagmark, 'reverseintgroup')) %to subtract the largest one from intersection of 2 or more groups
   indint = grouplist.groupindex{pgind(1)}; indnum(1) = numel(grouplist.groupindex{pgind(1)}); 
   for (ti = 2:numel(pgind))
       indnow = grouplist.groupindex{pgind(ti)}; indnum(ti) = numel(indnow);
       indint = intersect(indint, indnow);
   end
   [~, iij] = max(indnum); indall = grouplist.groupindex{pgind(iij)};
   ind = setdiff(indall, indint);
end

function data = keywordgroup(data, keyword, pgind, groupsel)
tagmark = strcat(lower(keyword), 'group');
if (strcmp(tagmark, 'combinegroup')) %to combine selected group into another group
   ind = [];
   for (ti = 1:numel(pgind))
        ind = union(ind, data.grouplist.groupindex{pgind(ti)}); 
   end
   data.grouplist.groupindex{groupsel} = ind;
elseif (strcmp(tagmark, 'subtractgroup')) %to subtract a small one from a large group
   ind1 = data.grouplist.groupindex{pgind(1)}; ind2 = data.grouplist.groupindex{pgind(2)};
   if (numel(ind1)>numel(ind2))
       ind = setdiff(ind1, ind2);
   else
       ind = setdiff(ind2, ind1);
   end
   data.grouplist.groupindex{groupsel} = ind;
elseif (strcmp(tagmark, 'intersectgroup')) %to choose the intersection of 2 or more groups
   ind = data.grouplist.groupindex{pgind(1)}; 
   for (ti = 2:numel(pgind))
        indnow = data.grouplist.groupindex{pgind(ti)}; ind = intersect(ind, indnow);
   end
   data.grouplist.groupindex{groupsel} = ind;
elseif (strcmp(tagmark, 'reverseintgroup')) %to choose the intersection of 2 or more groups
   indint = grouplist.groupindex{pgind(1)}; indnum(1) = numel(grouplist.groupindex{pgind(1)}); 
   for (ti = 2:numel(pgind))
       indnow = grouplist.groupindex{pgind(ti)}; indnum(ti) = numel(indnow);
       indint = intersect(indint, indnow);
   end
   [~, iij] = max(indnum); indall = grouplist.groupindex{pgind(iij)};
   ind = setdiff(indall, indint); 
%    indint = data.grouplist.groupindex{pgind(1)}; indall = data.grouplist.groupindex{pgind(1)}; 
%    for (ti = 2:numel(pgind))
%        indnow = data.grouplist.groupindex{pgind(ti)};
%        indint = intersect(indint, indnow); indall = union(indall, indnow);
%    end
%    ind = intersect(indall, indint);
   data.grouplist.groupindex{groupsel} = ind;
end

function groupind = autoselectnow_auto_import(pinfo, newgrouplist, pgroupind, crt, catok, varok)
ncrt = 1;
%%%autoselect within the selected group
groupsel = pgroupind; cellind = [];
for (tti = 1:numel(groupsel))
     cellind = union(cellind, newgrouplist.groupindex{groupsel(tti)});
end
parmfilenow = cell(1,ncrt);
for (kkk = 1:ncrt)
      parmfilenow{kkk} = DataManager_SearchSpikeDatabase(pinfo, crt{kkk}, catok, varok, cellind);
end
parmfile = findcommoncell(parmfilenow);
nspike = numel(pinfo.(catok).(varok)); spikesel = zeros(1, nspike); 
if (~isempty(parmfile))
          %mark spikes in pinfo that match the spikelist: only need to match one of the filenames: s1file, runfile or s2file
          for (i = 1:nspike)
              if (~isempty(find(strcmp(parmfile, pinfo.(catok).(varok){i}))))
                 spikesel(i) = 1;
              end
          end
end
groupind = find(spikesel == 1);

function data = autoselectnow_auto(pinfo, data, pgroupind, crt, groupind, catok, varok) 
ncrt = 1;
%%%autoselect within the selected group
groupsel = pgroupind; cellind = [];
for (tti = 1:numel(groupsel))
     cellind = union(cellind, data.grouplist.groupindex{groupsel(tti)});
end
parmfilenow = cell(1,ncrt);
for (kkk = 1:ncrt)
      parmfilenow{kkk} = DataManager_SearchSpikeDatabase(pinfo, crt{kkk}, catok, varok, cellind);
end
parmfile = findcommoncell(parmfilenow);
nspike = numel(pinfo.(catok).(varok)); spikesel = zeros(1, nspike); 
if (~isempty(parmfile))
          %mark spikes in pinfo that match the spikelist: only need to match one of the filenames: s1file, runfile or s2file
          for (i = 1:nspike)
              if (~isempty(find(strcmp(parmfile, pinfo.(catok).(varok){i}))))
                 spikesel(i) = 1;
              end
          end
end
spikeindnow = find(spikesel == 1);
data.grouplist.groupindex{groupind} = spikeindnow;

function data = autoselectnow(pinfo, data, pgroupind, oldcrit, groupind, catok, varok) 
ncrt = 1;
for (kkk = 1:ncrt)
    def = [];
      s1 = {'Criterion #1:'; 'Criterion #2:'; 'Criterion #3:'; 'Criterion #4:'; 'Criterion #5:';}; %select all spikes that satisfy all the criteria
      for (kt = 1:5)
          if (kt<=numel(oldcrit{kkk}))
             def{kt} = oldcrit{kkk}{kt};
          else
             def{kt} = '';
          end
      end
      input = inputdlg(s1, 'Criterion input', 5, def);
      for (i = 1:numel(input))
           if (~isempty(input{i})) 
              crt{kkk}{i} = input{i};
           end
       end
end 
%%%autoselect within the selected group
groupsel = pgroupind; cellind = [];
for (tti = 1:numel(groupsel))
           cellind = union(cellind, data.grouplist.groupindex{groupsel(tti)});
end
parmfilenow = cell(1,ncrt);
for (kkk = 1:ncrt)
           parmfilenow{kkk} = DataManager_SearchSpikeDatabase(pinfo, crt{kkk}, catok, varok, cellind);
end
parmfile = findcommoncell(parmfilenow);
disp(['-----> ', num2str(numel(parmfile)), ' cells match the criteria']);
nspike = numel(pinfo.(catok).(varok)); spikesel = zeros(1, nspike); 
if (~isempty(parmfile))
          %mark spikes in pinfo that match the spikelist: only need to match one of the filenames: s1file, runfile or s2file
          for (i = 1:nspike)
              if (~isempty(find(strcmp(parmfile, pinfo.(catok).(varok){i}))))
                 spikesel(i) = 1;
              end
          end
end
spikeindnow = find(spikesel == 1);
if (~isempty(spikeindnow))
             data.grouplist.groupindex{groupind} = spikeindnow;
             data.grouplist.groupcrit{groupind} = crt;
else
    disp('-------> no spikes selected')
end

function parmfile = findcommoncell(parmfilenow)
%%%find cells common in indexnow{1}, indexnow{2}, ...
nset = numel(parmfilenow);
parmfile = []; ncl = 0;
for (n = 1:numel(parmfilenow{1}))
    tt = 1;
    for (k = 1:nset)
        if (isempty(find(strcmp(parmfilenow{k}, parmfilenow{1}{n}))))
            tt = 0; break;
        end
    end
    if (tt == 1)
        ncl = ncl + 1; parmfile{ncl} = parmfilenow{1}{n};
    end
end

function savegrouplist(pinfo,data,groupsel, hfield, celllistvar)
%determine all possible field
fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle);
for (i = 1:nfield) subfield{i} = fieldnames(pinfo.(fieldtitle{i})); end
%determine selected groups and selected fields
numgroup = numel(groupsel);
for (xx = 1:numel(groupsel)) %if any groups selected
    spikeselectindex = data.grouplist.groupindex{groupsel(xx)};
    valuecol = []; valuetitle = []; ncol = 1;
    valuetitle{ncol} = data.grouplist.groupname{groupsel(xx)};
    valuecol{ncol} = pinfo.general.(celllistvar)(spikeselectindex); %first column = selected spike names
    for (i = 1:nfield) %search selected sub-fields in each big field
        fieldselection = getappdata(hfield(i), 'selection');
        fieldselectindex = find( fieldselection == 1);
        for (j = 1:numel(fieldselectindex) )
            nj = fieldselectindex(j);
            ncol = ncol + 1;
            valuetitle{ncol} = subfield{i}{nj};
            kkk = numel(pinfo.(fieldtitle{i}).(subfield{i}{nj}));
            for (mm = 1:numel(spikeselectindex))
                if (spikeselectindex(mm) > kkk)
                    valuecol{ncol}{mm} = [];
                elseif (iscell(pinfo.(fieldtitle{i}).(subfield{i}{nj})))
                    wield = pinfo.(fieldtitle{i}).(subfield{i}{nj}){spikeselectindex(mm)}; %this is a cell anyway
                    if ( isnumeric (wield) )
                       for (mnk = 1:numel(wield))
                           valuecol{ncol}{mm}{mnk} = num2str( wield(mnk), '%20.10f' );
                       end
                    else
                       valuecol{ncol}{mm} = wield;     
                    end
                else
                    valuecol{ncol}{mm} = pinfo.(fieldtitle{i}).(subfield{i}{nj})(spikeselectindex(mm));
                end
            end
        end
    end
    %re-arrange valuecol{ncol}{mm}{mnk} to fit textmsg{ncol}{mm}
    textmsg = rearrange(valuecol); 
    %%%re-arrange values into writeable strings
    nvar = numel(textmsg); nline = numel(textmsg{1});
    for (tti = 1:nvar)
        for (ttj = 1:nline)
            if (~isstr(textmsg{tti}{ttj}))
                textmsg{tti}{ttj} = num2str(textmsg{tti}{ttj}, '%20.10f');
            end
         end
    end
    
    %%%put columns together    
    tline = cell(nline+1, 1);
    for (tti = 1:nvar)
        tline{1} =[tline{1}, '    ', valuetitle{tti}];
        for (ttj = 1:nline)
            tline{ttj+1} = [tline{ttj+1}, '    ', textmsg{tti}{ttj}];
        end
    end
    %%%write
    filenow = fullfile(cd, strcat(data.grouplist.groupname{groupsel(xx)}, '.glist'));
    fid = fopen(filenow, 'wt');
    for (ttj = 1:nline+1)
       fprintf(fid, '%s \n', tline{ttj}); 
    end
    fclose(fid);
end

function textmsg = rearrange(valuecol)
%rearrnge valuecol{ncol}{mline}{nitem} to textmsg{nn}{nline}
ncol = numel(valuecol);
for (i = 1:ncol) %initial assignment
    textmsg{i} = [];
    nline(i) = 0;
    linenum(i) = numel(valuecol{i});
end
nspike = max(linenum); %all line number should be same(= spike number) for all columns, take the maximum anyway
for (i = 1:ncol)
for (j = 1:nspike)
    if ( j > numel(valuecol{i}) )
        itemnum(i,j) = 0;
    elseif (iscell(valuecol{i}{j}))
        itemnum(i,j) = numel(valuecol{i}{j});
    else 
        itemnum(i,j) = 1;
    end
end   
end

for (j = 1:nspike)
    maxitem(j) = max(itemnum(:,j));
end
for (i = 1:ncol)
    for (j = 1:nspike)
        for (k = 1:maxitem(j))
            nline(i) = nline(i) + 1;
            if (itemnum(i,j) == 0)
                textmsg{i}{nline(i)} = 'ND'; %stuffing
            elseif ( k > itemnum(i,j))
                textmsg{i}{nline(i)} = 'ND'; %stuffing
            elseif (iscell(valuecol{i}{j}))
                textmsg{i}{nline(i)} = valuecol{i}{j}{k};
            else
                textmsg{i}{nline(i)} = valuecol{i}{j};
            end
        end
    end
end

function textmsg = rearrangereal(valuecol)
%rearrnge valuecol{ncol}{mline}{nitem} to textmsg{nn}{nline}
ncol = numel(valuecol);
for (i = 1:ncol) %initial assignment
    textmsg{i} = [];
    nline(i) = 0;
    linenum(i) = numel(valuecol{i});
end
nspike = max(linenum); %all line number should be same(= spike number) for all columns, take the maximum anyway
for (i = 1:ncol)
for (j = 1:nspike)
    if ( j > numel(valuecol{i}) )
        itemnum(i,j) = 0;
    elseif (iscell(valuecol{i}{j}))
        itemnum(i,j) = numel(valuecol{i}{j});
    else 
        itemnum(i,j) = 1;
    end
end   
end

for (j = 1:nspike)
    maxitem(j) = max(itemnum(:,j));
end
for (i = 1:ncol)
    for (j = 1:nspike)
        for (k = 1:maxitem(j))
            nline(i) = nline(i) + 1;
            if (itemnum(i,j) == 0)
                textmsg{i}{nline(i)} = ' '; %stuffing
            elseif ( k > itemnum(i,j))
                textmsg{i}{nline(i)} = ' '; %stuffing
            elseif (iscell(valuecol{i}{j}))
                textmsg{i}{nline(i)} = valuecol{i}{j}{k};
            else
                textmsg{i}{nline(i)} = valuecol{i}{j};
            end
        end
    end
end

function outputevents(pinfo, data, spikeselection, hfield)
%%%output events identified by the selected time-related variables
ok = 0; value1 = []; value2 = []; valueR = []; nV = 0; %%%start, end, and ref values
allfinaldate = pinfo.general.finaldir(spikeselection);
ufinal = unique(allfinaldate); %%%this contains all final directories to write to: .\final\events\.evt files)
%%%%work out one paired variables to write events
fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle);
for (i = 1:nfield) subfield{i} = fieldnames(pinfo.(fieldtitle{i})); end
for (i = 1:nfield) %search selected sub-fields in each big field
    fieldselection = getappdata(hfield(i), 'selection');
    fieldselectindex = find( fieldselection == 1); ifreverse = 0;
    if numel(fieldselectindex) == 2 
        Fname = fieldtitle{i}; field1 = subfield{i}{fieldselectindex(1)}; field2 = subfield{i}{fieldselectindex(2)}; field3 = [];
        sel = questdlg(['Two variables selected are: ', field1, ' and ', field2, '. Reverse?']);
        if strcmp(sel, 'No')
            value1 = pinfo.(Fname).(field1); value2 = pinfo.(Fname).(field2); 
            nV = 2; ok = 1; break
        elseif strcmp(sel, 'Yes')
            value2 = pinfo.(Fname).(field1); value1 = pinfo.(Fname).(field2); ifreverse = 1;
            nV = 2; ok = 1; break
        end
    elseif numel(fieldselectindex) == 1
        Fname = fieldtitle{i}; field1 = subfield{i}{fieldselectindex(1)}; field2 = []; field3 = [];
        sel = questdlg(['One variable selected: ', field1, '. Continue?']);
        if strcmp(sel, 'Yes')
            value1 = pinfo.(Fname).(field1); nV = 1; ok = 1; break
        end
    elseif numel(fieldselectindex) == 3
        Fname = fieldtitle{i}; field1 = subfield{i}{fieldselectindex(1)}; 
        field2 = subfield{i}{fieldselectindex(2)}; field3 = subfield{i}{fieldselectindex(3)};
        fnn = {field1; field2; field3};
        [rsel, ok] = listdlg('ListString', fnn, 'PromptString', 'Which is the reference?');
        if ok 
            ii = setdiff([1:3], rsel);   
            value1 = pinfo.(Fname).(fnn{ii(1)}); value2 = pinfo.(Fname).(fnn{ii(2)}); valueR = pinfo.(Fname).(fnn{rsel});
            nV = 3; ok = 1; break
        end
    end
end
if ok
    qq = questdlg('Filter variables through particular events?');
    if strncmpi(qq, 'Yes', 1)
       input = inputdlg({'Event keyword'; 'Event type'; 'Event keyNOword'; 'Event NOtype'}, 'Event selection', 4, {'track'; 'run'; 'first'; ''}); 
       if (~isempty(input))
          evkeyword = input{1}; evkeytype = input{2}; evkeynoword = input{3}; evkeynotype = input{4};
       else
          ok = 0;
       end
    elseif strncmpi(qq, 'No', 1)
        evkeyword = []; evkeytype = [];
    else
        ok = 0;
    end
end
if ok
    disp(['-------> catogery name: ', Fname, '; fields selected: 1-', field1, '; 2-', field2, '; 3-', field3]);
    qq = questdlg(['Replace fieldname ' Fname ' by a new name tag?']);
    if strncmpi(qq, 'Yes', 1)
       input = inputdlg({'New name tag:'}, 'Event selection', 1, {'HeadTwitch'}); 
       if (~isempty(input))
           Fname = input{1}; 
       else
          ok = 0;
       end
    elseif strncmpi(qq, 'No', 1)
        ok = 1;
    else
        ok = 0;
    end
end
if ok   
    for (i = 1:numel(ufinal))
        disp(['-------> final dir now: ', ufinal{i}]);
        evT1 = []; evT2 = []; evTR = []; allsess = []; sessST = []; sessET = []; evTimes = [];
        spikeind = spikeselection(find(strcmp(allfinaldate, ufinal{i}))); %%%selected items belong to the final directory
        for (j = 1:numel(spikeind))
            if (~isempty(evkeyword)) || (~isempty(evkeytype)) || (~isempty(evkeynoword)) || (~isempty(evkeynotype))
                [~, ~, evTimes{j}] = filterevents(pinfo, data, evkeyword, evkeytype, evkeynoword, evkeynotype, spikeind(j));
            end
            if nV==2
               V1 = value1{spikeind(j)}; V2 = value2{spikeind(j)}; 
               V1 = reshape(V1, numel(V1), 1); V2 = reshape(V2, numel(V2), 1); VR = zeros(numel(V1), 1);
               %[a, b] = size(V1); if (a==1) V1 = V1'; end %%[a, b] = size(V2); if (a==1) V2 = V2'; end
            elseif nV ==1
               V1 = value1{spikeind(j)}; V1 = reshape(V1, numel(V1), 1); V2 = V1 + ones(numel(V1), 1); VR = zeros(numel(V1), 1);
            elseif nV==3
               V1 = value1{spikeind(j)}; V2 = value2{spikeind(j)}; VR = valueR{spikeind(j)};
               V1 = reshape(V1, numel(V1), 1); V2 = reshape(V2, numel(V2), 1); VR = reshape(VR, numel(V2), 1);
            end
            evT1 = [evT1; V1]; evT2 =[evT2; V2]; evTR = [evTR; VR];
            sessname = pinfo.general.sessname{spikeind(j)};
            if iscell(sessname)
                for (k = 1:numel(sessname))
                    if isempty(find(strcmp(allsess, sessname{k})))
                        allsess{numel(allsess)+1} = sessname{k}; 
                        sessST(numel(sessST)+1) = pinfo.general.sessstartT{spikeind(j)}(k);
                        sessET(numel(sessET)+1) = pinfo.general.sessendT{spikeind(j)}(k);
                    end
                end
            else
                if isempty(find(strcmp(allsess, sessname)))
                        allsess{numel(allsess)+1} = sessname; 
                        sessST(numel(sessST)+1) = pinfo.general.sessstartT{spikeind(j)}(1);
                        sessET(numel(sessET)+1) = pinfo.general.sessendT{spikeind(j)}(1);
                end
            end
        end
        for (j = 1:numel(allsess)) %%%for each identified session, write an event file
            if (numel(evT1) ~= numel(evT2))||(numel(evT1) ~= numel(evTR))||(numel(evT2) ~= numel(evTR))
                disp(['-----------> warning: identified event start/end/ref times do not match; no event generated for ', ufinal{i}, ': ', allsess{j}]);
            else
                iii = find ((evT1>=sessST(j)) & (evT1<=sessET(j)) & (evT2>=sessST(j)) & (evT2<=sessET(j)));
                ev1now = evT1(iii); ev2now = evT2(iii); evRnow = evTR(iii); 
                %ev1now = sort(unique(ev1now)); ev2now = sort(unique(ev2now)); evRnow = sort(unique(evRnow)); %%%unique and sort
                [ev1now, iii] = unique(ev1now); ev2now = ev2now(iii); evRnow = evRnow(iii);   %%%unique and sort according to ev1 since ev2 and ref matched durign generation
                [ev1now, iii] = sort(ev1now); ev2now = ev2now(iii); evRnow = evRnow(iii); 
                if ifreverse %%%if reverse two variables
                   ev1now = [sessST(j); ev1now]; ev2now = [ev2now; sessET(j)]; evRnow = [0; evRnow]; %%%this is possible only if 2 fields are selected 
                end
                if (~isempty(evkeyword)) || (~isempty(evkeytype))
                   [ev1now, ev2now, evRnow] = filtereventtimesnow(ev1now, ev2now, evRnow, evTimes, sessST(j), sessET(j));
                end
                if numel(ev1now)~=0 %%%&& (numel(ev2now)~=0) 
                    %if (numel(ev2now)==numel(ev1now)) 
                        D = ev1now-ev2now;
                        if (numel(find(D<=0)) == numel(ev1now))
                            writeevt(ufinal{i}, allsess{j}, Fname, field1, field2, field3, ev1now, ev2now, evRnow, evkeyword, evkeytype);
                        elseif (numel(find(D>=0)) == numel(ev1now))
                            writeevt(ufinal{i}, allsess{j}, Fname, field2, field1, field3, ev2now, ev1now, evRnow, evkeyword, evkeytype);
                        else
                            disp(['-----------> warning: identified event start/end times disordered; no event generated for ', ufinal{i}, ': ', allsess{j}]);
                        end
                else
                        disp(['-----------> no events found in ', ufinal{i}, ': ', allsess{j}]);
                    %end
                end
            end
        end
    end
end
function writeevt(ufinal, sessname, Fname, field1, field2, field3, ev1now, ev2now, evRnow, evkeyword, evkeytype)
if isempty(evkeyword) && isempty(evkeytype)
   filename = fullfile(ufinal, 'events', strcat(sessname, '_', Fname,'_',field1, field2, field3, '.evt'));
else
   filename = fullfile(ufinal, 'events', strcat(sessname, '_', Fname,'_',field1, field2, field3, '_', evkeyword, evkeytype, '.evt')); 
end
disp(['-----------> writing events to: ', filename]);
ev.start = ev1now; ev.ent = ev2now; ev.ref = evRnow; 
for (i = 1:numel(ev1now))
    ev.marker{i} = 'ND';
end
WriteEventFile(ev, filename);

function [evName, evType, evT] = filterevents(pinfo, data, evkeyword, evkeytype, evkeynoword, evkeynotype, neuronid)
evName = pinfo.general.eventname{neuronid}; evsel = ones(1,numel(evName));
if isfield(pinfo.parm, 'eventtype')
   evType = pinfo.parm.eventtype{neuronid};
else
   evType = pinfo.parm.eventType{neuronid};
end
if isfield(data, 'events')
   evT = data.events.eventtimes{neuronid}; 
else
   evT = data.event.eventtimes{neuronid}; 
end
evsel = checkifcorrectevents(evName, evType, evkeyword, evkeytype, evkeynoword, evkeynotype);
evpos = find(evsel == 1);
evName = evName(evpos); evType = evType(evpos); evT = evT(evpos);
function evsel = checkifcorrectevents(evName, evType, evkeyword, evkeytype, evkeynoword, evkeynotype)
evsel = ones(1,numel(evName));
for (i = 1:numel(evName))
    if (~isempty(evkeytype))||(~isempty(evkeyword)) %%%for inclusion
      if isempty(evkeytype)
         if ~contains(lower(evName{i}), lower(evkeyword)) evsel(i) = 0; end 
      elseif isempty(evkeyword)
         if ~strncmpi(evType{i}, evkeytype, 3) evsel(i) = 0; end 
      else
         if ~(contains(lower(evName{i}), lower(evkeyword)) && strncmpi(evType{i}, evkeytype, 3))
             evsel(i) = 0; 
         end 
      end
    end
    if (~isempty(evkeynotype))||(~isempty(evkeynoword)) %%%for exclusion
       if isempty(evkeynotype)
          if contains(lower(evName{i}), lower(evkeynoword)) evsel(i) = 0; end 
       elseif isempty(evkeynoword)
          if strncmpi(evType{i}, evkeynotype, 3) evsel(i) = 0; end
       else
          if contains(lower(evName{i}), lower(evkeynoword)) || strncmpi(evType{i}, evkeynotype, 3)
              evsel(i) = 0; 
          end   
       end
    end
end

function [ev1now, ev2now, evRnow] = filtereventtimesnow(ev1now, ev2now, evRnow, evTimes, sessST, sessET)%%%%This is too stringent, need to find overlaps
timeunit = 0.001; %time resolution
% iii = []; %%%ev1now, ev2now need to be sorted and match dimensions
% for (i = 1:numel(evTimes))
%     sT = evTimes{i}.start; eT = evTimes{i}.ent;
%     for (j = 1:numel(sT))
%         inow = find( ((ev1now>=sT(j)) & (ev1now<=eT(j))) &  ((ev2now>=sT(j)) & (ev2now<=eT(j))) );
%         iii = union(iii, inow);
%     end
% end
% ev1now = ev1now(iii); ev2now = ev2now(iii); evRnow = evRnow(iii);
%%%%% Finding everlap between [ev1now ev2now] and evTimes.start/ent is an infinitely complex problem
%%% Here use a comprehensive counting approach: judge point by point
evTstart = []; evTent = [];
for i = 1:numel(evTimes) %%%consolidate all filter through events in the SAME FINALDIR, COULD CONTAIN CROSS SESSIONS
    for j = 1:numel(evTimes{i})
        if ~isempty(evTimes{i}{j}.start)
           evTstart = [evTstart; evTimes{i}{j}.start]; evTent = [evTent; evTimes{i}{j}.ent];
        end
    end
end
%%%%filter evTimes through the session
iii = find( ((evTstart>=sessST) & (evTstart<=sessET)) &  ((evTent>=sessST) & (evTent<=sessET)) );
evTstart=evTstart(iii); evTent = evTent(iii);
[evTstart, iii] = sort(evTstart); evTent = evTent(iii);
%%%%now find the overlap between ev1/2now and evTstart/ent
nev = numel(ev1now); nT = numel(evTstart); 
%disp(['filter start end time ' num2str(evTstart(1)) '* end time ' num2str(evTent(numel(evTent)))]);
if (~isempty(ev1now)) && (~isempty(evTstart))
   st = ev1now(1); et = ev2now(nev); allT = st:timeunit:et; 
   s1=findpoints(allT, ev1now, ev2now); s2=findpoints(allT, evTstart, evTent);
   sel = s1&s2; [sind, eind] = parseindex(sel);
   ev1now = allT(sind); ev2now = allT(eind); evRnow = findref(ev1now, ev2now, evRnow);
else
   ev1now = []; ev2now = []; evRnow = [];
end
function sel = findpoints(allT, sT, eT)
sel = zeros(1,numel(allT)); iii = []; 
for (j = 1:numel(sT))
     inow = find( ((allT>=sT(j)) & (allT<=eT(j))) );
     iii = union(iii, inow);
end
sel(iii) = ones(1,numel(iii));
function R = findref(ev1now, ev2now, evRnow)
R=zeros(size(ev1now));
for i=1:numel(ev1now)
    for j=1:numel(evRnow)
        if (evRnow(j)>=ev1now(i)) && (evRnow(j)<=ev2now(i))
            R(i) = evRnow(j); break
        end
    end
end
function [sind, eind] = parseindex(sel)
sind = []; eind = [];
if numel(sel)>1
   if isempty(find(sel==0))
       sind = 1; eind = numel(sel);
   else
       D = diff(sel); sind = (find(D==1))+1; eind = find(D==-1);
       if sel(1)==1 sind = [1 sind]; end
       if sel(numel(sel))==1 eind = [eind numel(sel)]; end
   end
end
if isempty(sind) || isempty(eind)
    disp('-----------> warning: event starts or ends not found; no events generated')
    sind = []; eind = [];    
elseif (numel(sind) ~= numel(eind)) || (~isempty(find(sind-eind>0)))
    disp('-----------> warning: event starts and ends disordered; no events generated')
    sind = []; eind = [];
end
% if sel(1)==0 
%     if sel(numel(sel))==0
%        sind = (find(D==-1))+1; eind = find(D==1);
%     else
%        sind = (find(D==-1))+1; eind = find(D==1); eind = [eind numel(sel)];
%     end
% else
%     if sel(numel(sel))==0
%        sind = (find(D==-1))+1; sind = [1 sind];
%        eind = find(D==1); 
%     else
%        sind = (find(D==-1))+1;  sind = [1 sind];
%        eind = find(D==1); eind = [eind numel(sel)];
%     end
% end

