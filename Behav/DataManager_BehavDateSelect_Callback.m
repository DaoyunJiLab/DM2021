function DataManager_BehavDateSelect_Callback
%%select spikes from spike database according to input criteria
%%output selected spikes to a group
%%criteria must contain full field name, relational operater, and values in sequence with spaces in between:
%%    behavior.trajlength >= 10
%%    general.recarea = CA1
%%
hf = gcbf; pinfo = getappdata(hf, 'behav'); data = getappdata(hf, 'bhdata');
hspike = getappdata(hf, 'hspike'); hgroup= getappdata(hf, 'hgroup'); grouplistpos= getappdata(hf, 'grouplistpos');
spikeselection = getappdata(hspike, 'selection'); groupselection = getappdata(hgroup, 'selection');
nspike = numel(spikeselection); %number of spike (name)s
htext = getappdata(hspike, 'htext'); %every line (spike name) is a text object
ntext = numel(htext); %number of text lines displayed: not all spikes being displayed

tagmark = get(gcbo, 'Tag');
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
elseif (strcmp(tagmark, 'savegroup')) %to show the group spike list
   groupsel = find(groupselection == 1); hfield = getappdata(hf, 'hfield');
   savegrouplist(pinfo,data,groupsel, hfield);
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
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, 'bhdata', data);
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
   groupsel = find(groupselection == 1); ind = [];  ntype = 'Combine(';
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
       disp('-------> no sessions after the combination');
   else
       input = inputdlg({'Group name'}, 'Enter group name', 1);
       if (~isempty(input{1}))
           ngroup = numel(data.grouplist.groupname);
           data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = ind;
           data.grouplist.grouptype{ngroup+1} = ntype; data.grouplist.groupcrit{ngroup+1} = [];
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, 'bhdata', data);
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
elseif (strcmp(tagmark, 'subtractgroup')) %to subtract a small one from a large group
   groupsel = find(groupselection == 1); ind = [];
   if (numel(groupsel)~=2)
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
           disp('-------> no sessions after the subtraction');
       else
           input = inputdlg({'Group name'}, 'Enter group name', 1);
           if (~isempty(input{1}))
           ngroup = numel(data.grouplist.groupname);
           data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = ind;
           ntype = strcat('Subtract(', nam1, '_', nam2, ')');
           data.grouplist.grouptype{ngroup+1} = ntype; data.grouplist.groupcrit{ngroup+1} = [];
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, 'bhdata', data);
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
           disp('-------> no sessions after the reverse intersection');
       else
           input = inputdlg({'Group name'}, 'Enter group name', 1);
           if (~isempty(input{1}))
           ngroup = numel(data.grouplist.groupname);
           data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = ind;
           data.grouplist.grouptype{ngroup+1} = ntype; 
           data.grouplist.groupcrit{ngroup+1} = [];
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, 'bhdata', data);
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
elseif (strcmp(tagmark, 'reverseintgroup')) %to choose the intersection of 2 or more groups
   groupsel = find(groupselection == 1); ind = []; ntype = 'ReverseInt(';
   if (numel(groupsel)<2)
       disp('-------> less than 2 groups selected');
   else
       indint = data.grouplist.groupindex{groupsel(1)}; indall = data.grouplist.groupindex{groupsel(1)}; 
       ntype = strcat(ntype, data.grouplist.groupname{groupsel(1)}, '_');
       for (ti = 2:numel(groupsel))
           indnow = data.grouplist.groupindex{groupsel(ti)};
           indint = intersect(indint, indnow); indall = union(indall, indnow);
           if (ti == numel(groupsel))
               endstr = ')';
           else
               endstr = '_';
           end
           ntype = strcat(ntype, data.grouplist.groupname{groupsel(ti)}, endstr);
       end
       ind = intersect(indall, indint);
       if (isempty(ind))
           disp('-------> no sessions after the reverse intersection');
       else
           input = inputdlg({'Group name'}, 'Enter group name', 1);
           if (~isempty(input{1}))
           ngroup = numel(data.grouplist.groupname);
           data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = ind;
           data.grouplist.grouptype{ngroup+1} = ntype; 
           data.grouplist.groupcrit{ngroup+1} = [];
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, 'bhdata', data);
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
       if (~isempty(input{1}))
           ngroup = numel(data.grouplist.groupname);
           data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = spikeselectind;
           data.grouplist.grouptype{ngroup+1} = 'Manual'; data.grouplist.groupcrit{ngroup+1} = [];
           %%%%%update the group list field
           delete(hgroup);
           hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList');
           setappdata(hf, 'hgroup', hgroup); setappdata(hf, 'bhdata', data);
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
       disp('-------> no sessions selected')
   end
elseif (strcmp(tagmark, 'autoselect')) %if select spikes by criterion
   %get instruction to select spikes             %if first criterion is 'not': detele all spikes that satisfy all rest of the criteria
   input = inputdlg({'Number of criterion sets'}, 'Selection input', 1);
   ncrt = str2num(input{1});
   for (kkk = 1:ncrt)
       s1 = {'Criterion #1:'; 'Criterion #2:'; 'Criterion #3:'; 'Criterion #4:'; 'Criterion #5:';}; %select all spikes that satisfy all the criteria
       input = inputdlg(s1, 'Criterion input', 5);
       for (i = 1:numel(input))
           if (~isempty(input{i})) 
              crt{kkk}{i} = input{i};
           end
       end
   end
    
 %crt{1}{1} = 'general.date == 9_3_03'; crt{2}{1} = 'general.recarea == CA1';
 %  crt{1}{1} = 'general.recarea == CA*';
   %crt{2}{1} = 'parm.timeunit == 0.0001';
 %  crt{2}{1} = 'firing.runoverallrate < 10';
   %crt{2}{2} = 'firing.runmaxrate > 30';
%    crt{1}{1} = 'or';
%    crt{1}{2} = 'general.date == 8_21_03'; %quick
%    crt{1}{3} = 'general.date == 9_3_03';
%    crt{1}{4} = 'general.date == 9_21_03';
%    crt{1}{5} = 'general.date == 3_23_03'; %mongol
%    crt{1}{6} = 'general.date == 4_11_03';
%    crt{1}{7} = 'general.date == 4_22_03';
%    crt{1}{8} = 'general.date == 12_17_03';%pascalTS
%    crt{1}{9} = 'general.date == 11_30_03';
%    crt{1}{10} = 'general.date == 12_12_03';
%    %crt{1}{11} = 'general.date == 10_26_03';
%    
%    %crt{2}{1} = 'not'
%    %crt{2}{1} = 'general.recarea == V*';
%    crt{2}{1} = 'general.recarea == CA1';
%    %crt{2}{2} = 'firing.runoverallrate >= 0.2';
%    crt{2}{2} = 'firing.runoverallrate < 4';
%    %crt{2}{2} = 'firing.runcsi >= 5';
%    crt{2}{3} = 'waveform.runmaxhlwd >= 14';
%    crt{2}{4} = 'firing.runoverallrate >= 0.2';
%  ncrt = 2;
       groupsel = find(groupselection == 1); cellind = [];
       for (tti = 1:numel(groupsel))
           cellind = union(cellind, data.grouplist.groupindex{groupsel(tti)});
       end
       parmfilenow = cell(1,ncrt);
       for (kkk = 1:ncrt)
           parmfilenow{kkk} = DataManager_SearchSpikeDatabase(pinfo, crt{kkk}, 'general', 'sessID', cellind);
       end
       parmfile = findcommoncell(parmfilenow);
       disp(['-----> ', num2str(numel(parmfile)), ' sessions match the criteria']);
       spikesel = 0*spikeselection;
       if (~isempty(parmfile))
          %mark spikes in pinfo that match the spikelist: 
          for (i = 1:nspike)
              if (~isempty(find(strcmp(parmfile, pinfo.general.sessID{i}))))
                 spikesel(i) = 1;
              end
          end
       end
       spikeindnow = find(spikesel == 1);
       if (~isempty(spikeindnow))
          input = inputdlg({'Group name'}, 'Enter group name', 1);
          if (~isempty(input{1}))
             ngroup = numel(data.grouplist.groupname);
             data.grouplist.groupname{ngroup+1} = input{1}; data.grouplist.groupindex{ngroup+1} = spikeindnow;
             data.grouplist.grouptype{ngroup+1} = 'Auto'; data.grouplist.groupcrit{ngroup+1} = crt;
             %%%%%update the group list field
             delete(hgroup);
             hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList');
             setappdata(hf, 'hgroup', hgroup); setappdata(hf, 'bhdata', data);
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
          disp('-------> no sessions selected')
       end
end
%now re-set line (text) color for each spike (de-)selected
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

function parmfile = findcommoncell(parmfilenow)
%%%find cells common in indexnow{1}, indexnow{2}, ...
nset = numel(parmfilenow); parmfile = []; ncl = 0;
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

function savegrouplist(pinfo,data,groupsel, hfield)
%determine all possible field
fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle);
for (i = 1:nfield) subfield{i} = fieldnames(pinfo.(fieldtitle{i})); end
%determine selected groups and selected fields
numgroup = numel(groupsel);
for (xx = 1:numel(groupsel)) %if any groups selected
    spikeselectindex = data.grouplist.groupindex{groupsel(xx)};
    valuecol = []; valuetitle = []; ncol = 1;
    valuetitle{ncol} = data.grouplist.groupname{groupsel(xx)};
    valuecol{ncol} = pinfo.general.datedir(spikeselectindex); %first column = selected spike names
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
    nvar = numel(textmsg); 
    for (tti = 1:nvar)
        nline = numel(textmsg{tti});
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