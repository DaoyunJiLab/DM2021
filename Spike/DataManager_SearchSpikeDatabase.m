function outfile = DataManager_SearchSpikeDatabase(pinfo, crt, fID, subfID, cellind)
%%select spikes from spike database pinfo according to criteria crt
%%output spike full file names (s1file, runfile, s2file)
%%pinfo.(fID).(subfID) is used to identiy each spike
outfile = []; spikeindex = []; index = [];

%get field name, operator, and value from criteria
ncrt = numel(crt);
if (strcmp(crt{1}, 'not') | (strcmp(crt{1}, 'or')) ) %only find spike index that satisfyin the rest of the criteria
    for (i = 2:ncrt)
        [field{i}, rem] = strtok(crt{i}, ' '); %string before first space
        rem = rem(2:numel(rem));
        [opr{i}, behind] = strtok(rem, ' ');
        value{i} = behind(2:numel(behind));
        %%find spike indices according to the criterion
        spikeindex{i} = findspike(pinfo, field{i}, opr{i}, value{i});
    end
else
    for (i = 1:ncrt)
        [field{i}, rem] = strtok(crt{i}, ' '); %string before first space
        rem = rem(2:numel(rem));
        [opr{i}, behind] = strtok(rem, ' ');
        value{i} = behind(2:numel(behind));
        %%find spike indices according to the criterion
        spikeindex{i} = findspike(pinfo, field{i}, opr{i}, value{i});
    end
end

%select wanted spike index in variable index
if (strcmp(crt{1}, 'not')) %exclude common elements in all spikeindex found
   if (ncrt >= 2)
      commonindex = spikeindex{2};
      for (i = 3:ncrt)
          commonindex = intersect(commonindex, spikeindex{i});
      end
      allindex = 1:numel(pinfo.(fID).(subfID)); %all spike index
      index = setdiff(allindex, commonindex);
   end
elseif (strcmp(crt{1}, 'or'))
   if (ncrt >= 2)
      commonindex = spikeindex{2};
      for (i = 3:ncrt)
          commonindex = union(commonindex, spikeindex{i});
      end
      index = commonindex;
   end
else %%select common elements in all spikeindex found
   if (ncrt >= 1)
      index = spikeindex{1};
      for (i = 2:ncrt)
          index = intersect(index, spikeindex{i});
       end
   end
end  
index = intersect(index, cellind);
if (~isempty(index))
    outfile = pinfo.(fID).(subfID)(index);
end
%numel(index)


function spikeindex = findspike(pinfo, field, opr, value)
%select spike index from single criterion
spikeindex = []; feature = [];
%match input to pinfo fields
[fieldname, rem] = strtok(field, '.');
subfieldname = rem(2:numel(rem));
% fdname = fieldnames(pinfo);
% fdindex = find( strcmp(fdname, fieldname) );
if (isfield(pinfo, fieldname))
%     subfdname = fieldnames(pinfo.(fdname(fdindex(1))));
%     subfdindex = find( strcmp(subfdname, subfieldname) );
    if (isfield(pinfo.(fieldname), subfieldname))
        feature = pinfo.(fieldname).(subfieldname); %if find a field and subfield match, get all values
    else
        disp('-----> warning: no field match');
    end
end

%select index in feature by opr and value
%subfield values could be one of three types:
%%% (1) a vector of single values (a value for each spike)
%%% (2) a vector-like cell array of single strings (a string for each spike)
%%% (3) a vector-like cell array of cell arrays of single strings (a cell array for each spike)
%%% (4) a vector-like cell array of numeric vectors (a vector for each spike)

if (~isempty(feature))
    %disp('----in now');
    ttype = resolvetypeformultiple(feature); %%%feature contains all the data of the subfield
    switch ttype
        case 'vector'; %if not a cell - a vector of single values
            if (isnumeric(value)) value = num2str(value); end
            allowopr = {'>='; '<='; '>'; '<'; '=='; '~='};
            if (~isempty(find( strcmp(allowopr, opr) )))
               spikeindex = eval( ['find(feature', opr, value, ')'] ); %simply find index
            elseif (opr == '=') % = is also allowed
               spikeindex = eval( ['find(feature == ', value, ')'] ); %simply find index
            else
               disp('-----> warning: wrong operator');
            end
        case 'stringcell' %if a cell array of strings: only == and ~= allowed
            if (isnumeric(value)) value = num2str(value); end
            [all, sss, pre, post] = analyse(value);
            for (i = 1:numel(feature))
                tt = matchstr(feature{i}, all, sss, pre, post);
                if (strcmp(opr, '==') | strcmp(opr, '='))
                    if (tt) spikeindex = [spikeindex i]; end
                elseif (strcmp(opr, '~='))
                    if (~tt) spikeindex = [spikeindex i]; end
                end
            end
        case 'vectorcell'
            if (ischar(value)) value = str2num(value); end
            flag = askvectorkeyword(field, 'vector', feature);
            spikeindex = findspikevectorcell(feature, opr, value, flag);
        case 'cellcell'
            if (isnumeric(value)) value = num2str(value); end
            [all, sss, pre, post] = analyse(value);
            flag = askvectorkeyword(field, 'cell', feature);
            spikeindex = findspikecellcell(feature, opr, all, sss, pre, post, flag);
        case 'numericstrcell'
            for (k = 1:numel(feature))
                if (~isempty(feature{k}))
                   feature{k} = str2num(feature{k});
                end
            end
            if (ischar(value)) value = str2num(value); end
            flag = askvectorkeyword(field, 'vector', feature);
            spikeindex = findspikevectorcell(feature, opr, value, flag);
    end
else
    disp('-----> warning: feature is empty');
end

function ttype = resolvetypeformultiple(feature)
%disp(feature{1});
ttype = 'none';
if (isnumeric(feature))
    ttype = 'vector';
elseif (iscellstr(feature))
    ttype = 'stringcell';
elseif (iscell(feature)) %probe the real value in each individual element and assign array type
    for (i = 1:numel(feature))
        if (~isempty(feature{i}))
            if (iscellstr(feature{i}))
                ttype = 'cellcell';
            elseif (isnumeric(feature{i}))
                ttype = 'vectorcell';
            elseif (ischar(feature{i}))
                if (isempty(str2num(feature{i}))) %if not a numeric string
                   ttype = 'stringcell';
                else
                   ttype = 'numericstrcell';
                end
            elseif iscell(feature{i})
                for (j = 1:numel(feature{i}))
                    if (~isempty(feature{i}{j}))
                        if (isnumeric(feature{i}{j}))
                            ttype = 'vectorcell';
                        else
                            ttype = 'cellcell';
                        end
                        break
                    end
                end
            end
            break
        end
    end
end

function [all, sss, pre, post] = analyse(value)
%%analyze value string and find components relative to wild character '*'
sss = []; pre = []; post = []; all = [];
if (strcmp(value, '*'))
    all = 1; %if value == '*', all = 1
elseif (~isempty(value))
    all = 0;
    ind = find(value == '*');
    if (isempty(ind))
        sss = value;
    elseif (ind == 1)
        post = value(2:numel(value));
    elseif (ind == numel(value))
        pre = value(1:numel(value)-1);
    else
        pre = value(1:ind-1);
        post = value(ind+1:numel(value));
    end
end

function tt = matchstr(string, all, sss, pre, post)
tt = 0;
%%%%first make everything lower case: the match here is not case-sensitive
string = lower(string); all = lower(all); sss = lower(sss); pre = lower(pre); post = lower(post);
if (all)
    tt = 1;
else
    if (~isempty(sss))
        if ~isempty(strfind(string, sss)) tt = 1; end
    elseif (~isempty(pre)) & (~isempty(post))
        if (~isempty(strfind(string, pre))) & (~isempty(strfind(string, post)))  tt = 1; end
    elseif (~isempty(pre)) & (isempty(post))
        if (~isempty(strfind(string, pre))) tt = 1; end
    elseif (~isempty(post)) & (isempty(pre))
        if (~isempty(strfind(string, post))) tt = 1; end    
    end
end

% function iftrue = matchtrue(string, opr, all, sss, pre, post);
% %%match string to value (all, sss, pre, post) accorind to operator opr
% %%     =, == means at least one of the elements in the cell same as value
% %%     ~= means all elements in the cell different from value
% iftrue = 0;
% if (~isempty(all))  %if value non-empty
%    if ( (strcmp(opr, '==')) | (strcmp(opr, '=')) )
%        if (all == 1)
%            iftrue = 1;
%        elseif (~isempty(sss))
%            if (strcmp(string, sss)) iftrue = 1; end
%        elseif ( (isempty(pre)) & (~isempty(post)) ) %if only post match
%            npost = numel(post); nstring = numel(string);
%            if ( (nstring>=npost) )
%                if (strcmp( string(nstring-npost+1:nstring), post )) iftrue = 1; end
%            end
%        elseif ( (~isempty(pre)) & (isempty(post)) ) %if only pre match
%            npre = numel(pre);
%            if (strncmp(string, pre, npre)) iftrue = 1; end
%        elseif ( (~isempty(pre)) & (~isempty(post)) ) %if both pre and post match
%            npre = numel(pre); npost = numel(post); nstring = numel(string);
%            if (nstring >= npre + npost)
%                if ( (strncmp(string, pre, npre)) & (strcmp( string(nstring-npost+1:nstring), post )) ) iftrue = 1; end
%            end
%        end   
%    elseif (strcmp(opr, '~=')) 
%        if (isempty(string))
%            iftrue = 1;
%        elseif (~isempty(sss))
%            if (~strcmp(string, sss)) iftrue = 1; end
%        elseif ( (isempty(pre)) & (~isempty(post)) ) %if only post match
%            npost = numel(post); nstring = numel(string);
%            if ( (nstring>=npost) )
%                if (~strcmp( string(nstring-npost+1:nstring), post )) iftrue = 1; end
%            end
%        elseif ( (~isempty(pre)) & (isempty(post)) ) %if only pre match
%            npre = numel(pre);
%            if (~strncmp(string, pre, npre)) iftrue = 1; end
%        elseif ( (~isempty(pre)) & (~isempty(post)) ) %if both pre and post match
%            npre = numel(pre); npost = numel(post); nstring = numel(string);
%            if (nstring >= npre + npost)
%                if ( (~strncmp(string, pre, npre)) | (~strcmp( string(nstring-npost+1:nstring), post )) ) iftrue = 1; end
%            end
%        end   
%    else
%     disp('-----> warning: wrong operator');   
%    end
% end

function spikeindex = findspikecellcell(feature, opr, all, sss, pre, post, flag)
%%select spikes from a cell cell feature
%%operator much be =, ==, ~=
%%value must be single string, a single wild mark * allowed in value, but with some other characters
%%     =, == means at least one of the elements in the cell same as value
%%     ~= means all elements in the cell different from value
spikeindex = [];
if ( (strcmp(opr, '==')) | (strcmp(opr, '=')) )
    for (i = 1:numel(feature))
        iftrue = [];
        for (kk = 1:numel(feature{i}))
            iftrue(kk) = matchstr(feature{i}{kk}, all, sss, pre, post);
        end
        nt = numel(find(iftrue));
        if (strcmp(flag, 'all'))
           if (nt==numel(iftrue)) spikeindex = [spikeindex i]; end
        else
           if (nt>0) spikeindex = [spikeindex i]; end
        end
    end
elseif (strcmp(opr, '~='))
    for (i = 1:numel(feature))
        iftrue = [];
        for (kk = 1:numel(feature{i}))
            iftrue(kk) = matchstr(feature{i}{kk}, all, sss, pre, post);
        end
        nt = numel(find(iftrue));
        if (strcmp(flag, 'all'))
           if (nt==0) spikeindex = [spikeindex i]; end
        else
           if (nt<numel(iftrue)) spikeindex = [spikeindex i]; end
        end
    end
else
    disp('-----> warning: wrong operator');
end

function spikeindex = findspikevectorcell(feature, opr, value, flag)
%%select spikes from a vector cell feature
%%operator much be >, >=, <, <=, ==, =, ~=
%%value must be numeric
%%default interpretation of opr:
%%     >, >= means at least one of the elements in a vector bigger than value
%%     <, <= means at least one of the elements in a vector smaller than value
%%     =, == means at least one of the elements in a vector equal to value
%%     ~= means none of elecments in a vector equal tp value
%%flagged interpretation of opr
%%     min: minimum value of a vector
%%     max: maximum value of a vector
%%     mean: mean value of a vector
%%     median: median value of a vector
%%     none: use default interpretation
%%%%%%%%%need to pretreat feature{i}: either a numeric cell or a vector
for (i = 1:numel(feature))
  if (iscell(feature{i}))
    tt = []; n = 0;
    for (j = 1:numel(feature{i}))
        if (~isempty(feature{i}{j}))
            n = n + 1; tt(n) = feature{i}{j};
        end
    end
    feature{i} = tt;
  else
    feature{i} = feature{i}(~isnan(feature{i}));
  end
end

spikeindex = [];
switch opr
case '>'
    for (i = 1:numel(feature))
        if (~isempty(feature{i}))
           if (strcmp(flag, 'default') | strcmp(flag, 'one'))
               valnow = max(feature{i});
           elseif (strcmp(flag, 'all'))
               valnow = min(feature{i});
           else
               valnow = interpret(feature{i}, flag);
           end
           if (valnow > value) spikeindex = [spikeindex i]; end
        end
    end
case '>='
    for (i = 1:numel(feature))
        if (~isempty(feature{i}))
           if (strcmp(flag, 'default') | strcmp(flag, 'one'))
               valnow = max(feature{i});
           elseif (strcmp(flag, 'all'))
               valnow = min(feature{i});
           else
               valnow = interpret(feature{i}, flag);
           end
           if (valnow >= value) spikeindex = [spikeindex i]; end
        end
    end
case '<'
    for (i = 1:numel(feature))
        if (~isempty(feature{i}))
           if (strcmp(flag, 'default') | strcmp(flag, 'one'))
               valnow = min(feature{i});
           elseif (strcmp(flag, 'all'))
               valnow = max(feature{i});
           else
               valnow = interpret(feature{i}, flag);
           end
           if (valnow < value) spikeindex = [spikeindex i]; end
        end
    end
case '<='
    for (i = 1:numel(feature))
        if (~isempty(feature{i})) 
           if (strcmp(flag, 'default') | strcmp(flag, 'one'))
               valnow = min(feature{i});
           elseif (strcmp(flag, 'all'))
               valnow = max(feature{i});
           else
               valnow = interpret(feature{i}, flag);
           end
           if (valnow <= value) spikeindex = [spikeindex i]; end
        end
    end
case '=='
    for (i = 1:numel(feature))
        if (~isempty(feature{i}))
           nv = numel(find(feature{i} == value));
           if (strcmp(flag, 'all'))
              if (nv==numel(feature{i})) spikeindex = [spikeindex i]; end
           else
              if (nv>=1) spikeindex = [spikeindex i]; end
           end 
        end
    end
case '='
    for (i = 1:numel(feature))
        if (~isempty(feature{i}))
           nv = numel(find(feature{i} == value));
           if (strcmp(flag, 'all'))
              if (nv==numel(feature{i})) spikeindex = [spikeindex i]; end
           else
              if (nv>=1) spikeindex = [spikeindex i]; end
           end 
        end
    end
case '~='
    for (i = 1:numel(feature))
        if (~isempty(feature{i}))
           nv = numel(find(feature{i} == value));
           if (strcmp(flag, 'all'))
              if (nv==0) spikeindex = [spikeindex i]; end
           else
              if (nv<numel(feature{i})) spikeindex = [spikeindex i]; end
           end 
        end
    end
otherwise
    disp('-----> warning: wrong operator');
end

function valnow = interpret(feature, flag)
valnow = [];
switch flag
    case 'min'
        valnow = min(feature);
    case 'max'
        valnow = max(feature);
    case 'mean';
        valnow = mean(feature);
    case 'median'
        valnow = median(feature);
end

function flag = askvectorkeyword(field, type, feature)
flag = 'default'; str = []; ok = 0;
for (i = 1:numel(feature))
    if (numel(feature{i})>1) 
        ok = 1; break
    end
end

if ok
prompstr = ['Multiple values assigned to the variable: ', field];
if (strcmp(type, 'vector'))
    str = {'min'; 'max'; 'mean'; 'median'; 'all'; 'one'; 'default'};
elseif (strcmp(type, 'cell'))
    str = {'all'; 'one'; 'default'};
end
if (~isempty(str))
    [sel, ok] = listdlg('ListString', str, 'PromptString', prompstr, 'SelectionMode', 'single', 'InitialValue', numel(str), 'ListSize', [600 300]);
    if (ok)
        flag = str{sel};
    end
end
end



