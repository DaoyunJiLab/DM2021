function DataManager_CombineDatabase_Callback
%%combine multiple spike database files to a single spike database
%%remove the repeated cells, combine cells
%%%%%%% Attention %%%%%%%%%%%%%
%%%%%Only cell arrays and non-empty numeric matrix assignments are allowed as entries such as:
%%%%%  pinfo.general.celltype{i}  or   pinfo.general.fieldnumber(i)
%%%%%  set the database with most fields/subfields as the first db to combine with others 
%%%%%  data.grouplist.groupname & groupindex are reassigned to default List0

%get multiple .spikedb files
disp('Combine database files');
startpath = cd;
[dbfile, filename, filepath, fileext, ok] = FileDlg(startpath, '*.*db*');  %get a group of database files

nfile = numel(dbfile);
if (nfile <= 1)
    disp('-----> less than 2 database files selected; action aborted');
else
    %read first .spikedb file structure to the final output structure
    ok = 1;
    for (i = 2:nfile)
        if (~strcmp(fileext{i}, fileext{i-1})) 
            ok = 0; disp('--------------> file types do not match; combination not performed'); 
            break
        end
    end
    if ok
        disp('-----> adding data structure');
        [pinfo, data] = readdatabase(dbfile{1}); data = rmfield(data, 'grouplist');
        [checkcat, checkvar] = findcheckvar(fileext{1}, pinfo);
        for (i = 2:nfile)
            [nowpinfo, nowdata] = readdatabase(dbfile{i}); nowdata = rmfield(nowdata, 'grouplist');
            [pinfo, data] = CheckAndAddDB(pinfo, nowpinfo, data, nowdata, checkcat, checkvar);
        end
        %re-create the default grouplist
        data.grouplist.groupname{1} = 'List0'; data.grouplist.groupindex{1} = 1:numel(pinfo.general.datedir);
        data.grouplist.grouptype{1} = 'Manual'; data.grouplist.groupcrit{1} = [];
    
        disp('-----> writing data back to a file');
        [fname, pname] = uiputfile(fullfile(cd, fileext{1}), 'Write combined database to:');
        writefilename = fullfile(pname, fname);
        savedatabase(writefilename, pinfo, data, dbfile{1});
    end
end
disp('**********************************');

function [pinfo, data] = CheckAndAddDB(pinfo, nowpinfo, data, nowdata, checkcat, checkvar)
%add information in structure nowpinfo to structure pinfo
%if datepath differs, just add
nspike = numel(pinfo.(checkcat).(checkvar)); nspikenow = numel(nowpinfo.(checkcat).(checkvar));
for (i = 1:nspike)
    checkstr{i} = pinfo.(checkcat).(checkvar){i};
end
for (i = 1:nspikenow)
    checkstrnow = nowpinfo.(checkcat).(checkvar){i};
    if ( isempty(find(strcmp(checkstr, checkstrnow))) ) %if the now str is not included
        nspike = nspike + 1;
        checkstr{nspike} = checkstrnow;
        pinfo = addspike(pinfo, nowpinfo, i, nspike); data = addspike(data, nowdata, i, nspike);
    end
end

function pinfo = addspike(pinfo, nowpinfo, i, k)
%add the i-th spike in nowpinfo to pinfo as k-th spike
%%use all fields/subfields in pinfo, empty if not found in nowpinfo
%%%%%get all the fields in pinfo
fieldlist = fieldnames(pinfo);
nfield = numel(fieldlist);
for (iii = 1:nfield)
    if isstruct(pinfo.(fieldlist{iii}))
       subfieldlist = fieldnames(pinfo.(fieldlist{iii})); %%all the subfields
       for (jjj = 1:numel(subfieldlist))
        if (~isfield(nowpinfo, fieldlist{iii})) %%if field unavailable
            if (iscell(pinfo.(fieldlist{iii}).(subfieldlist{jjj})))
                pinfo.(fieldlist{iii}).(subfieldlist{jjj}){k} = [];
            else
                pinfo.(fieldlist{iii}).(subfieldlist{jjj})(k) = NaN;
            end
        elseif (~isfield(nowpinfo.(fieldlist{iii}), subfieldlist{jjj})) %%if subfield unavailable
            if (iscell(pinfo.(fieldlist{iii}).(subfieldlist{jjj})))
                pinfo.(fieldlist{iii}).(subfieldlist{jjj}){k} = [];
            else
                pinfo.(fieldlist{iii}).(subfieldlist{jjj})(k) = NaN;
            end
        else
            if (iscell(pinfo.(fieldlist{iii}).(subfieldlist{jjj})))
%                disp(i)
%                disp(fieldlist{iii})
%                disp(subfieldlist{jjj})
               pinfo.(fieldlist{iii}).(subfieldlist{jjj}){k} = nowpinfo.(fieldlist{iii}).(subfieldlist{jjj}){i};
            else
               pinfo.(fieldlist{iii}).(subfieldlist{jjj})(k) = nowpinfo.(fieldlist{iii}).(subfieldlist{jjj})(i);
            end
        end
       end
    end
end

function [checkcat, checkvar] = findcheckvar(fileext, pinfo)
checkcat = 'general'; checkvar = 'parmfile';
switch fileext
case '.spikedb'
     checkcat = 'general'; checkvar = 'parmfile';
case '.eegdb'
     checkcat = 'general'; checkvar = 'eegfile';
case '.behavdb'
     checkcat = 'general'; checkvar = 'sessID';
case '.seqdb'
     if contains(pinfo.parm.seqtype{1}, 'evtitemized')
        checkcat = 'general'; checkvar = 'eventname'; 
     else
        checkcat = 'general'; checkvar = 'tmpID';
     end
case '.crrdb'
     checkcat = 'general'; checkvar = 'parmfile';
end

function [pinfo, data] = readdatabase(dbfile)
[pp,nn,ee] = fileparts(dbfile); S = load(dbfile, '-mat');
pstr = 'pinfo'; dstr = 'data'; 
if strcmp(ee, '.eegdb')
     pstr = 'eeg'; dstr = 'eegdata';
elseif strcmp(ee, '.behavdb')
     pstr = 'behav'; dstr = 'bhdata';
end
pinfo = S.(pstr); data = S.(dstr);
S = [];
 
function savedatabase(writefilename, pinfo, data, dbfile)
[pp,nn,ee] = fileparts(dbfile); S = load(dbfile, '-mat');
if strcmp(ee, '.eegdb')
     eeg = pinfo; eegdata = data;
     save(writefilename, 'eeg', 'eegdata', '-mat', '-v7.3');
elseif strcmp(ee, '.behavdb')
     behav = pinfo; bhdata = data;
     save(writefilename, 'behav', 'bhdata', '-mat', '-v7.3'); 
else
     save(writefilename, 'pinfo', 'data', '-mat', '-v7.3');
end 


