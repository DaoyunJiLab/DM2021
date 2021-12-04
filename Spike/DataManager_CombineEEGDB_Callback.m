function DataManager_CombineEEGDB_Callback
%%combine multiple spike database files to a single spike database
%%remove the repeated cells, combine cells
%%%%%%% Attention %%%%%%%%%%%%%
%%%%%Only cell arrays and non-empty numeric matrix assignments are allowed as entries such as:
%%%%%  pinfo.general.celltype{i}  or   pinfo.general.fieldnumber(i)
%%%%%  set the database with most fields/subfields as the first db to combine with others 
%%%%%  data.grouplist.groupname & groupindex are reassigned to default List0

%get multiple .spikedb files
disp('Combine EEG database files');
startpath = cd;
[dbfile, filename, filepath, fileext, ok] = FileDlg(startpath, '.eegdb');  %get a group of database files

nfile = numel(dbfile);
if (nfile <= 1)
    disp('-----> less than 2 eeg database files selected; action aborted');
else
    %read first .spikedb file structure to the final output structure
    disp('-----> adding data structure');
    S = load(dbfile{1}, '-mat'); %load the first file
    pinfo = S.eeg; data = S.eegdata; S = []; data = rmfield(data, 'grouplist'); %load the first structure to the final (pinfo) structure

    %now load other files and check one by one to add to the main structure
    for (i = 2:nfile)
        S = load(dbfile{i}, '-mat');
        nowpinfo = S.eeg; nowdata = S.eegdata; S= []; nowdata = rmfield(nowdata, 'grouplist');
        [pinfo, data] = CheckAndAddDB(pinfo, nowpinfo, data, nowdata);
    end
    %re-create the default grouplist
    data.grouplist.groupname{1} = 'List0'; data.grouplist.groupindex{1} = 1:numel(pinfo.general.eegfile);
    data.grouplist.grouptype{1} = 'Manual'; data.grouplist.groupcrit{1} = [];
    
    disp('-----> writing data back to a file');
    [fname, pname] = uiputfile(fullfile(cd, '*.eegdb'), 'Write combined spike database to:');
    writefilename = fullfile(pname, fname);
    eeg = pinfo; eegdata = data; pinfo = []; data = [];
    save(writefilename, 'eeg', 'eegdata');
end
disp('**********************************');

function [pinfo, data] = CheckAndAddDB(pinfo, nowpinfo, data, nowdata)
%add information in structure nowpinfo to structure pinfo
%if datepath differs, just add
nspike = numel(pinfo.general.animalname); nspikenow = numel(nowpinfo.general.animalname);
for (i = 1:nspike)
    checkstr{i} = pinfo.general.eegfile{i};
end
for (i = 1:nspikenow)
    checkstrnow = nowpinfo.general.eegfile{i};
    if ( isempty(find(strcmp(checkstr, checkstrnow))) ) %if the now str is not included
        nspike = nspike + 1;
        checkstr{nspike} = checkstrnow;
        pinfo = addspike(pinfo, nowpinfo, i, nspike); data = addspike(data, nowdata, i, nspike);
    end
end

function pinfo = addspike(pinfo, nowpinfo, i, k);
%add the i-th spike in nowpinfo to pinfo as k-th spike
%%use all fields/subfields in pinfo, empty if not found in nowpinfo
%%%%%get all the fields in pinfo
fieldlist = fieldnames(pinfo);
nfield = numel(fieldlist);
for (iii = 1:nfield)
    subfieldlist = [];
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
               pinfo.(fieldlist{iii}).(subfieldlist{jjj}){k} = nowpinfo.(fieldlist{iii}).(subfieldlist{jjj}){i};
            else
               pinfo.(fieldlist{iii}).(subfieldlist{jjj})(k) = nowpinfo.(fieldlist{iii}).(subfieldlist{jjj})(i);
            end
        end
    end
end
