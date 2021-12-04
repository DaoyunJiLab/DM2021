function DM_Seq_BreakSequenceDB_Callback
%%display templates in DataAnimator menu

[fname, pname] = uigetfile(fullfile(cd, '*.seqdb'), 'Select a sequenceDB file');
tmpfile = fullfile(pname, fname); [pp, nn, ee] = fileparts(tmpfile);
S = load(tmpfile, '-mat'); pinfo = S.pinfo; data = S.data; S = [];
ok = 1;
disp('----> break a sequence database to individual dates');
aniname = pinfo.general.animalname; datedir = pinfo.general.datedir; anidate = cell(1, numel(aniname));
for (i = 1:numel(aniname))
    anidate{i} = strcat(aniname{i}, datedir{i});
end
unidate = unique(anidate); nam = numel(unidate);
if nam > 1
    pstr = [num2str(nam), ' dates found; break to individual dates?'];
    cc = questdlg(pstr); 
    if ~strcmp(cc, 'Yes')
        ok = 0; disp('-----> cancelled');
    end
end
if ok
    for (i = 1:nam)
        datenow = unidate{i};
        nt = find( strcmp(anidate, datenow) );
        if (~isempty(nt))
            tname = fullfile(pname, strcat(nn, '_', datenow, '.seqdb'));
            [pnow, dnow] = selectentriesindb(pinfo, data, nt);
            saveseqnow(tname, pnow, dnow);
        end
    end
end
disp('****************');
        
function saveseqnow(tname, pnow, dnow)
pinfo = pnow; data = dnow; file = tname; ok = 1;
if exist(tname, 'file') == 2
    str = [tname, ' already exist; overwrite?'];
    cc = questdlg(str);
    if strcmp(cc, 'Yes')
        file = tname;
    else
        [fname, pname] = uiputfile(fullfile(cd, '*.seqdb'), 'Save to a seqdb file');
        if fname ~= 0
           file = fullfile(pname, fname);
        else
            ok = 0;
        end
    end
end
save(file, 'pinfo', 'data', '-mat');

function [pnow, dnow] = selectentriesindb(pinfo, data, nt)
cat = fieldnames(pinfo); 
for (i = 1:numel(cat))
    subf = fieldnames(pinfo.(cat{i}));
    for (j = 1:numel(subf))
        pnow.(cat{i}).(subf{j}) = pinfo.(cat{i}).(subf{j})(nt);
    end
end
cat = fieldnames(data); 
for (i = 1:numel(cat))
    if strcmp(cat{i}, 'grouplist')
       dnow.grouplist.groupname{1} = 'List0'; dnow.grouplist.groupindex{1} = 1:numel(nt);
       dnow.grouplist.grouptype{1} = 'Manual'; dnow.grouplist.groupcrit{1} = []; dnow.grouplist.groupparents{1} = [];  
    elseif strcmp(cat{i}, 'parentfile')
       dnow.parentfile = data.parentfile; 
    else
       subf = fieldnames(data.(cat{i}));
       for (j = 1:numel(subf))
           dnow.(cat{i}).(subf{j}) = data.(cat{i}).(subf{j})(nt);
       end
    end
end
