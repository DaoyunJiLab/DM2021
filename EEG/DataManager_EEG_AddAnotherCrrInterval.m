function DataManager_EEG_AddAnotherCrrInterval
%%%%%%%%%%% add another interval for crr integration
hf = gcbf; hgroup = getappdata(hf, 'hgroup'); ok = 1;
plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
fname = get(hf, 'Name'); ftt = strfind(fname, '__'); currentfilename = fname(ftt+2:numel(fname));
[pp, nn, ee] = fileparts(currentfilename);
if ~strcmp(ee, '.crrdb')
    disp('-----> not a crrdb; aborted'); ok = 0;
else
    pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); plotparm = getappdata(hf, 'plotparm');
    hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); 
    grpind = find(groupselection == 1); ngroup = numel(grpind); cellind = [];
    for (i = 1:ngroup) cellind = union(cellind, data.grouplist.groupindex{grpind(i)}); end 
    input = inputdlg({'Interval start'; 'Interval end'}, 'Interval parameters', 2, {'-0.2'; '0.2'}); 
    if (~isempty(input))
       IntS = str2num(input{1}); IntE = str2num(input{2}); 
    else
       ok = 0;
    end
end
if ok
    if ~isfield(pinfo, 'crr')
        ok = 0;
    else
        allF = fieldnames(pinfo.crr);
        iii = find(strcmp(allF, 'sessIntcrr')); 
        if isempty(iii)
            sessF = 'sessIntcrr';
        else
            nc = zeros(1, numel(iii)); for (i = 1:numel(iii)) nc(i) = numel(allF{iii(i)}); end
            [~, ij] = max(nc); sessF = strcat(allF{iii(ij)}, '2');
        end
        iii = find(strcmp(allF, 'evIntcrr')); 
        if isempty(iii)
            evF = 'evIntcrr';
        else
            nc = zeros(1, numel(iii)); for (i = 1:numel(iii)) nc(i) = numel(allF{iii(i)}); end
            [~, ij] = max(nc); evF = strcat(allF{iii(ij)}, '2');
        end
        ncell = numel(pinfo.general.finaldir); 
        pinfo.crr.(sessF) = cell(1, ncell); pinfo.crr.(evF) = cell(1, ncell);
        allP = fieldnames(pinfo.parm);
        iii = find(strcmp(allP, 'intbin')); 
        if isempty(iii)
            intP = 'intbin';
        else
            nc = zeros(1, numel(iii)); for (i = 1:numel(iii)) nc(i) = numel(allP{iii(i)}); end
            [~, ij] = max(nc); intP = strcat(allP{iii(ij)}, '2');
        end
        for (i = 1:ncell) pinfo.parm.(intP){i} = [IntS IntE]; end
    end
end
if ok
    for (i = 1:numel(cellind))   
         cellnow = cellind(i); %%%this is the original pair index
         nsess = numel(pinfo.general.sessionname{cellnow});
         for (j = 1:nsess)
             pinfo.crr.(sessF){i}{j} = [];
             if ~isempty(data.crr.sesscrr{i}{j})
                 pinfo.crr.(sessF){i}{j} = findnewcrrvalue(data.crr.sesscrr{i}{j}, data.crr.timebin{i}, IntS, IntE);
             end
         end
         nev = numel(pinfo.general.eventname{cellnow});
         for (j = 1:nev)
             pinfo.crr.(evF){i}{j} = [];
             if ~isempty(data.crr.evcrr{i}{j})
                 pinfo.crr.(evF){i}{j} = findnewcrrvalue(data.crr.evcrr{i}{j}, data.crr.timebin{i}, IntS, IntE);
             end
         end
    end
end
if ok
    save(currentfilename, 'pinfo', 'data', '-mat');
end
disp('**********************');

function crr = findnewcrrvalue(evcrr, timebin, S, E)
crr = [];
if (~isempty(evcrr)) && (~isempty(timebin))
    iii = find( (timebin>=S) & (timebin<= E) );
    if ~isempty(iii)
        crr = mean(evcrr(iii));
    end
end




