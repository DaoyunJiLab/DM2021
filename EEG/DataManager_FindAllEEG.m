function [eeg,eegdata] = DataManager_FindAllEEG(datedir, eeg, varargin)
%Find all eegioral data under the animal directory
%Fields assigned here:
eeg.general.recarea = []; eeg.general.finaldir = []; eeg.general.datedir = []; eeg.general.eegfile = []; eeg.general.eegID = [];
eeg.general.eeggain = []; %%%gain in mV
eeg.general.freq = []; eeg.general.eegTT = []; 
eeg.general.sessname = []; eeg.general.sessstartT = []; eeg.general.sessendT = []; eeg.general.sesslength = [];
eeg.general.eventname = []; eegdata.event.eventtimes = []; %%%this is eegdata
eeg.general.animalname = []; %eeg.general.genotype, sex, age, etc.
nopvar = size(varargin, 2);
copyover = []; curdir = [];
if (nopvar == 1)
    copyover = varargin{1}; curdir = cd;
elseif (nopvar == 2)
    copyover = varargin{1}; curdir = varargin{2};
end
%%%eeg.parm.band, 
%%%        will be assigned in FindeegParm
% 
% finddddir = dir(animalpath);
% neegdir = 0; datedir = [];
% for (j = 1:numel(finddddir))
% if (finddddir(j).isdir == 1) && (~strcmp(finddddir(j).name, '.')) && (~strcmp(finddddir(j).name, '..'))
%     neegdir = neegdir + 1;
%     datedir{neegdir} = finddddir(j).name;
% end
% end

neeg = 0; neegdir = numel(datedir); %start count of eeg files found
default.animname = []; default.moreF = []; default.moreC = []; %%%% For assigning EMG animalname and properties etc.  
%disp(['-----> number of dates ', num2str(neegdir)]);
for (i = 1:neegdir)
    disp(['---------> date: ', datedir{i}]); %datename = []; [datename, tok] = strtok(datedir{i}, '_'); 
    %finalpath = strcat(animalpath, filesep, datedir{i}, filesep, 'final');
    finalpath = datedir{i};
    if (exist(finalpath, 'dir') == 7) %if final dir exist
        infofile = fullfile(finalpath, 'animal.info'); Ttype = []; Tfield = []; Tcontent = [];
        if (exist(infofile)== 2)
           disp(['------------> reading animal info file: ', infofile]);
           [Ttype, Tfield, Tcontent] = ReadAnimalInfoFile(infofile);
        end
        
        sessfile = fullfile(finalpath, 'session.evt');
        disp(['------------> reading session file: ', sessfile]);
        [id, sessions] = ReadEventFile(sessfile);
        
        evdir{1} = strcat(finalpath, filesep, 'events');
        disp(['------------> reading event files in ', evdir{1}]);
        [eventname, eventtimes] = GetAllEvents(evdir); nev = numel(eventname);
        
        eeghisfile = fullfile(finalpath, 'eeg.his');
        disp(['-------> reading EEG histology file: ', eeghisfile]);
        [EEGname, EEGarea] = ReadHisFile(eeghisfile);
        %%%%added to select only 1 EEGname per EEGarea, the others assigned
        %%%%to 'unknown' area
        eArea = unique(EEGarea); eName = cell(size(eArea)); 
        for (ti = 1:numel(eArea))
            iii = find(strcmp(EEGarea, eArea{ti})); eName{ti} = EEGname{iii(1)};
        end
        EEGname = eName; EEGarea = eArea;
        
        disp('------------> checking in all eeg data ');
        [EEGfile, EEGfilearea, eeggain, freq, TTnow, sess] = CheckInAllEEGfiles(finalpath, EEGname, EEGarea, copyover, curdir);  %%check in all eeg file names: no data imported
        for (tk = 1:numel(EEGfile))
             neeg = neeg + 1;
             eeg.general.finaldir{neeg} = finalpath; %{str}
             eeg.general.datedir{neeg} = getdatedir(finalpath); %{str}
             eeg.general.eegfile{neeg} = EEGfile{tk}; %{str}
             [pp, nn, ee] = fileparts(EEGfile{tk}); eeg.general.eegID{neeg} = nn;
             eeg.general.recarea{neeg} = EEGfilearea{tk}; %{str}
             eeg.general.eeggain{neeg} = 1000*eeggain{tk}; %{str} %gain in mV
             eeg.general.freq{neeg} = freq{tk}; %{str}
             eeg.general.eegTT{neeg} = TTnow{tk}; %{str}
             %%%session info
             eeg.general.sessname{neeg} = sess{tk}; eeg.general.sessstartT{neeg} = []; eeg.general.sessendT{neeg} = [];
             eeg.general.sesslength{neeg} = []; eeg.general.eventname{neeg} = []; eegdata.event.eventtimes{neeg} = []; 
             ii = find(strcmpi(sessions.marker, sess{tk}));
             if (numel(ii) == 1)
                eeg.general.sessstartT{neeg} = sessions.start(ii); eeg.general.sessendT{neeg} = sessions.ent(ii);
                eeg.general.sesslength{neeg} = sessions.ent(ii)-sessions.start(ii);
                %%%event info
                ntev = 0; enowname = []; enowtimes = [];
                for (ttk = 1:numel(eventtimes))
                    if ~isempty(eventtimes{ttk}.start)
                       if (min(eventtimes{ttk}.start) >= sessions.start(ii)) && (max(eventtimes{ttk}.ent) <= sessions.ent(ii))
                           ntev = ntev + 1; enowname{ntev} = eventname{ttk}; enowtimes{ntev} = eventtimes{ttk};
                       end
                    end
                end
                eeg.general.eventname{neeg} = enowname; eegdata.event.eventtimes{neeg} = enowtimes; 
             end
             %%%animal info
             [animname, moreF, moreC, default] = findanimalname(Ttype, Tfield, Tcontent, TTnow{tk}, default);
             if (~isempty(animname))
                 eeg.general.animalname{neeg} = animname; 
             else
                 eeg.general.animalname{neeg} = finalpath;    
             end
             for (ikk = 1:numel(moreF)) eeg.general.(moreF{ikk}){neeg} = moreC{ikk}; end
        end
    end
end

function ddd = getdatedir(finalpath)
ddd = []; kk = strfind(finalpath, filesep);
if (numel(kk)>=2)
    ddd = finalpath(kk(numel(kk)-1)+1 : kk(numel(kk))-1);
elseif (numel(kk) == 1)
    ddd = finalpath(1:kk-1);
end

function [animname, moreF, moreC, default] = findanimalname(Ttype, Tfield, Tcontent, TTnow, default)
animname = []; moreF = []; moreC = [];
ii = find( strcmpi(Ttype, 'tetrode') & strcmpi(Tfield, TTnow) );
if (numel(ii) == 1)
    animname = Tcontent{ii};
    kk = find( strcmp(Ttype, animname) );
    moreF = Tfield(kk); moreC = Tcontent(kk);
    default.animname = animname; default.moreF = moreF; default.moreC = moreC;
elseif strcmpi(TTnow, 'EMG')
    animname = default.animname; moreF = default.moreF; moreC = default.moreC;
end

function [eventname, eventtimes] = GetAllEvents(evdir)
eventname = []; eventtimes = []; nev = 0;
[filepath1, filename1] = GetAllFile(evdir, '', '.ep'); nep = numel(filename1);
[filepath2, filename2] = GetAllFile(evdir, '', '.evt'); nevt = numel(filename2);
for (i = 1:nep)
    nev = nev + 1;
    eventname{nev} = filename1{i};
    [id, sessions] = ReadEventFile(fullfile(filepath1{i}, filename1{i}));
    eventtimes{nev} = sessions;
end
for (i = 1:nevt)
    nev = nev + 1;
    eventname{nev} = filename2{i};
    [id, sessions] = ReadEventFile(fullfile(filepath2{i}, filename2{i}));
    eventtimes{nev} = sessions;
end

function [EEGfile, filearea, eeggain, freq, TTnow, sess] = CheckInAllEEGfiles(finalpath, EEGname, EEGarea, copyover, curdir) 
EEGfile = []; filearea = []; eeggain = []; nfile = 0; freq = []; TTnow = []; sess = [];
finddddir = dir(strcat(finalpath, filesep, 'eeg')); neegdir = 0; 
for (j = 1:numel(finddddir))
if (finddddir(j).isdir == 1) && (~strcmp(finddddir(j).name, '.')) && (~strcmp(finddddir(j).name, '..'))
    neegdir = neegdir + 1;
    areatypenow{neegdir} = finddddir(j).name;
end
end
for (i = 1:neegdir)   
    fp{1} = strcat(finalpath, filesep, 'eeg', filesep, areatypenow{i});
    [filepath, filename] = GetAllFile(fp, '', '.eeg'); 
    for(j = 1:numel(filename))
        nfile = nfile + 1; 
        [str, tok] = strtok(filename{j}, '_'); TTnow{nfile} = str; [str, tok] = strtok(tok, '_'); [str, tok] = strtok(tok, '_'); 
        if (~isempty(strfind(str, '.eeg')))
            sess{nfile} = str(1:numel(str)-4);
        else
            sess{nfile} = str;
        end
        EEGfile{nfile} = fullfile(filepath{j}, filename{j}); 
        if strcmp(copyover, 'Yes')
            newfile = constructnewfile(curdir, filepath{j}, filename{j}); 
            copyfile(EEGfile{nfile}, newfile);
            EEGfile{nfile} = newfile; 
        end
        [gain, fr] = readeeggain(EEGfile{nfile}); eeggain{nfile} = gain; freq{nfile} = fr;
        filearea{nfile} = 'unknown';
        for (kkt = 1:numel(EEGname))
            if strcmpi(TTnow{nfile}, EEGname{kkt})
                filearea{nfile} = EEGarea{kkt}; break
            end
        end  
    end
end

function [gain, freq] = readeeggain(filename)
fid = fopen(filename);  %first open a file for read
while 1               %find where header ends
          tline = fgets(fid);
          if contains(tline, 'ADBitVolts') %gain = voits per AD value
             [str, ent] = strtok(tline, ' '); ent = ent(2:numel(ent)); gain = str2num(ent);
          end
          if contains(tline, 'SamplingFrequency') %gain = voits per AD value
             [str, ent] = strtok(tline, ' '); ent = ent(2:numel(ent)); freq = str2num(ent);
          end
%           if (strncmpi(tline(2:numel(tline)), 'SamplingFrequency', 18)) %sampling frequency
%               tlinenow = tline(2:numel(tline)); [str, ent] = strtok(tlinenow, 'cy'); ent = ent(3:numel(ent)); freq = str2num(ent);
%           end
          if (strncmpi(tline, '%%ENDHEADER', 8)) break; end
end
fclose(fid);

function newfile = constructnewfile(curdir, filepath, filename)
%%%add the directory two levels higher than the final directory to the filename
filesepStr = strfind(filepath, filesep);
finalStr = strfind(filepath, strcat(filesep, 'final'));
ii = find(filesepStr == finalStr);
jj = max([1 ii-2]);
leadingStr = filepath(filesepStr(jj):numel(filepath));
newfile = fullfile(curdir, leadingStr, filename);
%%%%create directory if not exist
ij = strfind(leadingStr, filesep);
for (i = 1:numel(ij)-1)
    dirnow = strcat(curdir, leadingStr(1:(ij(i+1)-1)));
    if (exist(dirnow, 'dir')~=7) mkdir(dirnow); end
end
dirnow = strcat(curdir, leadingStr);
if (exist(dirnow, 'dir')~=7) mkdir(dirnow); end





