function [pinfo,data] = DataManager_FindAllSpike(datedir, pinfo)

%%%try some new
%Find all spikes (full or partial) under the animal directory
%Fields assigned here:
        pinfo.general.recarea = []; pinfo.general.finaldir = []; pinfo.general.datedir = []; pinfo.general.animalname = [];  %{str}
        %pinfo.general.genotype = []; %{str} %pinfo.general.age = []; %{str}
        pinfo.general.clname = []; pinfo.general.TTname = []; pinfo.general.parmfile = []; pinfo.general.wavefile = []; %{str};
        data.spike.spiketime = []; %{[nspike]}  
        data.events.eventtimes = []; %{{ev.start, ev.ent, ev.ref, ev.marker}}
        pinfo.general.gain = []; %{[4]}
        pinfo.general.sessionname = []; % {{session1; session2;}};
        pinfo.general.sessionstartT = []; pinfo.general.sessionendT = []; pinfo.general.sessionlength = []; % {[l1 l2]};
        pinfo.general.eventname = []; %{{str1; str2; ...}}
% finddddir = dir(animalpath); ndatedir = 0; datedir = [];
% for (j = 1:numel(finddddir))
% if (finddddir(j).isdir == 1) && (~strcmp(finddddir(j).name, '.')) && (~strcmp(finddddir(j).name, '..'))
%     ndatedir = ndatedir + 1; datedir{ndatedir} = finddddir(j).name;
% end
% end

nspike = 0; %start count of spikes found
%disp(['-----> number of dates ', num2str(numel(datedir))]);
for (i = 1:numel(datedir))
    disp(['---------> date: ', datedir{i}]); %datename = []; [datename, tok] = strtok(datedir{i}, '_'); 
    finalpath = datedir{i}; %strcat(animalpath, filesep, datedir{i}, filesep, 'final');
    if (exist(finalpath, 'dir') == 7) %if final dir exist
        hisfile = fullfile(finalpath, 'tetrode.his');
        disp(['------------> reading tetrode histology file: ', hisfile]);
        [TTname, TTarea] = ReadHisFile(hisfile);
        
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
        
        disp(['------------> checking in all spike data ']);
        [spikefilename, recarea, spikedata, spikegain] = DE_CheckInAllSpike(finalpath, TTname, TTarea);
        
        nspikenow = numel(spikefilename);
        for (n = 1:nspikenow)
            nspike = nspike + 1;
            pinfo.general.finaldir{nspike} = finalpath; %{str}
            datenownow = getdatedir(finalpath); %{str}
            pinfo.general.datedir{nspike} = datenownow;
            pinfo.general.gain{nspike} = spikegain{n}; %{[4]}
            pinfo.general.recarea{nspike} = recarea{n}; %{str}
            pinfo.general.parmfile{nspike} = spikefilename{n}; %{str};
            pinfo.general.wavefile{nspike} = strrep(spikefilename{n}, '.spm', '.spw'); %{str};
            data.spike.spiketime{nspike} = spikedata{n}; %{[nspike]}       
            [TTnow, clnow] = findttclnames(spikefilename{n});
            pinfo.general.clname{nspike} = strcat(datenownow, '_', TTnow, '_', clnow); %{str};
            %pinfo.general.clname{nspike} = strcat(TTnow, '_', clnow); %{str};
            pinfo.general.TTname{nspike} = TTnow; %{str}; 
            pinfo.general.eventname{nspike} = eventname;
            data.events.eventtimes{nspike} = eventtimes;
            [animname, moreF, moreC] = findanimalname(Ttype, Tfield, Tcontent, TTnow);
            if (~isempty(animname))
                pinfo.general.animalname{nspike} = animname; 
            else
                pinfo.general.animalname{nspike} = finalpath;    
            end
            for (ikk = 1:numel(moreF)) pinfo.general.(moreF{ikk}){nspike} = moreC{ikk}; end
            pinfo.general.sessionname{nspike} = sessions.marker;
            pinfo.general.sessionstartT{nspike} = sessions.start;
            pinfo.general.sessionendT{nspike} = sessions.ent;
            pinfo.general.sessionlength{nspike} = sessions.ent-sessions.start;
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

function [animname, moreF, moreC] = findanimalname(Ttype, Tfield, Tcontent, TTnow)
animname = []; moreF = []; moreC = [];
ii = find( strcmpi(Ttype, 'tetrode') & strcmpi(Tfield, TTnow) );
if (numel(ii) == 1);
    animname = Tcontent{ii};
    kk = find( strcmp(Ttype, animname) );
    moreF = Tfield(kk); moreC = Tcontent(kk);
end

function [TTnow, clnow] = findttclnames(filename)
[dd, nn, ee] = fileparts(filename);
clnow = nn;
ttstr = strcat(filesep, 'TT');
aa = strfind(dd, ttstr);
TTnow = dd(aa+1:numel(dd));

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





