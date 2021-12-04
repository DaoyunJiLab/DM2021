function [behav,bhdata] = DataManager_FindAllBehav(datedir, behav)
%Find all behavioral data under the animal directory
%Fields assigned here:
behav.general.datedir = []; %%%date dir
behav.general.finaldir = []; %{str}
behav.general.sessID = []; %{str} %%%%this is actually a session ID
%behav.general.genotype = []; %{str}
%behav.general.age = []; %{str}

behav.general.sessname = []; % {sess}{name};
behav.general.sessstartT = []; % {sess}{starttime};
behav.general.sessendT = []; % {sess}{endtime};
behav.general.sesslength = []; % {sess}{length}
%%%behav.parm.sessFmarker, behav.sessBmarker, behav.parm.sessPmarker, behav.parm.sessionType, 
%%%        behav.parm.eventType,behav.parm.eventSession, behav.parm.eventPosltr
%%%        will be assigned in FindBehavParm
behav.general.eventname = []; %{sess}{{name1}; {name2}; } in a session
bhdata.event.eventtimes = []; %{sess}{{ev.start, ev.ent, ev.ref, ev.marker}} in a session
bhdata.pos.postimestamp = []; %{sess}[] 
bhdata.pos.XX = []; bhdata.pos.YY = []; %{sess}{j}j = color
%bhdata.pos.marker = []; %%%{i} marker can be 'greenA1' for the green diode of animal 1
                       %%%          'frontA1' for the front diode of animal 1
behav.general.posMarker = []; %%%{} 
behav.general.posAnimalname = []; %%%{} same as bhdata.pos.animalname, but conjugate all names
behav.general.posLocation = []; %%%{} same as bhdata.pos.location, but conjugate all locations
bhdata.pos.ltrfilename = []; %%%{.ltr1, .ltr2 ...}
bhdata.pos.posltr = []; %%%{joints, joint2, ...}

% finddddir = dir(animalpath);
% ndatedir = 0;
% datedir = [];
% for (j = 1:numel(finddddir))
% if (finddddir(j).isdir == 1) && (~strcmp(finddddir(j).name, '.')) && (~strcmp(finddddir(j).name, '..'))
%     ndatedir = ndatedir + 1;
%     datedir{ndatedir} = finddddir(j).name;
% end
% end

nsess = 0; ndatedir = numel(datedir); %start count of spikes found
%disp(['-----> number of dates ', num2str(ndatedir)]);
for (i = 1:ndatedir)
    disp(['---------> date: ', datedir{i}]); %datename = []; [datename, tok] = strtok(datedir{i}, '_'); 
    %finalpath = strcat(animalpath, filesep, datedir{i}, filesep, 'final');
    finalpath = datedir{i};
    if (exist(finalpath, 'dir') == 7) %if final dir exist
        infofile = fullfile(finalpath, 'animal.info'); Ttype=[]; Tfield=[]; Tcontent = []; 
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

        disp('------------> checking in all position linearization data ');
        posdir{1} = strcat(finalpath, filesep, 'pos');
        [ltrpath, ltrname] = GetAllFile(posdir, '', '.ltr'); posltr = []; 
        for (ikk = 1:numel(ltrpath))
            posltr{ikk} = ReadPosltrFile(fullfile(ltrpath{ikk}, ltrname{ikk}));
        end
        
        disp('------------> checking in all position data ');
        for (ikk = 1:numel(sessions.start))
            nsess = nsess + 1;
            behav.general.finaldir{nsess} = finalpath; %{str}
            datenownow = getdatedir(finalpath); %{str
            behav.general.datedir{nsess} = datenownow;
            behav.general.sessname{nsess} = sessions.marker{ikk};
            behav.general.sessID{nsess} = strcat(datenownow, '_', sessions.marker{ikk});
            behav.general.sessstartT{nsess} = sessions.start(ikk);
            behav.general.sessendT{nsess} = sessions.ent(ikk);
            behav.general.sesslength{nsess} = sessions.ent(ikk)-sessions.start(ikk);
            sessionflag = sessions.marker{ikk};
            [postimestamp, XX, YY, poscolor] = DE_CheckInAllPosition(finalpath, sessionflag);
            bhdata.pos.postimestamp{nsess} = postimestamp;
            bhdata.pos.XX{nsess} = XX; bhdata.pos.YY{nsess} = YY; behav.general.posMarker{nsess} = poscolor;
            behav.general.posAnimalname{nsess} = []; behav.general.posLocation{nsess} = [];
            %%%find animal name and location of each postion file in each session
            for (ttk = 1:numel(poscolor))
                behav.general.posAnimalname{nsess}{ttk} = []; behav.general.posLocation{nsess}{ttk} = [];
                [animname, loca, moreF, moreC] = findanimlocation(Ttype, Tfield, Tcontent, sessionflag, poscolor{ttk});
                if (~isempty(animname))
                   behav.general.posAnimalname{nsess}{ttk} = animname; behav.general.posLocation{nsess}{ttk} = loca; 
                %else
                %   behav.general.posAnimalname{nsess}{ttk} = animalpath; behav.general.posLocation{nsess}{ttk} = 'front';
                end
                for (kik = 1:numel(moreF)) behav.general.(moreF{kik}){nsess}{ttk} = moreC{kik}; end
            end
            %%%%find event files within a session
            ntev = 0; enowname = []; enowtimes = []; enowposltr = [];
            for (ttk = 1:numel(eventtimes))
                if ~isempty(eventtimes{ttk}.start)
                    if (min(eventtimes{ttk}.start) >= sessions.start(ikk)) && (max(eventtimes{ttk}.ent) <= sessions.ent(ikk))
                        ntev = ntev + 1; enowname{ntev} = eventname{ttk}; enowtimes{ntev} = eventtimes{ttk};
                    end
                end
            end
            %disp(['------------> found ' num2str(ntev) ' events in session ' sessions.marker{ikk}]); 
            behav.general.eventname{nsess} = enowname; bhdata.event.eventtimes{nsess} = enowtimes; 
            bhdata.pos.ltrfilename{nsess} = ltrname; bhdata.pos.posltr{nsess} = posltr;
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

function [animname, loca, moreF, moreC] = findanimlocation(Ttype, Tfield, Tcotent, sessionflag, poscolor)
animname = []; loca = []; moreF = []; moreC = []; flag = strcat(sessionflag, '_', poscolor);
iii = find( strcmpi(Ttype, 'position') & strcmpi(Tfield, flag) );
if (~isempty(iii))
        ii = iii(1);
%if (numel(ii) == 1);
        [str,tok] = strtok(Tcotent{ii}, '_'); animname = str; 
        if (numel(tok)>1) loca = tok(2:numel(tok)); end
        kk = find( strcmpi(Ttype, animname) ); moreF = Tfield(kk); moreC = Tcotent(kk);
%end
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





