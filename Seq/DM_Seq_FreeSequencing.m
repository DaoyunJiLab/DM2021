function DM_Seq_FreeSequencing
%%Sequence analysis on a database (.spikedb), require a matching .tmp file
%%-- all free sequencing 
hf = gcbf; hgroup = getappdata(hf, 'hgroup'); ok = 1;
fname = get(hf, 'Name'); ftt = strfind(fname, '__'); currentfilename = fname(ftt+2:numel(fname));
[pp, nn, ee] = fileparts(currentfilename);
if ~strcmp(ee, '.spikedb')
    disp('-----> not a spikedb; aborted'); ok = 0;
else
    pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); behav = []; bhdata = [];
    plotparm = getappdata(hf, 'plotparm');
    hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
    for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end    
    seqparm.rankmode = 'OrderRank'; rrank = []; %%%not applicable at this time: %%%%%%%%Free sequences use rank order correlation;
                                                                                %%%%%%%%WithinEvents/SlidingWindows use pair-order similarity
    [MCroot, MCname, DAname, DEname, ANname, MAname] = CurrentVersion;
    rfile = fullfile(MCroot, 'DataExplorer', 'Sequence', 'rrankprob.mat'); 
    if (exist(rfile, 'file') == 2)
          S = load(rfile); rrank = S.rrank; S = [];
    else
          ok = 0; disp('-----> ranking probability file for selected rank mode not found; aborted');
    end
end

%%%%The following is to sequence the current .spikedb using multiple sets of templates: replace back to the single template file after done.
if ok
   [ttfile, ttfilename, ttfilepath, ttfileext, okk] = FileDlg(cd, '.tmp');  %get a group of template files
   if okk
       alltmpfile = ttfile;
   else
       ok = 0;
   end
end
% if ok
%     oTT = questdlg('Open a template file?', 'Template input', 'Yes');
%     if (strcmp(oTT, 'Yes')) %
%         cpathname = fullfile(cd, '*.tmp');
%         [fname, pname] = uigetfile(cpathname, 'Select a template file to open:');
%         if fname ~= 0
%            tmpfile = fullfile(pname, fname);
%            disp(['-----> template file: ', tmpfile]); seqparm.tmpfile = tmpfile;
%            S = load(tmpfile, '-mat'); tmp = S.pinfo; S = [];
%            disp(['--------> number of templates found: ', num2str(numel(tmp))]);
%         else
%             ok = 0;
%         end
%     elseif (strcmp(oTT, 'No')) 
%         seqparm.tmpfile = []; tmp = [];
%     else
%         ok = 0;
%     end
% end
if ok
    seqparm.sessionoption = 'FreeSequences'; seqparm.eventoption = 'FreeSequences';
    seqparm.sessiontimingoption = 'RatePeak_time'; seqparm.eventtimingoption = 'RatePeak_time'; 
    seqparm.timeres = 0; seqparm.spatres = 0; seqparm.timesigma = 0; seqparm.spatsigma = 0;
    seqparm.Nsigma = 5; seqparm.timePthre = 0; seqparm.spatPthre = 0;
    sss = {'RatePeak_time'; 'RatePeak_space'};
    [S, ok] = listdlg('ListString', sss, 'PromptString', 'Event timing option');
    if ok
          seqparm.eventtimingoption = sss{S};
          input = inputdlg({'Enter temporal resolution (s)'; 'Enter smoothing sigma (s)'; 'Enter peak threshld (std)'}, ...
                  'Temporal timing parameter Input', 3, {'0.005'; '0.5'; '0'}); 
          if (~isempty(input))
                 seqparm.timeres = str2num(input{1}); seqparm.timesigma = str2num(input{2}); seqparm.timePthre = str2num(input{3});
          else
                 ok = 0;
          end
          if ~isempty(strfind(seqparm.eventtimingoption, '_space')) %%%%if use peak locations of spatial rate curves
              if plotparm.linkbehav == 1
                 behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
              else
                 ok = 0; disp('-----> No linked behavdb; aborted');
              end
              input = inputdlg({'Enter smoothing sigma (pixels)'; 'Enter peak threshld (std)'}, ...
                  'Spatial timing parameter Input', 2, {'10'; '0.1'}); 
              if (~isempty(input))
                 seqparm.spatres = 10; seqparm.spatsigma = str2num(input{1}); seqparm.spatPthre = str2num(input{2});
              else
                 ok = 0;
              end
          end
    end
end
if ok && (~(strcmp(seqparm.eventtimingoption, 'RatePeak_space') || strcmp(seqparm.eventtimingoption, 'Rate1stPeak_space'))) 
    seqparm.setbacktime = 0;
    input = inputdlg({'Enter backtime (s)'}, 'Backtime input', 1, {'0'}); 
    if (~isempty(input))
        seqparm.setbacktime = str2num(input{1}); 
    else
        ok = 0;
    end
end
if ok 
   [fname, pname] = uiputfile(fullfile(cd, '.seqdb')); 
   if (fname ~= 0)
       seqfile = fullfile(pname, fname);
   else
       ok = 0;
   end
end
if ok
   input = inputdlg({'Enter max sequence length (x tmplate length)'; 'Enter track max sequence duration (s)';...
       'Enter open/sleep max sequence duration (s)'; 'Enter max gap (x template duration)'},...
           'Sequence leng options', 4, {'2'; '10'; '20'; '3'}); 
   if (~isempty(input))
           seqparm.maxtmpleng = str2num(input{1}); seqparm.maxtimeL = str2num(input{2}); seqparm.maxtimeO = str2num(input{3});
           seqparm.maxgap = str2num(input{4}); %%%parameters used to break down long sequences
   else
           ok = 0;
   end
end
if ok
    seqparm.Nshuffle = 100; seqparm.shufflemode = 'GlobalID'; seqparm.significancelevel = 0.05;
    disp(['---> number of sets of templates to sequence: ', num2str(numel(alltmpfile))]);
    for itt = 1:numel(alltmpfile)
        disp(['----> ']);
        disp(['----> template file now: ', alltmpfile{itt}]); 
        seqparm.tmpfile = alltmpfile{itt};
        S = load(alltmpfile{itt}, '-mat'); tmp = S.pinfo; S = [];
        [pp, nn, ee] = fileparts(alltmpfile{itt}); [ppp, nnn, eee] = fileparts(seqfile);
        filenamenow = fullfile(ppp, strcat([nnn, '_', nn], eee));
        disp(['--------> will save to seqdb file: ', filenamenow]);
        disp(['--------> number of templates found: ', num2str(numel(tmp))]);
        seq = []; seq.work = struct([]); 
        [seq, seqdata] = computeseqdb(seq, tmp, pinfo, data, behav, bhdata, cellind, seqparm, rrank);  
        if isfield(seq, 'general') && (~isempty(seq.general.tmpID))
           %%%%%group data
           seqdata.grouplist.groupname{1} = 'List0'; seqdata.grouplist.groupindex{1} = 1:numel(seq.general.tmpID);
           seqdata.grouplist.grouptype{1} = 'Manual'; seqdata.grouplist.groupcrit{1} = []; seqdata.grouplist.groupparents{1} = []; 
           seqdata.parentfile = currentfilename; 
           pppinfo = pinfo; dddata = data;
           pinfo = seq; data = seqdata; save(filenamenow, 'pinfo', 'data', '-mat');
           pinfo = []; pinfo = pppinfo; data = []; data = dddata;
        end
    end
end
disp('************************');

function [seq, seqdata] = computeseqdb(seq, tmp, pinfo, data, behav, bhdata, cellind, seqparm, rrank)
seqdata = []; nt = 0;
%%%%search for all the possible finaldir
allfinaldir = pinfo.general.finaldir(cellind); 
unifinaldir = unique(allfinaldir);
disp(['-----> number of final dirs: ', num2str(numel(unifinaldir))]);
for (i = 1:numel(unifinaldir))
    fdirnow = unifinaldir{i};
    disp(['-------> final dir now: ', fdirnow]);
    ind = find(strcmp(allfinaldir, unifinaldir{i}));
    neuronid = cellind(ind); %%%all cells in this final directory 
    cellparmfile = pinfo.general.parmfile(neuronid); 
    if (numel(cellparmfile)>62) %%%
        disp('----------> number of cells exceeds template capacity; try to be more selective; aborted');
    elseif (numel(cellparmfile)<4)
        disp('----------> too few cells; aborted');
    else
        animaldatenow = strcat(pinfo.general.animalname{neuronid(1)}, '_', pinfo.general.datedir{neuronid(1)});
        tmpnow = [];
        if isempty(tmp) %%%use default, blank template
            cellletter = DE_Seq_AssignLetter(numel(cellparmfile), 'allletters'); 
            tmpnow.tmp.tmpfileletter = cellletter; 
            tmpnow.tmp.tmpfilename = cellparmfile; tmpnow.tmp.tmpfinalseq{1} = cellletter; 
        else
            tmpnow = findtemplates(tmp, animaldatenow);
        end 
        %%%%check tmp length
        nlet = zeros(1, numel(tmpnow));
        for (tt = 1:numel(tmpnow))
            nlet(tt) = numel(tmpnow(tt).tmp.tmpfinalseq{1});
        end
        if (max(nlet)>=4)
           disp('----------> sequencing data ... ');
           [sessseq, evseq] = sequencedata_raw(pinfo, data, behav, bhdata, neuronid, seqparm, fdirnow); %sessseq.cellpeaktime{i=session}{j=cell}
           disp('----------> analyzing data sequences ... ');
           for (j = 1:numel(tmpnow)) %%%for each template
               if nlet(j) >= 4
                  nt = nt + 1; %%%entry of the seqdb items
                  [sessseq, evseq, tmpnow(j)] = analyzesequencedata(pinfo, data, sessseq, evseq, neuronid, tmpnow(j), seqparm, rrank);
                  tmprank = tmpnow(j).tmp.tmpfinalrank{1}; 
                  [seqnow, seqdatanow] = computeseqproperties(sessseq, evseq, seqparm, pinfo, data, neuronid(1), tmpnow(j), tmpnow, rrank, tmprank);
                  [seq, seqdata] = assignseqoutput(seq, seqdata, seqnow, seqdatanow, nt);
               else
                  disp(['-----------> template too short; skipping: ', tmpnow(j).evtfile, '(', tmpnow(j).tmp.tmpfinalseq{1}, ')(', num2str(tmpnow(j).tmp.finalrank{1}), ')']);
               end
           end
        else
           disp(['----------> no high-order templates found; skipping dir: ', fdirnow]); 
        end
    end
end

function [sessseq, evseq] = sequencedata_raw(pinfo, data, behav, bhdata, neuronid, seqparm, fdirnow) %sessseq.cellpeaktime{i=session}{j=cell}
%%%1. select out cells; 2. if temporal/spatial peaks, work out session
%%%smoothed rate curves; 3. pick the peak
%%%%%%%%output: sessseq.cellpeaktime{i}{j} i= nsess; j = ncell
aniname = pinfo.general.animalname{neuronid(1)};
sessName = pinfo.general.sessionname{neuronid(1)}; sessType = pinfo.parm.sessType{neuronid(1)}; 
sessST = pinfo.general.sessionstartT{neuronid(1)}; sessET = pinfo.general.sessionendT{neuronid(1)}; nsess = numel(sessName); ncell = numel(neuronid);
sessseq = cell(1, nsess);
for (i = 1:nsess) %%%for session data the only available timingoption is 'RatePeak_time'
    Tsigma = seqparm.timesigma; TPth = seqparm.timePthre;
    if strcmp(sessType{i}, 'open') || strcmp(sessType{i}, 'sleep') %%%for open/sleep session
        Tsigma = 2*seqparm.timesigma; TPth = 0; %seqparm.timePthre;
    end
    cellpeaktime = cell(1, ncell);
    %%%%%work out the firing rate
    timepoint = (sessST(i):seqparm.timeres:sessET(i));
    for (j = 1:ncell)
        ratedata{1} = getfiringrates_binspikes(data.spike.spiketime{neuronid(j)}, Tsigma, timepoint, seqparm.timeres, sessST(i), sessET(i));
        [timePeakTh, ratedata, ~, ~] = normalizefindthreshold(ratedata, TPth);
        [~, maxind] = FindLocalMinima(ratedata{1}); maxV = ratedata{1}(maxind);
        cellpeaktime{j} = timepoint(maxind(maxV > timePeakTh(1))); 
    end
    sessseq{i}.cellpeaktime = cellpeaktime;
end
evName = pinfo.general.eventname{neuronid(1)}; nev = numel(evName);
evseq= cell(1, nev); 
if strcmp(seqparm.eventtimingoption, 'RatePeak_space')
    evType = pinfo.parm.eventtype{neuronid(1)}; evT = data.events.eventtimes{neuronid(1)}; nev = numel(evName); 
    evsessid = zeros(1, nev); evsessName = cell(1, nev); allevtime = cell(1, nev);
    sigmanow = seqparm.spatsignma; Nsigmanow = seqparm.Nsigma; %=pinfo.fSmooth1DNsigma(neuronid(1));
    for (i = 1:nev)
        cellpeakloc = cell(1, ncell);
        if (strcmp(evType{i}, 'run'))
           [evsessName{i}, evsessid(i)] = identifysession(evT{i}, sessName, sessST, sessET);%%%%work out which session the event is
           for (j = 1:ncell)
              [spatrate{1}, spatX] = findspatialrate(data.spike.spiketime{neuronid(j)}, behav, bhdata,...
                  fdirnow, evsessName{i}, evName{i}, evT{i}, sigmanow, Nsigmanow);
              [spatPeakTh, spatrate, ~, ~] = normalizefindthreshold(spatrate, seqparm.spatPthre);
              for (k = 1:numel(evT{i}.start))
                  ratedata = spatrate{1}{k};
                  [~, maxind] = FindLocalMinima(ratedata); maxV = ratedata(maxind);
                  cellpeakloc{j}{k} = spatX(maxind(maxV > spatPeakTh(1)));
              end
           end
        end
        evseq{i}.cellpeakloc = cellpeakloc;
    end
end

function [sessseq, evseq, tmpnow] = analyzesequencedata(pinfo, data, sessseq, evseq, neuronid, tmpnow, seqparm, rrank)
%%%3. for each event, sessions, workout the ep scheme and break sequences accordingly
%%%%%%%%output: seq{ii,jj} for each event ep{ii,jj}
%%%%%first select which template sequence to use
cellletter = tmpnow.tmp.tmpfileletter; cellfilename = tmpnow.tmp.tmpfilename; 
alltmpseq = tmpnow.tmp.tmpfinalseq;
finalfiterr = []; if isfield(tmpnow.tmp, 'finalfiterr') finalfiterr = tmpnow.tmp.finalfiterr; end %%%select the best fit template to use
if ~isempty(finalfiterr)
   [minfiterr, sleind] = min(finalfiterr);
else
    ntmp = numel(alltmpseq); sleind = ceil(rand*ntmp); if (sleind ==0) sleind = 1; end 
end
tmpseq = alltmpseq{sleind}; tmprank = 1:numel(tmpseq);
if isfield(tmpnow.tmp, 'finalrank') 
    ccrank = tmpnow.tmp.finalrank{sleind}; cu = unique(ccrank);
    for (i = 1:numel(cu))
         iii = find(ccrank == cu(i));
         if (numel(iii) > 1)
             mm= mean(tmprank(iii)); tmprank(iii) = mm*ones(1, numel(iii)); 
         end
    end
    tmprank = ccrank; 
end 
%%%%%match letters in the sequence
nlet = numel(cellletter); intemp = zeros(1, nlet);
for (i = 1:nlet) 
    if ~isempty(find(tmpseq == cellletter(i))) intemp(i) = 1; end
end
iii = find(intemp); cellletter = cellletter(iii); cellfilename = cellfilename(iii);
nlet = numel(cellfilename); isok = zeros(1, nlet); indok = zeros(1, nlet);
for (i = 1:nlet)
    ti = find( strcmp(pinfo.general.parmfile(neuronid), cellfilename{i}) );
    if (numel(ti) == 1)
        isok(i) = 1; indok(i) = ti;
    end
end
iii = find(isok == 1); 
nlet = numel(iii); cellletter = cellletter(iii); cellfilename = cellfilename(iii); cellid = indok(iii); %%%this contains the cell indices in sesseq, evseq to analyze
tmpisok = zeros(1, numel(tmpseq));
for (i = 1:numel(tmpseq))
    ti = find(cellletter == tmpseq(i));
    if (numel(ti) == 1) tmpisok(i) = 1; end
end
kkk = find(tmpisok==1); tmpseq =tmpseq(kkk); tmprank = tmprank(kkk);
tmpnow.tmp.tmpfileletter = cellletter; tmpnow.tmp.tmpfilename = cellfilename; 
tmpnow.tmp.tmpfinalseq = []; tmpnow.tmp.tmpfinalrank = []; tmpnow.tmp.finalfiterr = [];
tmpnow.tmp.tmpfinalseq{1} = tmpseq; tmpnow.tmp.tmpfinalrank{1} = tmprank; tmpnow.tmp.finalfiterr{1} = minfiterr;
disp(['-------------> template now: ', tmpnow.evtfile, '(', tmpseq, ')(', num2str(tmprank), ')']);
%%%%%%%working on session sequences
sessName = pinfo.general.sessionname{neuronid(1)}; sessType = pinfo.parm.sessType{neuronid(1)};
sessST = pinfo.general.sessionstartT{neuronid(1)}; sessET = pinfo.general.sessionendT{neuronid(1)}; nsess = numel(sessName);
allletters= cell(1, nsess); allpeaks = cell(1, nsess);
for (i = 1:nsess) %%%for session data the only available timingoption is 'RatePeak_time'
    seq = []; seqstart = []; seqend = []; seqmarker = []; matchscore = []; posmatchprob = []; negmatchprob = []; seqtime = [];
    %%%%work out the ep structure (sliding windows)
    ep = []; nep = 1;
    for (ti = 1:nep)
             ep{ti,1}.start = sessST(i); ep{ti, 1}.ent = sessET(i); ep{ti,1}.marker = sessName{i};
    end
    %%%%work out the sequences
    [mm,nn] = size(ep);
    alllettersnow = []; allpeaksnow = [];
    for (j = 1:numel(cellid))
        peaknow = sessseq{i}.cellpeaktime{cellid(j)};
        allpeaksnow = [allpeaksnow peaknow]; alllettersnow = [alllettersnow cellletter(j)*ones(1, numel(peaknow))];
    end
    [allpeaks{i}, iii] = sort(allpeaksnow); allletters{i} = alllettersnow(iii);
        %%then break down the sequence to shoter ones with each event
    for (ii = 1:mm) %%%original event index
        for (jj = 1:nn)
            [seq{ii,jj}, seqstart{ii,jj}, seqend{ii,jj}, seqmarker{ii,jj}, matchscore{ii,jj}, posmatchprob{ii,jj}, negmatchprob{ii,jj}, seqtime{ii,jj}] = ...
               breakallsequences_nonselect(allletters{i}, allpeaks{i}, ep{ii,jj}, tmprank, rrank, tmpseq, seqparm, sessType{i});
        end
    end
    sessseq{i}.ep = ep; 
    sessseq{i}.seq = seq; sessseq{i}.seqstart = seqstart; sessseq{i}.seqend = seqend; sessseq{i}.seqmarker = seqmarker; sessseq{i}.seqtime = seqtime;
    sessseq{i}.matchscore = matchscore; sessseq{i}.posmatchprob = posmatchprob; sessseq{i}.negmatchprob = negmatchprob;
end
%%%%%%%working on sevent sequences
evName = pinfo.general.eventname{neuronid(1)}; evType = pinfo.parm.eventtype{neuronid(1)};
evT = data.events.eventtimes{neuronid(1)}; nev = numel(evName); 
evsessid = zeros(1, nev); evsessName = cell(1, nev); 
for (i = 1:nev) %%%for ev data the multiple available timingoption
    seq = []; seqstart = []; seqend = []; seqmarker = []; matchscore = []; posmatchprob = []; negmatchprob = []; seqtime = [];
    [evsessName{i}, evsessid(i)] = identifysession(evT{i}, sessName, sessST, sessET);%%%%work out which session the event is
    if isempty(strfind(seqparm.eventtimingoption, 'space'))
        evT{i}.start = evT{i}.start - seqparm.setbacktime; evT{i}.ent = evT{i}.ent + seqparm.setbacktime; 
    end
    %%%%work out the ep structure (sliding windows)
    evTime = evT{i}; ep = []; nep = numel(evTime.start);
    if ~isempty(strfind(seqparm.eventtimingoption, 'space')) %%%if spatial mode, evt changed to spatial windows
        if (~isempty(spatX))
            evTime.start = min(spatX)*ones(nep,1); evTime.ent = max(spatX)*ones(nep,1); 
        end
    end
    for (ti = 1:nep)
             ep{ti,1}.start = evTime.start(ti); ep{ti, 1}.ent = evTime.ent(ti); ep{ti,1}.marker = evTime.marker{ti};
    end
    %%%%work on the sequences: all different options
    [mm,nn] = size(ep);
    if strcmp(seqparm.eventoption, 'FreeSequences') && strcmp(seqparm.eventtimingoption, 'RatePeak_time')
        alltimenow = allpeaks{evsessid(i)};
        %%then break down the sequence to shoter ones with each event
        for (ii = 1:mm) %%%original event index
            %ik = find( (alltimenow>=evTime.start(ii)) & (alltimenow<=evTime.ent(ii)) ); 
            %allevletter = allletters{evsessid(i)}(ik); allevtime = alltimenow(ik);
        for (jj = 1:nn)
            [seq{ii,jj}, seqstart{ii,jj}, seqend{ii,jj}, seqmarker{ii,jj}, matchscore{ii,jj}, posmatchprob{ii,jj}, negmatchprob{ii,jj}, seqtime{ii,jj}] = ...
               breakfromsessionallsequences(sessseq{evsessid(i)}, ep{ii,jj});
        end
        end
    elseif strcmp(seqparm.eventoption, 'FreeSequences') && strcmp(seqparm.eventtimingoption, 'RatePeak_space')
        for (ii = 1:mm) %%%original event index
            alllevetters = []; allevpeaks = [];
            for (j = 1:numel(cellid))
                peaknow = evseq{i}.cellpeakloc{cellid(j)}{ii};
                allpeaksnow = [allpeaksnow peaknow]; alllettersnow = [alllettersnow cellletter(j)*ones(1, numel(peaknow))];
            end
            [allevpeaks, iii] = sort(allevpeaks); allevletters = allevlettersnow(iii);
        for (jj = 1:nn)
            [seq{ii,jj}, seqstart{ii,jj}, seqend{ii,jj}, seqmarker{ii,jj}, matchscore{ii,jj}, posmatchprob{ii,jj}, negmatchprob{ii,jj}, seqtime{ii,jj}] = ...
               breakallsequences_nonselect(allevletters, allevpeaks, ep{ii,jj}, tmprank, rrank, tmpseq, seqparm, 'linear');
        end
        end   
    end
    evseq{i}.ep = ep; 
    evseq{i}.seq = seq; evseq{i}.seqstart = seqstart; evseq{i}.seqend = seqend; evseq{i}.seqmarker = seqmarker; evseq{i}.seqtime = seqtime;
    evseq{i}.matchscore = matchscore; evseq{i}.posmatchprob = posmatchprob; evseq{i}.negmatchprob = negmatchprob;
end

function [seq, seqdata] = computeseqproperties(sessseq, evseq, seqparm, pinfo, data, neuronid, tmpnow, alltmp, rrank, tmprank)
%%%%now work on what sequence properties to compute: 1-template properties; 2-sequence properties; 3-significance
%%%%%%%use default parameters, can be recomputed later
%%%%%general variables
seq.general.tmpID = strcat(tmpnow.parm.animaldate, '_', tmpnow.evtfile); 
seq.general.recarea = pinfo.general.recarea{neuronid}; seq.general.finaldir = pinfo.general.finaldir{neuronid};
seq.general.datedir = pinfo.general.datedir{neuronid}; seq.general.animalname = pinfo.general.animalname{neuronid};
seq.general.sessionname = pinfo.general.sessionname{neuronid}; seq.general.sessionstartT = pinfo.general.sessionstartT{neuronid};
seq.general.sessionendT = pinfo.general.sessionendT{neuronid}; seq.general.sessionlength = pinfo.general.sessionlength{neuronid};
seq.general.eventname = pinfo.general.eventname{neuronid}; seq.general.genotype = pinfo.general.genotype{neuronid};
seq.general.sex = pinfo.general.sex{neuronid}; seq.general.age = pinfo.general.age{neuronid};
%%%%parameters
seq.parm = seqparm; seq.parm = rmfield(seq.parm, 'tmpfile'); seq.parm.seqtype = 'SeqMatching';
seq.parm.sessType = pinfo.parm.sessType{neuronid}; seq.parm.eventType = pinfo.parm.eventtype{neuronid};
%%%%%%%%%1.template properties - quality and similarities with other templates in the same finaldir
seq.tmp.tmpfile = seqparm.tmpfile; 
seq.tmp.tmpseq = tmpnow.tmp.tmpfinalseq{1}; seq.tmp.tmprank = tmpnow.tmp.tmpfinalrank{1};
seq.tmp.cellnum = numel(tmpnow.tmp.tmpfinalseq{1}); seq.tmp.length = max(seq.tmp.tmprank)-min(seq.tmp.tmprank);
seq.tmp.crrmode = tmpnow.parm.crrmode; seq.tmp.sessevname = tmpnow.evtfile;
seq.tmp.crrmode = tmpnow.parm.crrmode; 
if isfield(tmpnow.parm, 'sessev') 
    sessev = tmpnow.parm.sessev;
else
    sessev = 'unknown';
end
seq.tmp.sessev = sessev;
seq.tmp.filename = tmpnow.tmp.tmpfilename;  nlet = numel(tmpnow.tmp.tmpfileletter);
for (i = 1:nlet) seq.tmp.fileletter{i} = tmpnow.tmp.tmpfileletter(i); end
seq.tmp.pairscore = tmpnow.tmp.tmpfinalscore/numel(tmpnow.tmp.tmpfinalorderpair); 
seq.tmp.orderscore = tmpnow.tmp.tmpfinalscore/(nlet*(nlet-1)/2); 
seq.tmp.neighborscore = checkallpairs(tmpnow.tmp.tmpfinalseq{1}, tmpnow.tmp.tmpfinalorderpair);
seq.tmp.rankfiterr = tmpnow.tmp.finalfiterr;
[othertmpname,relationtype,matchscore, matchposprob, matchnegprob] = findrelationwithothertmps(tmpnow,alltmp, rrank);
seq.tmp.othertempname = othertmpname; seq.tmp.relationtype = relationtype;
seq.tmp.matchscore = matchscore; seq.tmp.posmatchprob = matchposprob; seq.tmp.negmatchprob = matchnegprob;
%%%%%%%%2.Sequence properties & significance
 
%%%%%%matchnumbers and theoretical significance 
nsess = numel(seq.general.sessionname); %%%%for sessions
posN = NaN*ones(1, nsess); posRate = NaN*ones(1, nsess); negN = NaN*ones(1, nsess); negRate = NaN*ones(1, nsess); 
for (i = 1:nsess)
    seqdata.data.sessseq{i} = sessseq{i};
    [posN(i), ~, ~, negN(i), ~, ~, ~, ~] = findtheorysig(sessseq{i}, seqparm.significancelevel, rrank, seqparm.sessionoption);
    posRate(i) = posN(i)/seq.general.sessionlength(i); negRate(i) = negN(i)/seq.general.sessionlength(i); 
end
seq.seq.sessPosMatchN = posN; seq.seq.sessPosRate = posRate; seq.seq.sessNegMatchN = negN; seq.seq.sessNegRate = negRate; 
nev = numel(seq.general.eventname); %%%%%for events
posN = NaN*ones(1, nev); negN = NaN*ones(1, nev); 
posTimeRate = NaN*ones(1, nev); negTimeRate = NaN*ones(1, nev); posLapRate = NaN*ones(1, nev); negLapRate = NaN*ones(1, nev);
evTimes = data.events.eventtimes{neuronid(1)}; seqdata.events.eventimes = evTimes;  
for (i = 1:nev)
    seqdata.data.evseq{i} = evseq{i};
    [posN(i), ~, ~, negN(i), ~, ~, ~, ~] = findtheorysig(evseq{i}, seqparm.significancelevel, rrank, seqparm.eventoption);
    totlap = numel(evTimes{i}.start); tottime = sum(evTimes{i}.ent-evTimes{i}.start);
    posTimeRate(i) = posN(i)/tottime; negTimeRate(i) = negN(i)/tottime;
    posLapRate(i) = posN(i)/totlap; negLapRate(i) = negN(i)/totlap;
end
seq.seq.evPosMatchN = posN; seq.seq.evPosTimeRate = posTimeRate; seq.seq.evPosLapRate = posLapRate; 
seq.seq.evNegMatchN = negN; seq.seq.evNegTimeRate = negTimeRate; seq.seq.evNegLapRate = negLapRate;
%%%%%shuffle significance
posP = NaN*ones(1, nsess); posE = NaN*ones(1, nsess); negP = NaN*ones(1, nsess); negE = NaN*ones(1, nsess); 
posRate = NaN*ones(1, nsess); negRate = NaN*ones(1, nsess); posSig = NaN*ones(1, nsess); negSig = NaN*ones(1, nsess);
for (i = 1:nsess) %for sessions
    [posP(i), posE(i), negP(i), negE(i), posSig(i), negSig(i)] = findshufflesig(seq.tmp.tmpseq, seq.seq.sessPosMatchN(i), seq.seq.sessNegMatchN(i),...
        seq.general.sessionstartT(i), seq.general.sessionendT(i), tmpnow.tmp.tmpfileletter,....
        sessseq{i}, rrank, seqparm, 'session', tmprank); %, alltimenow);
    posRate(i) = posE(i)/seq.general.sessionlength(i); negRate(i) = negE(i)/seq.general.sessionlength(i); 
end
seq.seq.sessPosMatchNshufsig = (seq.seq.sessPosMatchN - posE)./posSig; seq.seq.sessPosMatchShufP = posP; seq.seq.sessPosShufExpN = posE; 
seq.seq.sessPosShufExpSig = posSig; seq.seq.sessShufposRate = posRate;
seq.seq.sessNegMatchNshufsig = (seq.seq.sessNegMatchN - negE)./negSig; seq.seq.sessNegMatchShufP = negP; seq.seq.sessNegShufExpN = negE; 
seq.seq.sessNegShufExpSig = negSig; seq.seq.sessShufnegRate = negRate;
posP = NaN*ones(1, nev); posE = NaN*ones(1, nev); negP = NaN*ones(1, nev); negE = NaN*ones(1, nev); 
posTimeRate = NaN*ones(1, nev); negTimeRate = NaN*ones(1, nev); posLapRate = NaN*ones(1, nev); negLapRate = NaN*ones(1, nev); 
posSig = NaN*ones(1, nev); negSig = NaN*ones(1, nev);
for (i = 1:nev) %for events
    [~, evsessid] = identifysession(evTimes{i}, seq.general.sessionname, seq.general.sessionstartT, seq.general.sessionendT);%%%%work out which session the event is
    [posP(i), posE(i), negP(i), negE(i), posSig(i), negSig(i)] = findshufflesig(seq.tmp.tmpseq, seq.seq.evPosMatchN(i), seq.seq.evNegMatchN(i), ...
        seq.general.sessionstartT(evsessid), seq.general.sessionendT(evsessid), tmpnow.tmp.tmpfileletter,....
        evseq{i}, rrank, seqparm, 'event', tmprank); %, alltimenow);
    totlap = numel(evTimes{i}.start); tottime = sum(evTimes{i}.ent-evTimes{i}.start);
    posTimeRate(i) = posE(i)/tottime; negTimeRate(i) = negE(i)/tottime;
    posLapRate(i) = posE(i)/totlap; negLapRate(i) = negE(i)/totlap;
end
seq.seq.evPosMatchNshufsig = (seq.seq.evPosMatchN - posE)./posSig; seq.seq.evPosMatchShufP = posP; seq.seq.evPosShufExpN = posE; 
seq.seq.evPosShufExpSig = posSig; seq.seq.evShufposTimeRate = posTimeRate; seq.seq.evShufposLapRate = posLapRate;
seq.seq.evNegMatchNshufsig = (seq.seq.evNegMatchN - negE)./negSig; seq.seq.evNegMatchShufP = negP; seq.seq.evNegShufExpN = negE; 
seq.seq.evNegShufExpSig = negSig; seq.seq.evShufnegTImeRate = negTimeRate; seq.seq.evShufnegLapRate = negLapRate;

function [seq, seqstart, seqend, seqmarker, matchscore, posprob, negprob, seqtime] = breakfromsessionallsequences(sessseq, event)
iii = find( ((sessseq.seqstart{1,1}>=event.start)&(sessseq.seqstart{1,1}<=event.ent)) &  ((sessseq.seqend{1,1}>=event.start)&(sessseq.seqend{1,1}<=event.ent)));
seq = sessseq.seq{1,1}(iii); seqstart = sessseq.seqstart{1,1}(iii); seqend = sessseq.seqend{1,1}(iii); seqmarker = sessseq.seqmarker{1,1}(iii);
matchscore = sessseq.matchscore{1,1}(iii); posprob = sessseq.posmatchprob{1,1}(iii); negprob = sessseq.negmatchprob{1,1}(iii); seqtime = sessseq.seqtime{1,1}(iii);

function frac = checkallpairs(bestseq, orderpair)
frac = 0; nlet = numel(bestseq); 
if (nlet > 1)
   mm = zeros(1,nlet-1);
   for (i = 1:nlet-1)
        pairnow = bestseq(i:i+1); if (~isempty(find(strcmp(orderpair, pairnow)))) mm(i)=1; end
   end
   frac = numel(find(mm==1))/(nlet-1);
end

function [seq, seqstart, seqend, seqmarker, score, posprob, negprob, seqtime] = ...
               breakallsequences_nonselect(allletter, alltime, event, tmprank, rrank, tmpseq, seqparm, sessType)
Lmaxchecktime = seqparm.maxtimeL; Omaxchecktime = seqparm.maxtimeO; maxLeng = seqparm.maxtmpleng; maxgap = seqparm.maxgap;
%%%%%break long
iii = find( (alltime>=event.start) & (alltime<=event.ent) ); 
nletter = numel(iii); maxcheckleng = min([nletter maxLeng*numel(tmpseq)]); 
letternow = allletter(iii); timenow = alltime(iii); 
seq = []; seqstart = []; seqend = []; Ctime = []; nseq = 0;
%%%%work out the order rank
rankorder = 1:numel(tmprank); rankorder = resolvesimplerank(rankorder, tmprank);

if strcmp(sessType, 'linear')
    maxchecktime = Lmaxchecktime;
    for (leng = 4:maxcheckleng) %%%sequence length under consideration
    for (j = 1:nletter-leng+1)
        st = timenow(j); et = timenow(j+leng-1); T = timenow(j:j+leng-1);
        if (et-st<maxchecktime) && ((et-st)<=maxgap*(max(tmprank)-min(tmprank))) && (max(abs(diff(T)))<=maxgap*max(abs(diff(tmprank))))
            C = letternow(j:j+leng-1); 
            if (numel(unique(C)) >= 4)
               nseq = nseq + 1; seq{nseq} = C; Ctime{nseq} = T; seqstart(nseq) = st; seqend(nseq) = et;
            end
        end
    end
    end
else %%%%do the following for open/sleep sessions
    maxchecktime = Omaxchecktime;
    for (leng = 4:maxcheckleng) %%%sequence length under consideration
    for (j = 1:nletter-leng+1)
        st = timenow(j); et = timenow(j+leng-1);
        if (et-st<maxchecktime)
            nseq = nseq + 1; seq{nseq} = letternow(j:j+leng-1); Ctime{nseq} = timenow(j:j+leng-1); seqstart(nseq) = st; seqend(nseq) = et;
        end
    end
    end
end

nseq = numel(seq); score = NaN*ones(nseq,1); posprob = NaN*ones(nseq,1); negprob = NaN*ones(nseq,1); seqmarker = cell(nseq, 1);
if ~isempty(tmprank)
   for (i = 1:nseq)
       %[score(i), posprob(i), negprob(i)] = findseqscoreprob(seq{i}, tmpseq, rrank);
       [score(i), posprob(i), negprob(i)] = findseqscoreprob_simplerankorder(seq{i}, Ctime{i}, tmpseq, rankorder); %%%%
       %[score(i), posprob(i), negprob(i)] = findseqscoreprob_multiple(seq{i}, Ctime{i}, tmpseq, tmprank, rrank, seqparm); %%%%
   end
end
for (i = 1:nseq) seq{i} = char(seq{i}); seqmarker{i} = event.marker; end
seq = seq'; seqstart = seqstart'; seqend = seqend'; seqtime = Ctime';
function [score, posprob, negprob] = findseqscoreprob_simplerankorder(Cstr, Ctime, tmpseq, rankorder)%%%%
score = NaN; posprob = NaN; negprob = NaN; Cut = unique(Cstr); nlet = numel(Cut); neq = numel(Cstr);
if (nlet >= 4) %%%if unique letters more than 4
   %%%%work out str ranks
   Cstrtmp = zeros(1, neq); for (i = 1:neq) Cstrtmp(i) = rankorder(tmpseq == Cstr(i)); end
   Cstrrank = zeros(1, neq); [~, iii] = sort(Cstrtmp); Cstrrank(iii) = 1:neq; Cstrrank = resolvesimplerank(Cstrrank, Cstrtmp);
   Ctimerank = 1:numel(Ctime); Ctimerank = resolvesimplerank(Ctimerank, Ctime);
   [RR,PP] = corrcoef(Cstrrank, Ctimerank); 
   score = RR(1,2);
   if (score ==1) && (neq ==4)
           posprob = 0.0417;
   elseif (score == -1) && (neq == 4)
           negprob = 0.0417;
   elseif (score >= 0)
           posprob = PP(1,2); 
   else
           negprob = PP(1,2);
   end
end
function rankorder = resolvesimplerank(rankorder, tmprank)
%%%%rank the cells: not easy - find equal groups,
uniT = unique(tmprank); ngroup = numel(uniT); 
if (ngroup < numel(tmprank))
   for (i = 1:ngroup)
        iii = find(tmprank == uniT(i));
        if (numel(iii) > 1)
            rankorder(iii) = mean(rankorder(iii))*ones(1, numel(iii));
        end
   end
end

function [score, posprob, negprob] = findseqscoreprob_multiple(Cstr, Ctime, tmpseq, tmprank, rrank, seqparm) %%%%
%nshuffle = 1000;
neq = numel(Cstr); posprob = NaN; negprob = NaN;
%%%%work out str ranks
Cstrrank = zeros(1, neq); for (i = 1:neq) Cstrrank(i) = tmprank(tmpseq == Cstr(i)); end
% shufscore = NaN*ones(1, nshuffle); 
% %   if strcmp(seqparm.rankmode, 'PairRank') 
%        score = findnewrankscore(Cstrrank); %, Cstr, Ctime);
%        for (j = 1:nshuffle)
%            %iii = randperm(neq); Cstrranknow = Cstrrank(iii); Cstrnow = Cstr(iii); shufscore(j) = findnewrankscore(Cstrranknow, Cstrnow, Ctime);
%            shufscore(j) = findnewrankscore(Cstrrank(randperm(neq)));
%        end 
%        posprob = numel(find(shufscore>=score))/nshuffle;
%        negprob = numel(find(shufscore<=score))/nshuffle;
%   elseif strcmp(seqparm.rankmode, 'OrderRank')
       [RR,PP] = corrcoef(Cstrrank, Ctime); 
       score = RR(1,2); 
       if (score >= 0)
           posprob = PP(1,2); 
       else
           negprob = PP(1,2);
       end
       
%       score = DE_ComputeSimpleCrr(Cstrrank, timerank); %Ctime);
%        for (j = 1:nshuffle)
%            iii = randperm(neq); Cstrranknow = Cstrrank(iii); 
%            Cstrnow = Cstr(iii); shufscore(j) = DE_ComputeSimpleCrr(Cstrranknow, timerank); %Ctime);
%        end 
%   end
function score = findnewrankscore(Crank) %, Cstr, Ctime)
nq = numel(Crank); PP = 0; NN = 0; 
for (i = 1:nq-1)
     for (j = i+1:nq)
         %if (Cstr(j) ~= Cstr(i))
             %rawscore = (Ctime(j)-Ctime(i))*(Crank(j)-Crank(i));
             rawscore = Crank(j)-Crank(i);
             if (rawscore>0)
                    PP = PP +1;
             elseif (rawscore < 0)
                    NN = NN + 1;
             %elseif (Crank(j)==Crank(i)) && (abs(j-i)==1) %%%if ctime next each other for zeo-peak pairs 
             %       PP = PP + 1;
             %else
             %       NN = NN + 1;
             end
         %end
     end
end
score = (PP-NN)/(nq*(nq-1)/2);


function [allletter, alltime] = findfreesequences(spikedata, matchletter, PeakTh, timepoint) %%%first, map out entire
nspike = numel(spikedata); peakL = cell(1,nspike); peakT = cell(1,nspike);
for (i = 1:nspike)
    [~, maxind] = FindLocalMinima(spikedata{i}); maxV = spikedata{i}(maxind);
    iii = find( maxV > PeakTh(i) );  peakT{i} = timepoint(maxind(iii)); peakL{i} = matchletter(i)*ones(1,numel(iii));
end
allL = cell2mat(peakL); allT = cell2mat(peakT);
[alltime, iii] = sort(allT); allletter = allL(iii);

function tmpnow = findtemplates(tmp, animaldatenow)
int = zeros(1, numel(tmp));
for (i = 1:numel(tmp))
    if strcmp(tmp(i).parm.animaldate, animaldatenow) int(i) = 1; end
end
tmpnow = tmp(int==1);

function [seq, seqdata] = assignseqoutput(seq, seqdata, seqnow, seqdatanow, nt)
cat = fieldnames(seqnow); 
for (i = 1:numel(cat))
    subf = fieldnames(seqnow.(cat{i}));
    for (j = 1:numel(subf))
        seq.(cat{i}).(subf{j}){nt} = seqnow.(cat{i}).(subf{j});
    end
end
cat = fieldnames(seqdatanow); 
for (i = 1:numel(cat))
    subf = fieldnames(seqdatanow.(cat{i}));
    for (j = 1:numel(subf))
        seqdata.(cat{i}).(subf{j}){nt} = seqdatanow.(cat{i}).(subf{j});
    end
end

function data = getfiringrates_binspikes(spiketime, sigma, timepoint, binsize, sT, eT)
%%%%%%%%%computing rate using conventional bin counts
%%build the gaussian window
% disp(size(spiketime));
% disp(sT);
% disp(eT);
iii = find( (spiketime>=sT) & (spiketime<=eT) ); 
%disp(size(iii));
spiketime = spiketime(iii);
Nsigma = 5;
data = -NaN*ones(numel(timepoint),1); timepoint = timepoint'; %%%%column vector so output count is also clumn.
if (~isempty(timepoint))
   sig = round(sigma/binsize); %sigma now in binsize
   X = -Nsigma*sig:1:Nsigma*sig; %compute both side two-sigma away
   aa = normpdf(X, 0, sig); wind = aa/sum(aa);
   nw = numel(wind); %number of data points in a window: need to get rid of first and last (nw-1)/2 data poiints
   ybin = histc(spiketime, timepoint);
   ybin = conv(ybin, wind);
   data = (1/binsize)*ybin( (nw-1)/2+1 : numel(ybin)-(nw-1)/2 ); %%%final results in Hz
end

function [PeakTh, spikedata, ratemean, ratestd] = normalizefindthreshold(spikedata, Pthre)
nspike = numel(spikedata); PeakTh = 10000*ones(1, nspike); ratestd = NaN*ones(1, nspike); ratemean = NaN*ones(1, nspike);
for (i = 1:nspike)
     if ~iscell(spikedata{i})
         datanow = spikedata{i};
     else
         datanow = [];
         for (j = 1:numel(spikedata{i}))
             datanow = [datanow spikedata{i}{j}];
         end
     end
     sss = std(datanow); mmm = mean(datanow); ratestd(i) = sss; ratemean(i) = mmm;
     if (sss~=0) PeakTh(i) = Pthre; end
     if ~iscell(spikedata{i})
         spikedata{i} = spikedata{i} - mmm;
         if (sss~=0) spikedata{i} = spikedata{i}/sss; end
     else
         for (j = 1:numel(spikedata{i}))
             spikedata{i}{j} = spikedata{i}{j} - mmm;
             if (sss~=0) spikedata{i}{j} = spikedata{i}{j}/sss; end
         end
     end
end

function [score, posprob, negprob] = findseqscoreprob(Cstr, tmpseq, rrank) %%%%this works for straight template: 0123456A, etc.
score = NaN; posprob = NaN; negprob = NaN;
Cstr = 1*Cstr; tmpseq = 1*tmpseq; %%%%now both inputs are in numerical representation
nq = numel(Cstr); PP = 0; NN = 0; tmpind = zeros(1,nq);
for (i = 1:nq) 
    iii = find(tmpseq == Cstr(i));
    if (numel(iii)==1) tmpind(i) = iii; end
end
tmpind = tmpind(tmpind>0); nq = numel(tmpind);
if (numel(unique(Cstr)) > 3) && (nq>3)
       [RR,PP] = corrcoef(tmpind, 1:nq); 
       score = RR(1,2); 
       if (score >= 0)
           posprob = PP(1,2); 
       else
           negprob = PP(1,2);
       end

%    for (i = 1:nq-1)
%        for (j = i+1:nq)
%            matchnow = tmpind(j)-tmpind(i);
%            if (matchnow>0)
%                PP = PP +1;
%            elseif (matchnow<0)
%                NN = NN + 1;
%            end
%        end
%    end
%    score = (PP-NN)/(PP+NN);
%    epsilon = 0.001; posprob = 0; negprob = 0;
%    posindex = find(rrank.ratio{nq} >= score-epsilon);
%    negindex = find(rrank.ratio{nq} <= score+epsilon);
%    if (~isempty(posindex))
%        posprob = sum(rrank.prob{nq}(posindex));
%    end
%    if (~isempty(negindex))
%        negprob = sum(rrank.prob{nq}(negindex));
%    end
end

function [evSess, evsessid] = identifysession(evTime, sessionname, startT, endT)
evSess = []; evsessid = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii}; evsessid = iii;
end

function [spatrate, spatX] = findspatialrate(spiketime, behav, bhdata, finaldirnow, evsessName, evName, evTime, sigma, Nsigma)
spatrate = []; spatX = []; posid = []; evid = [];
if (~isempty(evsessName))
    posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evsessName) );
end
if numel(posid) == 1
       if (isfield(behav.general, 'eventname'))
           evid = find(strcmp(behav.general.eventname{posid}, evName));
       else
           evid = find(strcmp(behav.behavior.eventname{posid}, evName));
       end
end
if (numel(posid)~=1)||(numel(evid)~=1)
       disp(['-------------> 1D field not computed: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName]);
else
       spatX = bhdata.event.Xbin{posid}{evid};
       nlap = numel(evTime.start); spatrate = cell(1, nlap);
       for (jnow = 1:nlap)
            occutimenow = bhdata.event.LapOccuptime{posid}{evid}{jnow}; 
            lappostime{1} = bhdata.event.LapAllPostimestamp{posid}{evid}{jnow}; lapx{1} = bhdata.event.LapAllX{posid}{evid}{jnow};
            evok.start = evTime.start(jnow); evok.ent = evTime.ent(jnow);
            [spatrate{jnow}, ~] = getlinearmap(spiketime, evok, lappostime, lapx, spatX, occutimenow, sigma, Nsigma);
       end 
end
function [D1rate, nact] = getlinearmap(spiketime, evTime, lappostime, lapx, xbin, occutime, sigma, Nsigma)
nx = numel(xbin); D1rate = zeros(nx,1); nlap = numel(lappostime); count = zeros(nx,1); spikex = []; nact = NaN;
if (nx > 1)
    binsize = xbin(2) - xbin(1); sigma = round(sigma/binsize); %%sigma now in number of bins
%%%%%count spikes in xbin lap bu lap
for (i = 1:nlap)
    counti = zeros(nx,1);
    spikenow = sort( spiketime( (spiketime>=evTime.start(i)) & (spiketime<=evTime.ent(i)) ) );
    spikex = NaN*ones(size(spikenow)); [laptimenow, iii] = sort(lappostime{i}); allxnow = lapx{i}(iii);
    lastpoint = 1;  %this only for saving time
    for (j = 1:numel(spikenow))
         for (k = lastpoint:numel(laptimenow)-1) %find corresponding time in position data
             if (laptimenow(k) <= spikenow(j)) && (laptimenow(k+1) > spikenow(j)) 
                 spikex(j) = allxnow(k); lastpoint = k; break; 
             end
         end
    end
    for (j = 1:nx-1) counti(j) = numel(find((spikex>=xbin(j)) & (spikex<xbin(j+1)))); end
    counti(nx) = numel(find(spikex>=xbin(nx)));
    count = count + counti;
end
%%%%%compute firing rate
for (i = 1:nx)
     if (occutime(i) ~= 0) D1rate(i) = count(i)/occutime(i); end
end
nact = numel(find(count>0))/numel(find(occutime>0));
%%now do a 1d smoothing
if (strcmpi(smParm.d1sm, 'yes'))
    %D1rate = OneDSmooth_new(D1rate, occutime, sigma, smParm.d1Nsig);
    XX = [-Nsigma*sigma:1:Nsigma*sigma]; 
    wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
    D1rate = conv(D1rate, wind); D1rate = D1rate((nw-1)/2+1:numel(D1rate)-(nw-1)/2);
end
end

function [othertmpname,relationtype,matchscore, matchposprob, matchnegprob] = findrelationwithothertmps(tmpnow, alltmp, rrank)
allanimdate = cell(1, numel(alltmp));
for (i = 1:numel(alltmp)) allanimdate{i} = alltmp(i).parm.animaldate; end
iii = find(strcmp(allanimdate, tmpnow.parm.animaldate));
okevfile = cell(1, numel(iii));
for (i = 1:numel(iii)) okevfile{i} = alltmp(iii(i)).evtfile; end
ii = find(strcmp(okevfile, tmpnow.evtfile)); iii = iii(setdiff(1:numel(iii), ii));
othertmpname = okevfile(iii); nothertmp = numel(iii);
relationtype = cell(1, nothertmp); matchscore = NaN*ones(1, nothertmp); 
matchposprob = NaN*ones(1, nothertmp); matchnegprob = NaN*ones(1, nothertmp);
for (i = 1:nothertmp)
    evthen = alltmp(iii(i)).evtfile; seqthen = alltmp(iii(i)).tmp.tmpfinalseq{1}; sessevthen = alltmp(iii(i)).parm.sessev;
    relationtype{i} = findrelationtype(evthen, sessevthen, tmpnow.evtfile, tmpnow.parm.sessev);
    [matchscore(i), matchposprob(i), matchnegprob(i)] = findseqscoreprob(seqthen, tmpnow.tmp.tmpfinalseq{1}, rrank);
end
function relationtype = findrelationtype(evthen, sessevthen, evnow, sessevnow)
relationtype = 'unknown';
if strcmp(sessevthen, sessevnow)
   if strcmp(sessevthen, 'sessions')
      nn= min([numel(evthen) numel(evnow)]);
      if strncmpi(evthen, evnow, nn-1)
          relationtype = 'sameTask';
      else
          relationtype = 'diffTask';
      end
   elseif strcmp(sessevthen, 'events')
      [str1, tok1] = strtok(evthen, '_'); [str2, tok2] = strtok(evnow, '_');
      if strcmpi(str1, str2)
          relationtype = 'crossSess';
      elseif strcmpi(tok1, tok2)
          relationtype = 'crossTraj';
      end
   end
end

function [npos, nneg] = findmatchnum(evposprob, evnegprob, evseqstart, evseqend, slevel)
%%%%%%here seq{i,j} corresponds to ep{i,j}
%%%%%found sig matches need to resolve timing overlaps within laps (cross all window size if sliding windows)
npos = 0; nneg = 0;
[mm,nn] = size(evposprob); 
for (i = 1:mm)
    seq = []; posprob = []; negprob = []; seqstart = []; seqend = [];
    for (j = 1:nn)
        posprob =[posprob; evposprob{i,j}]; negprob = [negprob; evnegprob{i,j}];
        seqstart =[seqstart; evseqstart{i,j}]; seqend = [seqend; evseqend{i,j}];
    end
    iii = (find(posprob<slevel));%%%%positive matches
    [newstart, ~, ~, ~] = resolveoverlap(posprob(iii), seqstart(iii), seqend(iii)); npos = npos + numel(newstart);
    iii = find(negprob<slevel);%%%%negative matches
    [newstart, ~, ~, ~] = resolveoverlap(negprob(iii), seqstart(iii), seqend(iii)); nneg = nneg + numel(newstart);
end
function [npos, posP, posmm, nneg, negP, negmm, posSig, negSig] = findtheorysig(evseq, slevel, rrank, eventoption)
%%%%%%here seq{i,j} corresponds to ep{i,j}
%%%%%found sig matches need to resolve timing overlaps within laps (cross all window size if sliding windows)
%%%%%theory significance produced only if eventoption = WithinEvents
posP = NaN; posmm = NaN; negP = NaN; negmm = NaN; npos = 0; nneg = 0; posSig = NaN; negSig = NaN;
[mm,nn] = size(evseq.seq); 
for (i = 1:mm)
    seq = []; posprob = []; negprob = []; mscore = []; seqstart = []; seqend = [];
    for (j = 1:nn)
        seq = [seq; evseq.seq{i,j}]; posprob =[posprob; evseq.posmatchprob{i,j}]; negprob = [negprob; evseq.negmatchprob{i,j}];
        mscore = [mscore; evseq.matchscore{i,j}]; seqstart =[seqstart; evseq.seqstart{i,j}]; seqend = [seqend; evseq.seqend{i,j}];
    end
    iii = (find(posprob<slevel));%%%%positive matches
    [newstart, ~, ~, ~] = resolveoverlap(posprob(iii), seqstart(iii), seqend(iii)); npos = npos + numel(newstart);
    iii = find(negprob<slevel);%%%%negative matches
    [newstart, ~, ~, ~] = resolveoverlap(negprob(iii), seqstart(iii), seqend(iii)); nneg = nneg + numel(newstart);
end
if strcmp(eventoption, 'WithinEvents') %%%not available to 'FreeSequences', dangerous to 'Sliding windows'
   npos = numel(find(posprob<slevel)); nneg = numel(find(negprob<slevel)); %%%%positive/negative matches
   %%%%%expected distribution of sig seq
   [posexpprob, negexpprob] = findexpprob(rrank, slevel);
   seqleng = zeros(1, numel(seq));
   for (i = 1:numel(seq)) seqleng(i) = numel(seq{i}); end
   Npot = zeros(size(posexpprob)); %%number of potential seqs for each seq length
   for (i = 1:numel(Npot)) Npot(i) = numel(find(seqleng == i)); end
   posmm = Npot*posexpprob'; posP = 1-poisscdf(npos-1, posmm); posSig = sqrt(posmm);
   negmm = Npot*negexpprob'; negP = 1-poisscdf(nneg-1, negmm); negSig = sqrt(negmm); %%%%assuming poisson distribution, sigma = sqrt(mean)
end
function [posexpprob, negexpprob] = findexpprob(rrank, slevel)
maxncell = numel(rrank.ratio);
for (i = 1:maxncell)
    nratio(i) = numel(rrank.ratio{i});
end
posexpprob = zeros(1, maxncell); negexpprob = zeros(1, maxncell); %row vector
for (kkk = 2:maxncell)
     ii = numel(rrank.ratio{kkk}); jj = 1;
     for (xxx = ii:-1:1)
          if (sum(rrank.prob{kkk}(xxx:ii)) > slevel)
              jj = xxx+1; break;
          end
     end
     if (jj<=ii) posexpprob(kkk) = sum(rrank.prob{kkk}(jj:ii)); end
end
for (kkk = 2:maxncell)
     ii = numel(rrank.ratio{kkk}); jj = ii;
     for (xxx = 1:ii)
          if (sum(rrank.prob{kkk}(1:xxx)) > slevel)
             jj = xxx-1; break;
          end
     end
     if (jj>=1) negexpprob(kkk) = sum(rrank.prob{kkk}(1:jj)); end
end
function [newstart, newend, newscore, jjj] = resolveoverlap(score, sigstart, sigend)
newstart = []; newend = []; newscore = []; jjj = [];
for (i = 1:numel(sigstart))
    iii = findoverlap(newstart, newend, sigstart(i), sigend(i)); 
    if isempty(iii)
        newstart = [newstart sigstart(i)]; newend = [newend sigend(i)]; newscore = [newscore score(i)];
        jjj = [jjj i];
    else
        if (score(i)>max(newscore(iii)))
            kkk = setdiff( (1:numel(newstart)), iii );
            newstart = newstart(kkk); newend = newend(kkk); newscore = newscore(kkk); jjj = jjj(kkk);
            newstart = [newstart sigstart(i)]; newend = [newend sigend(i)]; newscore = [newscore score(i)];
            jjj = [jjj i];
        end
    end
end
[newstart, kkk] = sort(newstart); newend = newend(kkk); newscore = newscore(kkk); jjj = jjj(kkk);
function ind = findoverlap(newstart, newend, sigstart, sigend)
ind = [];
iii = find( (newstart<sigstart) & (newend>sigstart) ); ind = union(ind, iii);
iii = find( (newstart<sigend) & (newend>sigend) ); ind = union(ind, iii);
iii = find( (newstart>=sigstart) & (newend<=sigend) ); ind = union(ind, iii);

function [posP, posE, negP, negE, posSig, negSig] = findshufflesig(tmpseq, npos, nneg, sessST, sessET, fileletter, evseq, rrank, seqparm, sessev, tmprank) %, alltime)
posP = NaN; posE = NaN; negP = NaN; negE = NaN; posSig = NaN; negSig = NaN;
%%%%%%here seq{i,j} corresponds to ep{i,j}
%%%%first generate shuffles
ok = 1; nshuffle = seqparm.Nshuffle; slevel = seqparm.significancelevel;
if (strcmp(sessev, 'session'))
    eventoption = seqparm.sessionoption; timingoption = seqparm.sessiontimingoption;
else
    eventoption = seqparm.eventoption; timingoption = seqparm.eventtimingoption;
end
rankorder = 1:numel(tmprank); rankorder = resolvesimplerank(rankorder, tmprank);
if strcmp(seqparm.shufflemode, 'GlobalID')
   [~, ~, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_gloabalID(tmpseq, fileletter, evseq, seqparm, rrank, rankorder, eventoption);%, alltime);
elseif strcmp(seqparm.shufflemode, 'LocalID')
   [~, ~, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_localID(tmpseq, fileletter, evseq, seqparm, rrank, rankorder, eventoption);%, alltime);
elseif strcmp(seqparm.shufflemode, 'CircleSlide')
   [~, ~, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_circleslide(tmpseq, fileletter, evseq, seqparm, ...
       sessST, sessET, rrank, tmprank, eventoption, timingoption);
else
    ok = 0; disp('----------> can not generate shuffled data or shuffle mode not available; shuffle significance not computed');
end
if ok
   shufnpos = zeros(1,nshuffle); shufnneg = zeros(1,nshuffle);
   for (k = 1:nshuffle)
     [shufnpos(k), shufnneg(k)] = findmatchnum(posmatchprob{k}, negmatchprob{k}, seqstart{k}, seqend{k}, slevel); 
   end
   posP = numel(find(shufnpos>=npos))/nshuffle; negP = numel(find(shufnneg>=nneg))/nshuffle;
   posE = mean(shufnpos); negE = mean(shufnneg); posSig = std(shufnpos); negSig = std(shufnneg);
end
function [seq, matchscore, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_gloabalID(tmpseq, fileletter, evseq, seqparm, rrank, rankorder, eventoption) %alltime
nshuffle = seqparm.Nshuffle; %%%%this version is different from the DE version, but doing the same thing: shuffle template instead of sequence
seq = cell(1,nshuffle); matchscore = cell(1, nshuffle); posmatchprob = cell(1,nshuffle); negmatchprob = cell(1, nshuffle);
seqstart = cell(1, nshuffle); seqend = cell(1, nshuffle); 
%%%%%%%%%%%%%%%%%still evseq.seq{i,j} corresponds to ep{i,j}
%%%%%%%%%%for global ID shuffling: just need to randomly swap cell IDs and then replace original sequences 
sequence = evseq.seq; ncell = numel(tmpseq); [mm,nn] = size(sequence);
for (i = 1:nshuffle)
   %tt = cputime;
    %%%%%these variable do not change
    seqstart{i} = evseq.seqstart; seqend{i} = evseq.seqend; 
    seq{i} = cell(mm,nn); matchscore{i} = cell(mm,nn); posmatchprob{i} = cell(mm,nn); negmatchprob{i} = cell(mm,nn); 
    tmpseqnow = tmpseq(randperm(ncell)); %mletternow = fileletter(tmpind(IX));
    for (ii = 1:mm)
        for (jj = 1:nn)
            for (j = 1:numel(sequence{ii,jj}))
                [matchscore{i}{ii,jj}(j), posmatchprob{i}{ii,jj}(j), negmatchprob{i}{ii,jj}(j)] = findseqscoreprob_simplerankorder(...
                   sequence{ii,jj}{j}, evseq.seqtime{ii,jj}{j}, tmpseqnow, rankorder); %%%%
                % [matchscore{i}{ii,jj}(j), posmatchprob{i}{ii,jj}(j), negmatchprob{i}{ii,jj}(j)] = findseqscoreprob_multiple(...
                %           sequence{ii,jj}{j}, evseq.seqtime{ii,jj}{j}, tmpseqnow, tmprank, rrank, seqparm);
            end
        end
    end
    %disp(['----------> single shuffle time: ', num2str(cputime-tt)]);
end
function [seq, matchscore, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_localID(tmpseq, fileletter, evseq, seqparm, rrank, rankorder, eventoption) %, alltime)
nshuffle = seqparm.Nshuffle;
seq = cell(1,nshuffle); matchscore = cell(1, nshuffle); posmatchprob = cell(1,nshuffle); negmatchprob = cell(1, nshuffle);
seqstart = cell(1, nshuffle); seqend = cell(1, nshuffle); 
%%%%%%%%%%%%%%%%%still evseq.seq{i,j} corresponds to ep{i,j}
%%%%%%%%%%for global ID shuffling: just need to randomly swap cell IDs and then replace original sequences 
sequence = evseq.seq; ncell = numel(tmpseq); [mm,nn] = size(sequence);
tmpind = zeros(1, ncell); for (i = 1:ncell) tmpind(i) = find(fileletter==tmpseq(i)); end
for (i = 1:nshuffle)
    %%%%%these variable do not change
    seqstart{i} = evseq.seqstart; seqend{i} = evseq.seqend; 
    seq{i} = cell(mm,nn); matchscore{i} = cell(mm,nn); posmatchprob{i} = cell(mm,nn); negmatchprob{i} = cell(mm,nn); 
    for (ii = 1:mm)
        for (jj = 1:nn)
            for (j = 1:numel(sequence{ii,jj}))
                 IX = randperm(ncell); mletternow = fileletter(tmpind(IX));
                 seqnow = sequence{ii,jj}{j}; nletter = numel(seqnow); ind = zeros(1, nletter);
                 for (k = 1:nletter)
                      ind(k) = find(tmpseq == seqnow(k));
                 end
                 if isempty(find(ind == 0))
                    seq{i}{ii,jj}{j} = mletternow(ind);
                    Ctime = evseq.seqtime{ii,jj}{j}; %alltime( (alltime>=seqstart{i}{ii,jj}(j)) & (alltime<=seqend{i}{ii,jj}(j)) );
                    [matchscore{i}{ii,jj}(j), posmatchprob{i}{ii,jj}(j), negmatchprob{i}{ii,jj}(j)] = findseqscoreprob_simplerankorder(...
                            sequence{ii,jj}{j}, Ctime, tmpseq, rankorder);
                    %[matchscore{i}{ii,jj}(j), posmatchprob{i}{ii,jj}(j), negmatchprob{i}{ii,jj}(j)] = findseqscoreprob_multiple(...
                    %       seq{i}{ii,jj}{j}, Ctime, tmpseq, tmprank); %, rrank, seqparm);
                 end
            end
        end
    end
end
function [seq, matchscore, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_circleslide(tmpseq, fileletter, evseq, seqparm, sessST, sessET, rrank, tmprank, eventoption, timingoption)
nshuffle = seqparm.Nshuffle; %%%%%%%now this need to change to just slide peaktimes against each other -- not done yet
seq = cell(1,nshuffle); matchscore = cell(1, nshuffle); posmatchprob = cell(1,nshuffle); negmatchprob = cell(1, nshuffle);
seqstart = cell(1, nshuffle); seqend = cell(1, nshuffle); seqmarker = cell(1, nshuffle);
%%%%%%%%%%for circular slide shuffling: randomly slide and circularly wrap every rate curves (spike trains) and re-compute
matchletter = fileletter; PeakTh = evseq.peakth;
spikedata = evseq.seqdata; ep = evseq.ep; timepoint = evseq.timepoint;   
[mm, nn] = size(ep);
for (i = 1:nshuffle) 
    datanow = circleslideratecurve(spikedata, timingoption, sessST, sessET); %%%circle slide: already taking into account the spatial {nlap} case
    %%%%re-compute sequences
    seqstart{i} = cell(mm,nn); seqend{i} = cell(mm,nn); 
    seq{i} = cell(mm,nn); matchscore{i} = cell(mm,nn); posmatchprob{i} = cell(mm,nn); negmatchprob{i} = cell(mm,nn); 
    if ~strcmp(eventoption, 'FreeSequences')
       for (ii = 1:mm)
       for (jj = 1:nn)
             [seq{i}{ii,jj}, seqstart{i}{ii,jj}, seqend{i}{ii,jj}, ~, matchscore{i}{ii,jj}, posmatchprob{i}{ii,jj}, negmatchprob{i}{ii,jj}] = ...
                          findallsequence(datanow, matchletter, timingoption, ep{ii,jj}, PeakTh, timepoint, rrank, tmpseq);
       end
       end    
    else %%%%for free sequencing
       [allletter, alltime] = findfreesequences(datanow, matchletter, PeakTh, timepoint); %%%first, map out entire
       %%then break down the sequence to shoter ones with each event
        for (ii = 1:mm) %%%original event index
        for (jj = 1:nn)
            [seq{i}{ii,jj}, seqstart{i}{ii,jj}, seqend{i}{ii,jj}, ~, matchscore{i}{ii,jj}, posmatchprob{i}{ii,jj}, negmatchprob{i}{ii,jj}, ~] = ...
               breakallsequences_nonselect(allletter, alltime, ep{ii,jj}, tmprank, rrank, tmpseq, seqparm, sessType);
%            [seq{i}{ii,jj}, seqstart{i}{ii,jj}, seqend{i}{ii,jj}, ~, matchscore{i}{ii,jj},...
%                posmatchprob{i}{ii,jj}, negmatchprob{i}{ii,jj}] = ...
%                breakallsequences(allletter, alltime, ep{ii,jj}, rrank, tmpseq);
        end
        end
    end
end
function datanow = circleslideratecurve(spikedata, timingoption, sessstart, sessend)
datanow = spikedata; %%%a single cell data: long vectors for non-spatial-rate data
if ~isempty(strfind(timingoption, '_time')) %%%if peak time%%%if a temporal rate curve
    %%%%%%%%%%need to get rid of trailing (pre and post) zeros, otherwise zeros dominates since the session start/end data are lost and use 0 as start time instead lost
    %%%%%%%%%%%%%%%%%evt start/end data have to be available, otherwise the sliding may slide into other sessions*****
    nbin = numel(spikedata); IX = mod((1:nbin)+ ceil(rand*nbin), nbin) + 1; datanow = spikedata(IX);
elseif ~isempty(strfind(timingoption, '_space')) %%%if peak time%%%if a spatial rate curve {nlap}
    nev = numel(spikedata);
    for (i = 1:nev)
        nbin = numel(spikedata{i}); IX = mod((1:nbin)+ ceil(rand*nbin), nbin) + 1; datanow{i} = spikedata{i}(IX);
    end
else  %%%%if a spike train: circle slide inside the trains
    tintv = sessend-sessstart; tmin = sessstart; tmax = sessend; %%%%this is session timing even for event shuffles since same spike trains are used.
    datanow = spikedata-tmin; datanow = datanow + rand*(tmax-tmin);
    datanow = mod(datanow, tintv) + tmin;
end





