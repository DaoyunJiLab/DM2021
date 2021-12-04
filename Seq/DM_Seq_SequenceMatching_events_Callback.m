function DM_Seq_SequenceMatching_events_Callback
%%Sequence analysis on a database (.spikedb), require a matching .tmp file
%%%%%%%%%%%%%%%ONLY DO IT WITHIN EVENTS SPECIFIED BY A KEYWORD%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%work on selected cells (groups), group cells into days, and work on each day
%%%%%%%% find if templates exist and then work on each template (if no template found for a date, use default/balnk template - natural order)
%%%%%%%%%%%%% sequence through every sessions and every events
%%%%%%%%%%%%%%%%%%%for each session/event: output sequences and result variables
%%%%event parsing options: free sequencing; single sequence per event; sliding windows
%%%%sequence timing options: first spike, spike center, rate main peak; rate first peak, spatial rate first peak, spatial rate main peak 
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
    seqparm.rankmode = 'OrderRank/rrank'; rrank = []; %%%not applicable at this time: %%%%%%%%Free sequences use rank order correlation; will be automatically changed later
                                                                                %%%%%%%%WithinEvents/SlidingWindows use pairrank similarity
    input = inputdlg({'Event keyword';'Event type'}, 'Event slelection parms', 2, {'Track1_leftright'; 'run'}); 
    if (~isempty(input))
        seqparm.evkeyword = input{1}; seqparm.evkeytype = input{2};
        [MCroot, MCname, DAname, DEname, ANname, MAname] = CurrentVersion;
        rfile = fullfile(MCroot, 'DataExplorer', 'Sequence', 'rrankprob.mat'); 
        if (exist(rfile, 'file') == 2)
               S = load(rfile); rrank = S.rrank; S = [];
        else
               ok = 0; disp('-----> ranking probability file for selected rank mode not found; aborted');
        end
    else
        ok = 0;
    end
end
if ok
    oTT = questdlg('Open a template file?', 'Template input', 'Yes');
    if (strcmp(oTT, 'Yes')) %
        cpathname = fullfile(cd, '*.tmp');
        [fname, pname] = uigetfile(cpathname, 'Select a template file to open:');
        if fname ~= 0
           tmpfile = fullfile(pname, fname);
           disp(['-----> template file: ', tmpfile]); seqparm.tmpfile = tmpfile;
           S = load(tmpfile, '-mat'); tmp = S.pinfo; S = [];
           disp(['--------> number of templates found: ', num2str(numel(tmp))]);
        else
            ok = 0;
        end
    elseif (strcmp(oTT, 'No')) 
        seqparm.tmpfile = []; tmp = [];
    else
        ok = 0;
    end
end
if ok
   seqparm.eventoption = 'FreeSequences';
   evsss = {'FreeSequences'; 'WithinEvents'; 'SlidingWindows'};
   [sss, ok] = listdlg('ListString', evsss, 'PromptString', 'Event sequencing option');
   if ok seqparm.eventoption = evsss{sss};  end
end
if ok
    seqparm.eventtimingoption = 'RatePeak_time'; 
    seqparm.timeres = 0; seqparm.spatres = 0; seqparm.timesigma = 0; seqparm.spatsigma = 0;
    seqparm.Nsigma = 5; seqparm.timePthre = 0; seqparm.spatPthre = 0;
    if strcmp(seqparm.eventoption, 'FreeSequences')
      sss = {'RatePeak_time'; 'RatePeak_space'}; seqparm.rankmode = 'OrderRank';
    else
      sss = {'FirstSpike'; 'SpikeTrainCenter'; 'RatePeak_time'; 'RatePeak_space'; 'Rate1stPeak_time'; 'Rate1stPeak_space'};
      seqparm.rankmode = 'PairRank';
    end
    [S, ok] = listdlg('ListString', sss, 'PromptString', 'Event timing option');
    if ok
          seqparm.eventtimingoption = sss{S};
          if ~isempty(strfind(seqparm.eventtimingoption, '_time')) %%%%if use peak times of temporal rate curves
              input = inputdlg({'Temporal rate resolution (s)'; 'Smoothing option'; 'Smoothing sigma (s)'; 'Enter peak threshld (std)'}, ...
                  'Temporal timing parameter Input', 3, {'0.005'; 'Yes'; '0.5'; '0.05'}); 
              if (~isempty(input))
                 seqparm.timeres = str2num(input{1}); seqparm.timesmoothopt = input{2}; seqparm.timesigma = str2num(input{3}); seqparm.timePthre = str2num(input{4});
              else
                 ok = 0;
              end
          end
          if ~isempty(strfind(seqparm.eventtimingoption, '_space')) %%%%if use peak locations of spatial rate curves
              if plotparm.linkbehav == 1
                 behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
              else
                 ok = 0; disp('-----> No linked behavdb; aborted');
              end
              input = inputdlg({'Spatial rate resolution (cm)'; 'Spatial smoothing option'; 'Smoothing sigma (cm)'; 'Peak threshld (std)'}, ...
                  'Spatial timing parameter Input', 4, {'4'; 'Yes'; '10'; '0.05'}); 
              if (~isempty(input))
                 seqparm.spatres = str2num(input{1}); seqparm.spatsmoothopt = input{2}; seqparm.spatsigma = str2num(input{3}); seqparm.spatPthre = str2num(input{4});
              else
                 ok = 0;
              end
          end
    end
end
if ok
      if strcmp(seqparm.eventoption, 'SlidingWindows')%%events: use sliding windows
          if (~strcmp(seqparm.evtimingoption, 'RatePeak_space'))&&(~strcmp(seqparm.evtimingoption, 'Rate1stPeak_space'))
              input = inputdlg({'Event min/max window sizes (s)'; 'Window size step (s)'; 'Window shifting time'},...
                  'Sliding window parameters', 3, {'1 20'; '1'; '1'}); 
              if (~isempty(input))
                    aa = str2num(input{1}); bb = str2num(input{2}); windowsize = aa(1):bb:aa(numel(aa)); %various wind sizes
                    seqparm.eventminwindowsize = aa(1); seqparm.eventmaxwindowsize = aa(numel(aa)); 
                    seqparm.eventwindowsizestep = bb;
                    seqparm.eventwindowshifttime = str2num(input{3}); %%%window shiftting time
              else
                    ok = 0;
              end
          else
              input = inputdlg({'Event min/max window sizes (cm)'; 'Window size step (cm)'; 'Window shifting distance'},...
                  'Sliding window parameters', 3, {'10 100'; '10'; '5'}); 
              if (~isempty(input))
                   aa = str2num(input{1}); bb = str2num(input{2}); windowsize = aa(1):bb:aa(numel(aa)); %various wind sizes
                   seqparm.eventminwindowsize = aa(1); seqparm.eventmaxwindowsize = aa(numel(aa)); 
                   seqparm.eventwindowsizestep = bb;
                   seqparm.eventwindowshifttime = str2num(input{3}); %%%window shiftting time
              else
                   ok = 0;
              end
           end
   end
end
if ok && strcmp(seqparm.eventoption, 'FreeSequences')
       input = inputdlg({'Max sequence length (x tmplate length)'; 'Max sequence duration (s)'},...
                  'Sequence leng options', 2, {'3'; '60'}); 
       if (~isempty(input))
               seqparm.maxtmpleng = str2num(input{1}); seqparm.maxtime = str2num(input{2}); %%%parameters used to break down long sequences
       else
               ok = 0;
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
    seqparm.minactivecell = 4;
    seqparm.Nshuffle = 200; seqparm.Sigshuffle = 1000; seqparm.shufflemode = 'GlobalID'; seqparm.significancelevel = 0.05;
    [seq, seqdata] = computeseqdb(tmp, pinfo, data, behav, bhdata, cellind, seqparm, rrank);  
    if (~isempty(seq))
       %%%%%group data
       if (~isempty(seq)) seq.work = struct([]); end
       seqdata.grouplist.groupname{1} = 'List0'; seqdata.grouplist.groupindex{1} = 1:numel(seq.general.tmpID);
       seqdata.grouplist.grouptype{1} = 'Manual'; seqdata.grouplist.groupcrit{1} = []; seqdata.grouplist.groupparents{1} = []; 
       seqdata.parentfile = currentfilename; 
       %%%save a stripped version
       pinfo = seq; 
       seqfilestripped = strrep(seqfile, '.seqdb', '_stripped.seqdb');
       data = []; data.grouplist = seqdata.grouplist; 
       save(seqfilestripped, 'pinfo', 'data', '-mat');
       %%%save a full version
       data = seqdata; save(seqfile, 'pinfo', 'data', '-mat');
    end
end
disp('************************');

function [seq, seqdata] = computeseqdb(tmp, pinfo, data, behav, bhdata, cellind, seqparm, rrank)
seqdata = []; nt = 0; seq = [];
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
    if (numel(cellparmfile)<4)
        disp('----------> too few cells; aborted');
    else
        animaldatenow = strcat(pinfo.general.animalname{neuronid(1)}, '_', pinfo.general.datedir{neuronid(1)});
        tmpnow = [];
        if isempty(tmp) %%%use default, blank template
           if (numel(cellparmfile)>62) %%%
               disp('----------> number of cells exceeds template capacity; try to be more selective; aborted');
           else
               cellletter = DE_Seq_AssignLetter(numel(cellparmfile), 'allletters'); 
               tmpnow.tmp.tmpfileletter = cellletter; 
               tmpnow.tmp.tmpfilename = cellparmfile; tmpnow.tmp.tmpfinalseq{1} = cellletter;
           end
        else
            tmpnow = findtemplates(tmp, animaldatenow);
        end    
        for (j = 1:numel(tmpnow)) %%%for each template
            disp(['----------> template now: ', tmpnow(j).tmp.tmpfinalseq{1}, '; ', tmpnow(j).evtfile]);
            nt = nt + 1; %%%entry of the seqdb items
            [sessseq, evseq, tmpnow(j), evName, evType, evT] = sequencedata(pinfo, data, behav, bhdata, neuronid, tmpnow(j), seqparm, rrank, fdirnow);
            tmprank = tmpnow(j).tmp.tmpfinalrank{1}; 
            [seqnow, seqdatanow] = computeseqproperties(sessseq, evseq, seqparm, pinfo, data, neuronid(1), tmpnow(j), tmpnow, rrank, tmprank, evName, evType, evT);
            [seq, seqdata] = assignseqoutput(seq, seqdata, seqnow, seqdatanow, nt);
        end
    end
end

function [sessseq, evseq, tmpnow, evName, evType, evT] = sequencedata(pinfo, data, behav, bhdata, neuronid, tmpnow, seqparm, rrank, finaldirnow)
%%%1. select out cells; 2. if temporal/spatial peaks, work out session smoothed rate curves
%%%3. for each event, sessions, workout the ep scheme and sequence accordingly
%%%%%%%%output: seq{ii,jj} for each event ep{ii,jj}
sessseq = []; evseq = []; 
%%%%%first select which template sequence to use
cellletter = tmpnow.tmp.tmpfileletter; cellfilename = tmpnow.tmp.tmpfilename; 
alltmpseq = tmpnow.tmp.tmpfinalseq;
finalfiterr = []; if isfield(tmpnow.tmp, 'finalfiterr') finalfiterr = tmpnow.tmp.finalfiterr; end %%%select the best fit template to use
minfiterr = NaN;
if ~isempty(finalfiterr)
   [minfiterr, sleind] = min(finalfiterr);
else
    ntmp = numel(alltmpseq); sleind = ceil(rand*ntmp); if (sleind ==0) sleind = 1; end 
end
tmpseq = alltmpseq{sleind}; tmprank = 1:numel(tmpseq);
if isfield(tmpnow.tmp, 'finalrank') tmprank = tmpnow.tmp.finalrank{sleind}; end 
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
nlet = numel(iii); cellletter = cellletter(iii); cellfilename = cellfilename(iii); indok = indok(iii); spikedata = cell(1, nlet);
for (i = 1:nlet) spikedata{i} = data.spike.spiketime{neuronid(indok(i))}; end %%%%This is the rawdata to analyze
tmpisok = zeros(1, numel(tmpseq));
for (i = 1:numel(tmpseq))
    ti = find(cellletter == tmpseq(i));
    if (numel(ti) == 1) tmpisok(i) = 1; end
end
kkk = find(tmpisok==1); tmpseq =tmpseq(kkk); tmprank = tmprank(kkk);
tmpnow.tmp.tmpfileletter = cellletter; tmpnow.tmp.tmpfilename = cellfilename; 
tmpnow.tmp.tmpfinalseq = []; tmpnow.tmp.tmpfinalrank = []; tmpnow.tmp.finalfiterr = [];
tmpnow.tmp.tmpfinalseq{1} = tmpseq; tmpnow.tmp.tmpfinalrank{1} = tmprank; tmpnow.tmp.finalfiterr{1} = minfiterr;
%%%%%%%working on session sequences
sessName = pinfo.general.sessionname{neuronid(1)}; 
sessST = pinfo.general.sessionstartT{neuronid(1)}; sessET = pinfo.general.sessionendT{neuronid(1)}; nsess = numel(sessName);
ratedata = cell(1, nsess); timepoint = cell(1, nsess); timePeakTh = 10000*ones(1,nlet); 
allsessletter = cell(1, nsess); allsesstime = cell(1, nsess);
sessseq = cell(1, nsess);
for (i = 1:nsess) %%%for session data the only available timingoption is 'RatePeak_time'
    seq = []; seqstart = []; seqend = []; seqmarker = []; matchscore = []; posmatchprob = []; negmatchprob = []; seqtime = []; ratemean = []; ratestd = []; ep = [];
    if (~isempty(strfind(seqparm.eventtimingoption, 'Rate'))) %%do session sequencing only if eventtimeoption is "FreeSequencing"
       timepoint{i} = (sessST(i):seqparm.timeres:sessET(i)); ratedata{i}=cell(1, nlet);
       for (j = 1:nlet)
            ratedata{i}{j} = getfiringrates_binspikes(spikedata{j}, seqparm.timesmoothopt, seqparm.timesigma, ...
                timepoint{i}, seqparm.timeres, sessST(i), sessET(i));
       end
       [timePeakTh, ratedata{i}, ratemean, ratestd] = normalizefindthreshold(ratedata{i}, seqparm.timePthre);
        %%%%work out the ep structure (sliding windows)
       ep = []; nep = 1;
       for (ti = 1:nep)
             ep{ti,1}.start = sessST(i); ep{ti, 1}.ent = sessET(i); ep{ti,1}.marker = sessName{i};
       end
    end
    if strcmp(seqparm.eventoption, 'FreeSequences') %%do session sequencing only if eventtimeoption is "FreeSequencing"
       %%%%work out the sequences
       [mm,nn] = size(ep);
       [allsessletter{i}, allsesstime{i}] = findfreesequences(ratedata{i}, cellletter, timePeakTh, timepoint{i}); %%%first, map out entire
       %%then break down the sequence to shoter ones with each event
        for (ii = 1:mm) %%%original event index
        for (jj = 1:nn)
            [seq{ii,jj}, seqstart{ii,jj}, seqend{ii,jj}, seqmarker{ii,jj}, matchscore{ii,jj}, posmatchprob{ii,jj}, negmatchprob{ii,jj}, seqtime{ii,jj}] = ...
               breakallsequences_nonselect(allsessletter{i}, allsesstime{i}, ep{ii,jj}, tmprank, rrank, tmpseq, seqparm);
%            [seq{ii,jj}, seqstart{ii,jj}, seqend{ii,jj}, seqmarker{ii,jj}, matchscore{ii,jj},...
%                posmatchprob{ii,jj}, negmatchprob{ii,jj}] = ...
%                breakallsequences(allsessletter{i}, allsesstime{i}, ep{ii,jj}, rrank, tmpseq);
        end
        end
    end
    sessseq{i}.seqdata = ratedata{i}; sessseq{i}.timepoint = timepoint{i}; 
    sessseq{i}.timepeakth = timePeakTh; sessseq{i}.ratestd = ratestd; sessseq{i}.ratemean = ratemean;
    sessseq{i}.ep = ep; sessseq{i}.allletter = allsessletter{i}; sessseq{i}.alltime = allsesstime{i}; 
    sessseq{i}.seq = seq; sessseq{i}.seqstart = seqstart; sessseq{i}.seqend = seqend; sessseq{i}.seqmarker = seqmarker; sessseq{i}.seqtime = seqtime;
    sessseq{i}.matchscore = matchscore; sessseq{i}.posmatchprob = posmatchprob; sessseq{i}.negmatchprob = negmatchprob;
end
%%%%%%%working on sevent sequences
[evName, evType, evT] = filterevents(pinfo, data, seqparm.evkeyword, seqparm.evkeytype, neuronid);
nev = numel(evName); 
evsessid = zeros(1, nev); evsessName = cell(1, nev); allevletter = cell(1, nev); allevtime = cell(1, nev);
evseq = cell(1, nev);
for (i = 1:nev) %%%for ev data the multiple available timingoption
    seq = []; seqstart = []; seqend = []; seqmarker = []; matchscore = []; posmatchprob = []; negmatchprob = []; seqtime = [];
    [evsessName{i}, evsessid(i)] = identifysession(evT{i}, sessName, sessST, sessET);%%%%work out which session the event is
    if isempty(strfind(seqparm.eventtimingoption, 'space'))
        evT{i}.start = evT{i}.start - seqparm.setbacktime; evT{i}.ent = evT{i}.ent + seqparm.setbacktime; 
    end
    timePeakTh = sessseq{evsessid(i)}.timepeakth;
    %%%%compute spatial firing rate if in spatial mode; also need to recompute timepoint and ep{ii,jj}
    spatX = []; spatpoint = []; spatrate = cell(1, nlet); spatPeakTh = NaN*ones(1, nlet); spatratestd = NaN*ones(1, nlet); spatratemean = NaN*ones(1, nlet);
    if ~isempty(strfind(seqparm.eventtimingoption, 'space'))
        smopt = seqparm.spatsmoothopt;
        sigmanow = seqparm.spatsignma; %=pinfo.parm.fSmooth1DSigma(neuronid(1)); 
        Nsigmanow = seqparm.Nsigma; %=pinfo.fSmooth1DNsigma(neuronid(1));
        if (strcmp(evType{i}, 'run'))
           for (j = 1:nlet)
            [spatrate{j}, spatX] = findspatialrate(spikedata{j}, behav, bhdata, finaldirnow, evsessName{i}, evName{i}, evT{i}, smopt, sigmanow, Nsigmanow);
           end
           spatpoint = min(spatX):seqparm.spatres:max(spatX); [spatPeakTh, spatrate, spatratemean, spatratestd] = normalizefindthreshold(spatrate, seqparm.spatPthre);
        end
    end    
    %%%%work out the ep structure (sliding windows)
    evTime = evT{i}; ep = []; nep = numel(evTime.start);
    if ~isempty(strfind(seqparm.eventtimingoption, 'space')) %%%if spatial mode, evt changed to spatial windows
        if (~isempty(spatX))
            evTime.start = min(spatX)*ones(nep,1); evTime.ent = max(spatX)*ones(nep,1); 
        end
    end
    if strcmp(seqparm.eventoption, 'SlidingWindows') %%if sliding windows
        ep = []; mm = numel(evTime.start); nn = numel(seqparm.eventwindowsize);  ep = cell(mm,nn);
        for (ii = 1:mm)
             winshift = evTime.start(ii):seqparm.eventwindowshifttime:evTime.ent(ii); pp = numel(winshift);
             for (jj = 1:nn)
                   ep{ii,jj}.start = winshift - seqparm.eventwindowsize(jj)/2;
                   ep{ii,jj}.ent = winshift + seqparm.eventwindowsize(jj)/2;
                   for (k = 1:pp) ep{ii,jj}.marker{k} = evTime.marker{ii}; end
             end
        end
    else
        for (ti = 1:nep)
             ep{ti,1}.start = evTime.start(ti); ep{ti, 1}.ent = evTime.ent(ti); ep{ti,1}.marker = evTime.marker{ti};
        end
    end
    %%%%work on the sequences: all different options
    [mm,nn] = size(ep);
    if strcmp(seqparm.eventoption, 'FreeSequences') && strcmp(seqparm.eventtimingoption, 'RatePeak_time')
        alltimenow = allsesstime{evsessid(i)};
        %%then break down the sequence to shoter ones with each event
        for (ii = 1:mm) %%%original event index
            ik = find( (alltimenow>=evTime.start(ii)) & (alltimenow<=evTime.ent(ii)) ); 
            allevletter{i}{ii} = allsessletter{evsessid(i)}(ik); allevtime{i}{ii} = alltimenow(ik);
        for (jj = 1:nn)
            [seq{ii,jj}, seqstart{ii,jj}, seqend{ii,jj}, seqmarker{ii,jj}, matchscore{ii,jj}, posmatchprob{ii,jj}, negmatchprob{ii,jj}, seqtime{ii,jj}] = ...
               breakfromsessionallsequences(sessseq{evsessid(i)}, ep{ii,jj});
        end
        end
        datanow = ratedata{evsessid(i)}; timepointnow = timepoint{evsessid(i)};
    elseif strcmp(seqparm.eventoption, 'FreeSequences') && strcmp(seqparm.eventtimingoption, 'RatePeak_space')
        for (ii = 1:mm) %%%original event index
           [allevletter{i}{ii}, allevtime{i}{ii}] = findfreesequences_space(spatrate, cellletter, spatPeakTh, spatpoint, ii); %%%first, map out entire lap
        for (jj = 1:nn)
            [seq{ii,jj}, seqstart{ii,jj}, seqend{ii,jj}, seqmarker{ii,jj}, matchscore{ii,jj}, posmatchprob{ii,jj}, negmatchprob{ii,jj}, seqtime{ii,jj}] = ...
               breakallsequences_nonselect(allevletter{i}{ii}, allevtime{i}{ii}, ep{ii,jj}, tmprank, rrank, tmpseq, seqparm);
        end
        end
        datanow = spatrate; timepointnow = spatpoint;
    else %%%%this works for both SLidingWindows and WinthinEvents
        if (~isempty(strfind(seqparm.eventtimingoption, 'space')))
            datanow = spatrate; timepointnow = spatpoint;
        elseif (~isempty(strfind(seqparm.eventtimingoption, 'time')))
            datanow = ratedata{evsessid(i)}; timepointnow = timepoint{evsessid(i)};
        else %%%use spike train
            datanow = cell(1, nlet); timepointnow = [];
            for (j = 1:nlet) [datanow{j}, ~] = SpikeEventFilter(spikedata{j}, evT{i}); end
        end
        for (ii = 1:mm) %%%original event index
        for (jj = 1:nn)
            [seq{ii,jj}, seqstart{ii,jj}, seqend{ii,jj}, seqmarker{ii,jj}, matchscore{ii,jj}, ...
                posmatchprob{ii,jj}, negmatchprob{ii,jj}] = ...
                findallsequence(datanow, cellletter, seqparm.eventtimingoption, ...
                ep{ii,jj}, timePeakTh, timepointnow, rrank, tmpseq);
        end
        end
    end
    nactivecell = cell(mm,nn);
    for (ii = 1:mm) %%%original event index
    for (jj = 1:nn)
        for (ttni = 1:numel(seq{ii,jj}))
            nactivecell{ii,jj}(ttni) = numel(seq{ii,jj}{ttni});
        end
    end
    end
    evseq{i}.nactivecell = nactivecell;
    evseq{i}.seqdata = datanow; evseq{i}.timepoint = timepointnow; 
    evseq{i}.timepeakth = timePeakTh; evseq{i}.ratestd = ratestd; evseq{i}.spatratemean = ratemean;
    evseq{i}.spatpeakth = spatPeakTh; evseq{i}.spatratestd = spatratestd; evseq{i}.spatratemean = spatratemean;
    evseq{i}.ep = ep; evseq{i}.allletter = allevletter{i}; evseq{i}.alltime = allevtime{i};
    evseq{i}.seq = seq; evseq{i}.seqstart = seqstart; evseq{i}.seqend = seqend; evseq{i}.seqmarker = seqmarker; evseq{i}.seqtime = seqtime;
    evseq{i}.matchscore = matchscore; evseq{i}.posmatchprob = posmatchprob; evseq{i}.negmatchprob = negmatchprob;
end

function [seq, seqdata] = computeseqproperties(sessseq, evseq, seqparm, pinfo, data, neuronid, tmpnow, alltmp, rrank, tmprank, evName, evType, evTimes)
%%%%now work on what sequence properties to compute: 1-template properties; 2-sequence properties; 3-significance
%%%%%%%use default parameters, can be recomputed later
%%%%%general variables
seq.general.tmpID = strcat(tmpnow.parm.animaldate, '_', tmpnow.evtfile); 
seq.general.recarea = pinfo.general.recarea{neuronid}; seq.general.finaldir = pinfo.general.finaldir{neuronid};
seq.general.datedir = pinfo.general.datedir{neuronid}; seq.general.animalname = pinfo.general.animalname{neuronid};
seq.general.sessionname = pinfo.general.sessionname{neuronid}; seq.general.sessionstartT = pinfo.general.sessionstartT{neuronid};
seq.general.sessionendT = pinfo.general.sessionendT{neuronid}; seq.general.sessionlength = pinfo.general.sessionlength{neuronid};
seq.general.eventname = evName; seq.general.genotype = pinfo.general.genotype{neuronid};
seq.general.sex = pinfo.general.sex{neuronid}; seq.general.age = pinfo.general.age{neuronid};
%%%%parameters
seq.parm = seqparm; seq.parm = rmfield(seq.parm, 'tmpfile'); seq.parm.seqtype = 'SeqMatching';
seq.parm.sessType = pinfo.parm.sessType{neuronid}; seq.parm.eventType = evType;
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
if (~strcmp(tmpnow.tmp.tmptype, 'pfratetmp'))
   seq.tmp.pairscore = tmpnow.tmp.tmpfinalscore/numel(tmpnow.tmp.tmpfinalorderpair); seq.tmp.orderscore = tmpnow.tmp.tmpfinalscore/(nlet*(nlet-1)/2); 
   seq.tmp.rankfiterr = tmpnow.tmp.finalfiterr;
end
[othertmpname,relationtype,matchscore, matchposprob, matchnegprob] = findrelationwithothertmps(tmpnow,alltmp, rrank);
seq.tmp.othertempname = othertmpname; seq.tmp.relationtype = relationtype;
seq.tmp.matchscore = matchscore; seq.tmp.posmatchprob = matchposprob; seq.tmp.negmatchprob = matchnegprob;
%%%%%%%%2.Sequence properties & significance
 
%%%%%%matchnumbers and theoretical significance 
nev = numel(evName); %%%%%for events
%%%%%candidate eventssignificant positive/negative events
totalN = NaN*ones(1, nev); candN = NaN*ones(1, nev); candLapRate = NaN*ones(1, nev); 
evTimeRate = NaN*ones(1, nev); candTimeRate = NaN*ones(1, nev);
posN = NaN*ones(1, nev); posP = NaN*ones(1, nev); posE = NaN*ones(1, nev); posSig = NaN*ones(1, nev); posNSig = NaN*ones(1, nev);
negN = NaN*ones(1, nev); negP = NaN*ones(1, nev); negE = NaN*ones(1, nev); negSig = NaN*ones(1, nev); negNSig = NaN*ones(1, nev);
posLapRate = NaN*ones(1, nev); negLapRate = NaN*ones(1, nev); posTimeRate = NaN*ones(1, nev); negTimeRate = NaN*ones(1, nev);
posReplayRatio = NaN*ones(1, nev); negReplayRatio = NaN*ones(1, nev);
seqdata.events.eventimes = evTimes;  
for (i = 1:nev)
    seqdata.data.evseq{i} = evseq{i};
    [posN(i), posP(i), posE(i), negN(i), negP(i), negE(i), posSig(i), negSig(i),candN(i)] = findtheorysig(evseq{i}, seqparm.significancelevel, seqparm.minactivecell,rrank, seqparm.eventoption);
    totlap = numel(evTimes{i}.start); totalN(i) = totlap;
    %tottime = sum(evTimes{i}.ent-evTimes{i}.start);
    evTimes{i}.start = evTimes{i}.start - seqparm.setbacktime; evTimes{i}.ent = evTimes{i}.ent + seqparm.setbacktime;
    [~, evsessid] = identifysession(evTimes{i}, seq.general.sessionname, seq.general.sessionstartT, seq.general.sessionendT);%%%%work 
    tottime = seq.general.sessionlength(evsessid); 
    posTimeRate(i) = posN(i)/tottime; negTimeRate(i) = negN(i)/tottime;
    posLapRate(i) = posN(i)/totlap; negLapRate(i) = negN(i)/totlap;
    evTimeRate(i) = totalN(i)/tottime; candTimeRate(i) = candN(i)/tottime; candLapRate(i) = candN(i)/totlap; 
    posReplayRatio(i) = posN(i)/candN(i); negReplayRatio(i) = negN(i)/candN(i);
    posNSig(i) = (posN(i)-posE(i))/posSig(i); negNSig(i) = (negN(i)-negE(i))/negSig(i);
end
seq.seq.evTotalN = totalN; seq.seq.evCandN = candN; 
seq.seq.evTotalTimeRate = evTimeRate; seq.seq.evCandTimeRate = candTimeRate; 
seq.seq.evCandLapRate = candLapRate;
seq.seq.evPosMatchN = posN; seq.seq.evPosMatchZ = posNSig; seq.seq.evPosRatio = posReplayRatio;seq.seq.evPosMatchTheoP = posP; seq.seq.evPosTheoExpN = posE; seq.seq.evPosTheoExpSig = posSig;
seq.seq.evPosTimeRate = posTimeRate; seq.seq.evPosLapRate = posLapRate; 
seq.seq.evNegMatchN = negN; seq.seq.evNegMatchZ = negNSig; seq.seq.evNegRatio = negReplayRatio; seq.seq.evNegMatchTheoP = negP; seq.seq.evNegTheoExpN = negE; seq.seq.evNegTheoExpSig = negSig;
seq.seq.evNegTimeRate = negTimeRate; seq.seq.evNegLapRate = negLapRate;
%%%%%shuffle significance
posP = NaN*ones(1, nev); posE = NaN*ones(1, nev); negP = NaN*ones(1, nev); negE = NaN*ones(1, nev); 
posTimeRate = NaN*ones(1, nev); negTimeRate = NaN*ones(1, nev); posLapRate = NaN*ones(1, nev); negLapRate = NaN*ones(1, nev); 
posSig = NaN*ones(1, nev); negSig = NaN*ones(1, nev);
for (i = 1:nev) %for events
    alltimenow = evseq{i}.alltime;
    [~, evsessid] = identifysession(evTimes{i}, seq.general.sessionname, seq.general.sessionstartT, seq.general.sessionendT);%%%%work out which session the event is
    [posP(i), posE(i), negP(i), negE(i), posSig(i), negSig(i)] = findshufflesig(seq.tmp.tmpseq, seq.seq.evPosMatchN(i), seq.seq.evNegMatchN(i), ...
        seq.general.sessionstartT(evsessid), seq.general.sessionendT(evsessid), tmpnow.tmp.tmpfileletter,....
        evseq{i}, rrank, seqparm, 'event', tmprank, seqparm.eventoption); %, alltimenow);
    totlap = numel(evTimes{i}.start); tottime = sum(evTimes{i}.ent-evTimes{i}.start);
    posTimeRate(i) = posE(i)/tottime; negTimeRate(i) = negE(i)/tottime;
    posLapRate(i) = posE(i)/totlap; negLapRate(i) = negE(i)/totlap;
end
seq.seq.evPosMatchShufZ = (seq.seq.evPosMatchN - posE)./posSig; seq.seq.evPosMatchShufP = posP; seq.seq.evPosShufExpN = posE; 
seq.seq.evPosShufExpSig = posSig; seq.seq.evShufposTimeRate = posTimeRate; seq.seq.evShufposLapRate = posLapRate;
seq.seq.evNegMatchShufZ = (seq.seq.evNegMatchN - negE)./negSig; seq.seq.evNegMatchShufP = negP; seq.seq.evNegShufExpN = negE; 
seq.seq.evNegShufExpSig = negSig; seq.seq.evShufnegTImeRate = negTimeRate; seq.seq.evShufnegLapRate = negLapRate;

function [seq, seqstart, seqend, seqmarker, matchscore, posprob, negprob, seqtime] = breakfromsessionallsequences(sessseq, event)
iii = find( ((sessseq.seqstart{1,1}>=event.start)&(sessseq.seqstart{1,1}<=event.ent)) &  ((sessseq.seqend{1,1}>=event.start)&(sessseq.seqend{1,1}<=event.ent)));
seq = sessseq.seq{1,1}(iii); seqstart = sessseq.seqstart{1,1}(iii); seqend = sessseq.seqend{1,1}(iii); seqmarker = sessseq.seqmarker{1,1}(iii);
matchscore = sessseq.matchscore{1,1}(iii); posprob = sessseq.posmatchprob{1,1}(iii); negprob = sessseq.negmatchprob{1,1}(iii); seqtime = sessseq.seqtime{1,1}(iii);

function [seq, seqstart, seqend, seqmarker, score, posprob, negprob, seqtime] = ...
               breakallsequences_nonselect(allletter, alltime, event, tmprank, rrank, tmpseq, seqparm)
maxTleng = seqparm.maxtmpleng; maxchecktime = seqparm.maxtime; %20 s
%%%%%break long
iii = find( (alltime>=event.start) & (alltime<=event.ent) ); 
nletter = numel(iii); maxcheckleng = min([nletter maxTleng*numel(tmpseq)]);
letternow = allletter(iii); timenow = alltime(iii); 
seq = []; seqstart = []; seqend = []; Ctime = []; nseq = 0;
for (leng = 4:maxcheckleng) %%%sequence length under consideration
    for (j = 1:nletter-leng+1)
        st = timenow(j); et = timenow(j+leng-1);
        if (et-st<maxchecktime)
            nseq = nseq + 1; seq{nseq} = letternow(j:j+leng-1); Ctime{nseq} = timenow(j:j+leng-1);
            seqstart(nseq) = st; seqend(nseq) = et;
        end
    end
end
nseq = numel(seq); score = NaN*ones(nseq,1); posprob = NaN*ones(nseq,1); negprob = NaN*ones(nseq,1); seqmarker = cell(nseq, 1);
if ~isempty(tmprank)
   for (i = 1:nseq)
       %[score(i), posprob(i), negprob(i)] = findseqscoreprob(seq{i}, tmpseq, rrank);
       [score(i), posprob(i), negprob(i)] = findseqscoreprob_multiple(seq{i}, Ctime{i}, tmpseq, tmprank); %, rrank, seqparm); %%%%
   end
end
for (i = 1:nseq) seq{i} = char(seq{i}); seqmarker{i} = event.marker; end
seq = seq'; seqstart = seqstart'; seqend = seqend'; seqtime = Ctime';
function [score, posprob, negprob] = findseqscoreprob_multiple(Cstr, Ctime, tmpseq, tmprank) %, rrank, seqparm) %%%%
%nshuffle = seqparm.Sigshuffle; %default shuffle times to prob significance
score = NaN; posprob = NaN; negprob = NaN; Cut = unique(Cstr); nlet = numel(Cut); neq = numel(Cstr);
if (nlet >= 4) %%%if unique letters more than 4
   %%%%work out str ranks
       Cstrrank = zeros(1, neq); for (i = 1:neq) Cstrrank(i) = tmprank(tmpseq == Cstr(i)); end
       [RR,PP] = corrcoef(Cstrrank, Ctime); 
       score = RR(1,2); 
       if (score >= 0)
           posprob = PP(1,2); 
       else
           negprob = PP(1,2);
       end
end

function [allletter, alltime] = findfreesequences(spikedata, matchletter, PeakTh, timepoint) %%%first, map out entire
nspike = numel(spikedata); peakL = cell(1,nspike); peakT = cell(1,nspike);
for (i = 1:nspike)
    [~, maxind] = FindLocalMinima(spikedata{i}); maxV = spikedata{i}(maxind);
    iii = find( maxV > PeakTh(i) );  peakT{i} = timepoint(maxind(iii)); peakL{i} = matchletter(i)*ones(1,numel(iii));
end
allL = cell2mat(peakL); allT = cell2mat(peakT);
[alltime, iii] = sort(allT); allletter = allL(iii);

function [allletter, alltime] = findfreesequences_space(spikedata, matchletter, PeakTh, timepoint, lapnum)
nspike = numel(spikedata); peakL = cell(1,nspike); peakT = cell(1,nspike);
for (i = 1:nspike)
    [~, maxind] = FindLocalMinima(spikedata{i}{lapnum}); maxV = spikedata{i}(maxind);
    iii = find( maxV > PeakTh(i) );  peakT{i} = timepoint(maxind(iii)); peakL{i} = matchletter(i)*ones(1,numel(iii));
end
allL = cell2mat(peakL); allT = cell2mat(peakT);
[alltime, iii] = sort(allT); allletter = allL(iii);

function tmpnow = findtemplates(tmp, animaldatenow)
int = zeros(1, numel(tmp));
for (i = 1:numel(tmp))
    if strcmp(tmp(i).parm.animaldate, animaldatenow) 
        int(i) = 1; 
    else
        animdatethen = rearrangeanimdate(tmp(i).parm.animaldate);
        if strcmp(animdatethen, animaldatenow) 
           int(i) = 1;
        end
    end
end
tmpnow = tmp(int==1);
function animdatethen = rearrangeanimdate(animaldate) %%%translate DE animaldate to DM animname_date
[str, tok] = strtok(animaldate, ':\'); animaldate = tok(3:numel(tok));
[str, tok] = strtok(animaldate, '\'); animnow = str; datenow = tok(2:numel(tok)); 
animdatethen = strcat(animnow, '_', datenow);

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

function data = getfiringrates_binspikes(spiketime, smoothopt, sigma, timepoint, binsize, sT, eT)
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
   if strcmp(smoothopt, 'Yes')
      ybin = conv(ybin, wind);
      data = (1/binsize)*ybin( (nw-1)/2+1 : numel(ybin)-(nw-1)/2 ); %%%final results in Hz
   else
      data = ybin;
   end
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

function [seq, seqstart, seqend, seqmarker, score, posprob, negprob] = ...
    findallsequence(spikedata, letter, timingoption, event, PeakTh, timepoint, rrank, tmpseq)
%%%%find cell spiking sequence
%%%%Input: spikedata{i} = i-th spike timestamps in timeunit second
%%%%       letter(i) = single character letter assigned to i-th cell
%%%%       option = timing option (0 = first spike; 1 = mass center; 2 = peak)
%%%%       event file (ep.start, ep.ent, ep.marker)
%%%%       all times in second
%%%%Output:
%%%%       seq{i} = i-th spiking sequence found (could be singleton, doublet, tripplet, or higher order)
%%%%       seqstart(i) = i-th sequence starting time in second
%%%%       seqend(i) = i-th sequence ending time in second

%%%%%%%%%%randomize to eliminate the bias toward the canonical sequence (01234...)
%%%%%%%%%%       ----when peak times are the same, the sort function does not change the indices.
nspike = numel(spikedata); 
IXX = randperm(nspike); spikedata = spikedata(IXX); letter = letter(IXX); PeakTh = PeakTh(IXX); 
nev = numel(event.start); seq = cell(nev, 1); evstart = event.start; evend = event.ent;
score = NaN*ones(nev,1); posprob = NaN*ones(nev,1); negprob = NaN*ones(nev,1); 
seqstart = NaN*ones(nev, 1); seqend = NaN*ones(nev,1); seqmarker = cell(nev,1);
if iscell(event.marker)
   for (i = 1:nev) seqmarker{i} = event.marker{i}; end
else
   for (i = 1:nev) seqmarker{i} = event.marker; end
end
for (i = 1:nev)
    timenow = zeros(1, nspike);
    if strcmp(timingoption, 'RatePeak_time') %%%if peak time
        tid = find( (timepoint>=evstart(i)) & (timepoint<=evend(i)) ); TT = timepoint(tid); 
        for (j = 1:nspike)
            RR = spikedata{j}(tid); [PP,ii] = max(RR);
            if (PP>0)&&(PP > PeakTh(j)) timenow(j) = TT(ii); end
        end
    elseif strcmp(timingoption, 'Rate1stPeak_time') %%%if peak time
        tid = find( (timepoint>=evstart(i)) & (timepoint<=evend(i)) ); TT = timepoint(tid); 
        for (j = 1:nspike)
            dat = spikedata{j}(tid); 
            [~, maxind] = FindLocalMinima(dat); maxV = dat(maxind);
            iii = find( (maxV > PeakTh(j)) & (maxV > 0) );  peakT = TT(maxind(iii)); 
            if ~isempty(peakT) timenow(j) = min(peakT); end
        end
    elseif strcmp(timingoption, 'RatePeak_space') %%%if peak space
        for (j = 1:nspike) spikedatathen{j} = spikedata{j}{i}; end
        tid = find( (timepoint>=evstart(i)) & (timepoint<=evend(i)) ); TT = timepoint(tid);
        for (j = 1:nspike)
            RR = spikedatathen{j}(tid); [PP,ii] = max(RR);
            if (PP>0)&&(PP > PeakTh(j)) timenow(j) = TT(ii); end
        end
    elseif strcmp(timingoption, 'Rate1stPeak_space') %%%if peak space
        for (j = 1:nspike) spikedatathen{j} = spikedata{j}{i}; end
        tid = find( (timepoint>=evstart(i)) & (timepoint<=evend(i)) ); TT = timepoint(tid);
        for (j = 1:nspike)
            dat = spikedatathen{j}(tid);
            [~, maxind] = FindLocalMinima(dat); maxV = dat(maxind);
            iii = find( (maxV > PeakTh(j)) & (maxV > 0) );  peakT = TT(maxind(iii)); 
            if ~isempty(peakT) timenow(j) = min(peakT); end
        end
    else
        for (j = 1:nspike)
            sid = find( (spikedata{j}>=evstart(i)) & (spikedata{j}<=evend(i)) );
            if (~isempty(sid))
                if strcmp(timingoption, 'FirstSpike')
                    timenow(j) = min(spikedata{j}(sid)); 
                elseif strcmp(timingoption, 'SpikeTrainCenter')
                    timenow(j) = mean(spikedata{j}(sid)); 
                end
            end
        end
    end
    %%%%%%%output sequence
    iii = find(timenow ~= 0); 
    if (~isempty(iii))
       timethen = timenow(iii); letterthen = letter(iii); 
       [~, IX] = sort(timethen); letterthen = letterthen(IX); seq{i} = char(letterthen);
       %if ~isempty(timethen) 
           seqstart(i) = min(timethen); seqend(i) = max(timethen);
       %end
       if ~isempty(rrank)
          [score(i), posprob(i), negprob(i)] = findseqscoreprob(letterthen, tmpseq, rrank); %%%%this works for straight template: 0123456A, etc.
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
if (numel(unique(Cstr)) > 1)
   for (i = 1:nq-1)
       for (j = i+1:nq)
           matchnow = tmpind(j)-tmpind(i);
           if (matchnow>0)
               PP = PP +1;
           elseif (matchnow<0)
               NN = NN + 1;
           end
       end
   end
   score = (PP-NN)/(PP+NN);
   epsilon = 0.001; posprob = 0; negprob = 0;
   posindex = find(rrank.ratio{nq} >= score-epsilon);
   negindex = find(rrank.ratio{nq} <= score+epsilon);
   if (~isempty(posindex))
       posprob = sum(rrank.prob{nq}(posindex));
   end
   if (~isempty(negindex))
       negprob = sum(rrank.prob{nq}(negindex));
   end
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
       spatX = bhdata.event.Xbin{posid}{evid}; framerate = behav.parm.framerate(posid);
       nlap = numel(evTime.start); spatrate = cell(1, nlap);
       for (jnow = 1:nlap)
            occutimenow = bhdata.event.LapOccuptime{posid}{evid}{jnow}; 
            lappostime{1} = bhdata.event.LapAllPostimestamp{posid}{evid}{jnow}; lapx{1} = bhdata.event.LapAllX{posid}{evid}{jnow};
            evok.start = evTime.start(jnow); evok.ent = evTime.ent(jnow);
            [spatrate, ~] = getlinearmap(spiketime, evok, lappostime, lapx, spatX, occutimenow, sigma, Nsigma, 1/framerate);
       end 
end
function [D1rate, nact] = getlinearmap(spiketime, evTime, lappostime, lapx, xbin, occutime, smoothopt, sigma, Nsigma, frametime)
nx = numel(xbin); D1rate = zeros(nx,1); nlap = numel(lappostime); count = zeros(nx,1); spikex = []; nact = NaN;
frametime = 1.03*frametime; %%%slightly expand frametime to assure detection, due to slightly variable frame-to-frame time
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
             %if (laptimenow(k) <= spikenow(j)) && (laptimenow(k+1) > spikenow(j)) 
             if abs(laptimenow(k) - spikenow(j)) <= frametime 
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
if (strcmp(smoothopt, 'Yes'))
%%now do a 1d smoothing
if (strcmpi(smParm.d1sm, 'yes'))
    %D1rate = OneDSmooth_new(D1rate, occutime, sigma, smParm.d1Nsig);
    XX = [-Nsigma*sigma:1:Nsigma*sigma]; 
    wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
    D1rate = conv(D1rate, wind); D1rate = D1rate((nw-1)/2+1:numel(D1rate)-(nw-1)/2);
end
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
    evthen = alltmp(iii(i)).evtfile; seqthen = alltmp(iii(i)).tmp.tmpfinalseq{1}; 
    if isfield(alltmp(iii(i)).parm, 'sessev')
        sessevthen = alltmp(iii(i)).parm.sessev;
    else
        [str, tok] = strtok(evthen, '_'); sessevthen = str;
    end
    if isfield(tmpnow.parm, 'sessev')
        sessevnow = tmpnow.parm.sessev;
    else
        [str, tok] = strtok(tmpnow.evtfile, '_'); sessevnow = str;
    end
    relationtype{i} = findrelationtype(evthen, sessevthen, tmpnow.evtfile, sessevnow);
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

function [npos, nneg] = findmatchnum(evposprob, evnegprob, evseqstart, evseqend, slevel, eventoption)
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
    if ~strcmp(eventoption, 'WithinEvents')
       [newstart, ~, ~, ~] = resolveoverlap(posprob(iii), seqstart(iii), seqend(iii)); 
    else
        newstart = seqstart(iii);
    end
    npos = npos + numel(newstart);
    iii = find(negprob<slevel);%%%%negative matches
    if ~strcmp(eventoption, 'WithinEvents')
        [newstart, ~, ~, ~] = resolveoverlap(negprob(iii), seqstart(iii), seqend(iii)); 
    else
        newstart = seqstart(iii);
    end
    nneg = nneg + numel(newstart);
end
function [npos, posP, posmm, nneg, negP, negmm, posSig, negSig, candN] = findtheorysig(evseq, slevel, minactivecell, rrank, eventoption)
%%%%%%here seq{i,j} corresponds to ep{i,j}
%%%%%found sig matches need to resolve timing overlaps within laps (cross all window size if sliding windows)
%%%%%theory significance produced only if eventoption = WithinEvents
posP = NaN; posmm = NaN; negP = NaN; negmm = NaN; npos = 0; nneg = 0; posSig = NaN; negSig = NaN; 
[mm,nn] = size(evseq.seq); seq = []; candN = 0; %%%seq variable only useful for "WithinEvents'
for (i = 1:mm)
     posprob = []; negprob = []; mscore = []; seqstart = []; seqend = [];
    for (j = 1:nn)
        seq = [seq; evseq.seq{i,j}]; posprob =[posprob; evseq.posmatchprob{i,j}]; negprob = [negprob; evseq.negmatchprob{i,j}];
        mscore = [mscore; evseq.matchscore{i,j}]; seqstart =[seqstart; evseq.seqstart{i,j}]; seqend = [seqend; evseq.seqend{i,j}];
    end
    iii = (find(posprob<slevel));%%%%positive matches
    if ~strcmp(eventoption, 'WithinEvents')
       [newstart, ~, ~, ~] = resolveoverlap(posprob(iii), seqstart(iii), seqend(iii)); 
    else
        newstart = seqstart(iii);
    end
    npos = npos + numel(newstart);
    iii = find(negprob<slevel);%%%%negative matches
    if ~strcmp(eventoption, 'WithinEvents')
        [newstart, ~, ~, ~] = resolveoverlap(negprob(iii), seqstart(iii), seqend(iii)); 
    else
        newstart = seqstart(iii);
    end
    nneg = nneg + numel(newstart);
end
if strcmp(eventoption, 'WithinEvents') %%%not available to 'FreeSequences', dangerous to 'Sliding windows'
   %npos = numel(find(posprob<slevel)); nneg = numel(find(negprob<slevel)); %%%%positive/negative matches
   %%%%%expected distribution of sig seq
   [posexpprob, negexpprob] = findexpprob(rrank, slevel);
   seqleng = zeros(1, numel(seq));
   for (i = 1:numel(seq)) seqleng(i) = numel(seq{i}); end
   candN = numel(find(seqleng>=minactivecell));
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

function [posP, posE, negP, negE, posSig, negSig] = findshufflesig(tmpseq, npos, nneg, sessST, sessET, fileletter, evseq, rrank, seqparm, sessev, tmprank, eventoption) %, alltime)
posP = NaN; posE = NaN; negP = NaN; negE = NaN; posSig = NaN; negSig = NaN;
%%%%%%here seq{i,j} corresponds to ep{i,j}
%%%%first generate shuffles
ok = 1; nshuffle = seqparm.Nshuffle; slevel = seqparm.significancelevel;
if (strcmp(sessev, 'session'))
    eventoption = seqparm.sessionoption; timingoption = seqparm.sessiontimingoption;
else
    eventoption = seqparm.eventoption; timingoption = seqparm.eventtimingoption;
end
if strcmp(seqparm.shufflemode, 'GlobalID')
   [~, ~, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_gloabalID(tmpseq, fileletter, evseq, seqparm, rrank, tmprank, eventoption);%, alltime);
elseif strcmp(seqparm.shufflemode, 'LocalID')
   [~, ~, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_localID(tmpseq, fileletter, evseq, seqparm, rrank, tmprank, eventoption);%, alltime);
elseif strcmp(seqparm.shufflemode, 'CircleSlide')
   [~, ~, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_circleslide(tmpseq, fileletter, evseq, seqparm, ...
       sessST, sessET, rrank, tmprank, eventoption, timingoption);
else
    ok = 0; disp('----------> can not generate shuffled data or shuffle mode not available; shuffle significance not computed');
end
if ok
   shufnpos = zeros(1,nshuffle); shufnneg = zeros(1,nshuffle);
   for (k = 1:nshuffle)
     [shufnpos(k), shufnneg(k)] = findmatchnum(posmatchprob{k}, negmatchprob{k}, seqstart{k}, seqend{k}, slevel, eventoption); 
   end
   posP = numel(find(shufnpos>=npos))/nshuffle; negP = numel(find(shufnneg>=nneg))/nshuffle;
   posE = mean(shufnpos); negE = mean(shufnneg); posSig = std(shufnpos); negSig = std(shufnneg);
end
function [seq, matchscore, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_gloabalID(tmpseq, fileletter, evseq, seqparm, rrank, tmprank, eventoption) %alltime
nshuffle = seqparm.Nshuffle;
seq = cell(1,nshuffle); matchscore = cell(1, nshuffle); posmatchprob = cell(1,nshuffle); negmatchprob = cell(1, nshuffle);
seqstart = cell(1, nshuffle); seqend = cell(1, nshuffle); 
%%%%%%%%%%%%%%%%%still evseq.seq{i,j} corresponds to ep{i,j}
%%%%%%%%%%for global ID shuffling: just need to randomly swap cell IDs and then replace original sequences 
sequence = evseq.seq; ncell = numel(tmpseq); [mm,nn] = size(sequence);
tmpind = zeros(1, ncell); for (i = 1:ncell) tmpind(i) = find(fileletter==tmpseq(i)); end
for (i = 1:nshuffle)
   %tt = cputime;
    %%%%%these variable do not change
    seqstart{i} = evseq.seqstart; seqend{i} = evseq.seqend; 
    seq{i} = cell(mm,nn); matchscore{i} = cell(mm,nn); posmatchprob{i} = cell(mm,nn); negmatchprob{i} = cell(mm,nn); 
    IX = randperm(ncell); mletternow = fileletter(tmpind(IX));
    for (ii = 1:mm)
        for (jj = 1:nn)
            for (j = 1:numel(sequence{ii,jj}))
                 seqnow = sequence{ii,jj}{j}; nletter = numel(seqnow); ind = zeros(1, nletter);
                 for (k = 1:nletter)
                      ind(k) = find(tmpseq == seqnow(k));
                 end
                 if isempty(find(ind == 0))
                    seq{i}{ii,jj}{j} = mletternow(ind);
                    if strcmp(eventoption, 'FreeSequences')
                       Ctime = evseq.seqtime{ii,jj}{j}; 
                       [matchscore{i}{ii,jj}(j), posmatchprob{i}{ii,jj}(j), negmatchprob{i}{ii,jj}(j)] = findseqscoreprob_multiple(...
                           seq{i}{ii,jj}{j}, Ctime, tmpseq, tmprank); %, rrank, seqparm);
                    else
                       [matchscore{i}{ii,jj}(j), posmatchprob{i}{ii,jj}(j), negmatchprob{i}{ii,jj}(j)] =...
                         findseqscoreprob(seq{i}{ii,jj}{j}, tmpseq, rrank);
                    end
                 end
            end
        end
    end
    %disp(['----------> single shuffle time: ', num2str(cputime-tt)]);
end
function [seq, matchscore, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_localID(tmpseq, fileletter, evseq, seqparm, rrank, tmprank, eventoption) %, alltime)
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
                    if strcmp(eventoption, 'FreeSequences')
                       Ctime = evseq.seqtime{ii,jj}{j}; %alltime( (alltime>=seqstart{i}{ii,jj}(j)) & (alltime<=seqend{i}{ii,jj}(j)) );
                       [matchscore{i}{ii,jj}(j), posmatchprob{i}{ii,jj}(j), negmatchprob{i}{ii,jj}(j)] = findseqscoreprob_multiple(...
                           seq{i}{ii,jj}{j}, Ctime, tmpseq, tmprank); %, rrank, seqparm);
                    else
                       [matchscore{i}{ii,jj}(j), posmatchprob{i}{ii,jj}(j), negmatchprob{i}{ii,jj}(j)] =...
                         findseqscoreprob(seq{i}{ii,jj}{j}, tmpseq, rrank);
                    end
                 end
            end
        end
    end
end
function [seq, matchscore, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_circleslide(tmpseq, fileletter, evseq, seqparm, sessST, sessET, rrank, tmprank, eventoption, timingoption)
nshuffle = seqparm.Nshuffle; 
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
               breakallsequences_nonselect(allletter, alltime, ep{ii,jj}, tmprank, rrank, tmpseq, seqparm);
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

function [evName, evType, evT] = filterevents(pinfo, data, evkeyword, evkeytype, neuronid)
evName = pinfo.general.eventname{neuronid(1)}; evsel = ones(1,numel(evName));
evType = pinfo.parm.eventtype{neuronid(1)};
evT = data.events.eventtimes{neuronid(1)}; 
if ~isempty(evkeyword)
   for (i = 1:numel(evName))
       if isempty(strfind(lower(evName{i}), lower(evkeyword))) evsel(i) = 0; end 
   end
end
if ~isempty(evkeytype)
   for (i = 1:numel(evName))
       if isempty(strfind(lower(evType{i}), lower(evkeytype))) evsel(i) = 0; end 
   end
end
evpos = find(evsel == 1);
evName = evName(evpos); evType = evType(evpos); evT = evT(evpos);




