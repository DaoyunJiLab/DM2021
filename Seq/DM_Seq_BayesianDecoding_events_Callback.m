function DM_Seq_BayesianDecoding_events_Callback
%%Do Bayesian decoding analysis on a database (.spikedb), require a matching place-field-generated .tmp file
%%This one ONLY DECODEs A 1D SPACE!!
%%%%%%%%%%%%%%%ONLY DO IT WITHIN EVENTS SPECIFIED BY A KEYWORD%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%work on selected cells (groups), group cells into days, and work on each day
%%%%%%%% find if templates exist and then work on each template (if no template found for a date, warn and skip)
%%%%%%%%%%%%%%%%%%%for each event: output decoded prob and replay quantifications
%%%%event parsing options: within events
hf = gcbf; hgroup = getappdata(hf, 'hgroup'); ok = 1; disp('-----> 1D Bayesian decoding: using templates as items');
fname = get(hf, 'Name'); ftt = strfind(fname, '__'); currentfilename = fname(ftt+2:numel(fname));
[~, ~, ee] = fileparts(currentfilename);
if ~strcmp(ee, '.spikedb')
    disp('--------> not a spikedb; aborted'); ok = 0;
else
    pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); behav = []; bhdata = [];
    plotparm = getappdata(hf, 'plotparm');
    hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
    for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end    
    if (plotparm.linkbehav == 0)
       disp(['--------> no behav data linked']); ok = 0;
    else
       behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
    end
end
if ok
   seqparm.seqtype = 'Bayesian'; seqparm.eventoption = 'WithinEvents'; 
   evsss = {'WithinEvents'; 'FreeSequences'};
   [sss, ok] = listdlg('ListString', evsss, 'PromptString', 'Event parsing option');
   if ok seqparm.eventoption = evsss{sss};  end
end
if ok
    input = inputdlg({'Event keyword'; 'Event type'; 'Event keyNOword'; 'Event NOtype'}, 'Event selection', 4, {'track'; 'run'; 'first'; ''}); 
    if (~isempty(input))
        seqparm.evkeyword = input{1}; seqparm.evkeytype = input{2}; seqparm.evkeynoword = input{3}; seqparm.evkeynotype = input{4};
    else
        ok = 0;
    end
end
if ok
    input = inputdlg({'sessioin/events?'; 'Event keyword'; 'Event type'},...
        'Total time options', 3, {'session'; 'SWS'; 'sws'}); 
    if (~isempty(input))
       seqparm.tottimeoption = input{1}; 
       seqparm.tottimekeyword = input{2}; seqparm.tottimekeytype = input{3};
    else
       ok = 0;
    end
end
if ok
    input = inputdlg({'Peak prob threshold(x)'; 'Min active cell #'; 'Window size (s)'; 'Window shift time (s)'; 'Setback time (s)'},...
        'Decoding parameters', 5,  {'1'; '4'; '0.3'; '0.15'; '0.02'}); % {'1'; '4'; '0.3'; '0.03'; '0'}); 
    if (~isempty(input))
       seqparm.probthres = str2num(input{1}); seqparm.minactivecell = str2num(input{2});
       seqparm.windowsize = str2num(input{3}); seqparm.shifttime = str2num(input{4}); seqparm.setbacktime = str2num(input{5}); 
    else
       ok = 0;
    end
end
if ok
    if strcmp(seqparm.eventoption, 'WithinEvents')
       input = inputdlg({'Significance level'; 'Number of local shuffles'; 'Number of global shuffles'; 'Gloabal shuffle mode (GlobalID/CircularSlide)'; 'Shuffle data save (raw/quant) rates (1/#)'},...
        'Significance parameters', 5, {'0.05'; '20'; '5'; 'GlobalID'; '5 5'}); 
       if (~isempty(input))
         seqparm.significancelevel = str2num(input{1}); seqparm.Nwithineventshuffle = str2num(input{2}); 
         seqparm.Nshuffle = str2num(input{3}); seqparm.shufflemode = input{4}; T = str2num(input{5});
         if (numel(T)>=2) seqparm.saveshufrawdatarate = T(1); seqparm.saveshufquandatarate = T(2);
         elseif (numel(T)==1) seqparm.saveshufrawdatarate = T(1); seqparm.saveshufquandatarate = T(1);
         else
           ok=0;
         end
       else
         ok = 0;
       end
    elseif strcmp(seqparm.eventoption, 'FreeSequences')
       input = inputdlg({'Significance level'; 'Number of global shuffles'; 'Gloabal shuffle mode (GlobalID/CircularSlide)'; 'Shuffle data save (raw/quant) rates (1/#)'},...
        'Significance parameters', 5, {'0.05'; '5'; 'GlobalID'; '5 5'}); 
       if (~isempty(input))
         seqparm.significancelevel = str2num(input{1}); seqparm.Nwithineventshuffle = 0;  
         seqparm.Nshuffle = str2num(input{2}); seqparm.shufflemode = input{3}; T = str2num(input{4});
         if (numel(T)>=2) seqparm.saveshufrawdatarate = T(1); seqparm.saveshufquandatarate = T(2);
         elseif (numel(T)==1) seqparm.saveshufrawdatarate = T(1); seqparm.saveshufquandatarate = T(1);
         else
           ok = 0;
         end
         if ok
              input = inputdlg({'Free sequence min time points'; 'max time points'}, 'Leng parameters', 2, {'30'; '10000'}); 
              if (~isempty(input))
                 seqparm.freeMinTimepoints = str2num(input{1}); seqparm.freeMaxTimepoints = str2num(input{2}); 
              else
                 ok=0;
              end
         end
       else
         ok = 0;
       end
    end
end
if ok
    cpathname = fullfile(cd, '*.tmp');
    [fname, pname] = uigetfile(cpathname, 'Select a template file to open:');
    if fname ~= 0
       tmpfile = fullfile(pname, fname);
       disp(['--------> template file: ', tmpfile]); seqparm.tmpfile = tmpfile;
       S = load(tmpfile, '-mat'); tmp = S.pinfo; S = [];
       if strcmp(tmp(1).tmp.tmptype, 'pfratetmp')
          disp(['--------> number of templates found: ', num2str(numel(tmp))]);
       else
          ok = 0; disp('--------> templates not built on 1D spaces; aborted.');
       end
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
    %%%%output results to database structure
    seq.work = struct([]);
    [seq, seqdata] = computeseqdb(seq, tmp, pinfo, data, behav, bhdata, cellind, seqparm);  
    if (~isempty(seqdata))
       %%%%%group data
       seqdata.grouplist.groupname{1} = 'List0'; seqdata.grouplist.groupindex{1} = 1:numel(seq.general.tmpID);
       seqdata.grouplist.grouptype{1} = 'Manual'; seqdata.grouplist.groupcrit{1} = []; seqdata.grouplist.groupparents{1} = []; 
       seqdata.parentfile = currentfilename; seqdata.seqmethod = 'Bayesian';
       pinfo = seq; 
       seqfilestripped = strrep(seqfile, '.seqdb', '_stripped.seqdb');
       data = []; data.grouplist = seqdata.grouplist; 
       save(seqfilestripped, 'pinfo', 'data', '-mat');
       %%%save a full version
       data = seqdata; save(seqfile, 'pinfo', 'data', '-mat', '-v7.3'); %%% This allows to save files bigger than 2GB on 64-bit systems
    end
end
disp('************************');

function [seq, seqdata] = computeseqdb(seq, tmp, pinfo, data, behav, bhdata, cellind, seqparm)
seqdata = []; nt = 0;
%%%%search for all the possible finaldir
allfinaldir = pinfo.general.finaldir(cellind); 
unifinaldir = unique(allfinaldir);
disp(['--------> number of final dirs: ', num2str(numel(unifinaldir))]);
for (i = 1:numel(unifinaldir))
    fdirnow = unifinaldir{i};
    disp(['-----------> final dir now: ', fdirnow]);
    ind = find(strcmp(allfinaldir, unifinaldir{i}));
    neuronid = cellind(ind); %%%all cells in this final directory 
    cellparmfile = pinfo.general.parmfile(neuronid); 
    if (numel(cellparmfile)<4)
        disp(['-----------> warning: less than 4 cells: ', fdirnow]);
    end
    animaldatenow = strcat(pinfo.general.animalname{neuronid(1)}, '_', pinfo.general.datedir{neuronid(1)});
    tmpnow = findtemplates(tmp, animaldatenow);
  if ~isempty(tmpnow)  
    for (j = 1:numel(tmpnow)) %%%for each template
        disp(['--------------> template now: ', tmpnow(j).tmp.tmpfinalseq{1}, '; ', tmpnow(j).evtfile]);
        nt = nt + 1; %%%entry of the seqdb items
        t0 = cputime; %disp('-------------> start cputime');
        [evseq, shufevseq, tmpnow(j), evName, evType, evT, tottime] = sequencedataandshuffle(pinfo, data, behav, bhdata, neuronid, tmpnow(j), seqparm); %evseq{nevent file}.field{nevents, 1}(1); 
        %t1 = cputime; disp(['-------------> decoding cpu time: ', num2str(t1 - t0)]);
        [seqnow, seqdatanow] = computeassignseqproperties(evseq, shufevseq, seqparm, pinfo, neuronid, tmpnow(j), evName, evType, evT, tottime); %shufevseq{nshuf, nevent file}.field.(morefield){nevents, 1}(1)
        t2 = cputime; disp(['-------------------> decoding/shuffling/quantification cpu time: ', num2str(t2 - t0)]);
        %t2 = cputime; disp(['-------------> quantification cpu time: ', num2str(t2 - t1)]);
        [seq, seqdata] = assignseqoutput(seq, seqdata, seqnow, seqdatanow, nt);
    end
  end
end

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

function [seq, seqdata] = computeassignseqproperties(evseq, shufevseq, seqparm, pinfo, allneuronid, tmpnow, evName, evType, evTimes, tottime) 
% shufevseq{nshuf, nevent file}.field.(morefield){nevents, 1}(1)
%%%%Sequence properties to compute: 1-template propert; 2-decoded position properties; 3-global significance
%%%%%general variables
neuronid = allneuronid(1);
seq.general.tmpID = strcat(tmpnow.parm.animaldate, '_', tmpnow.evtfile); 
seq.general.recarea = pinfo.general.recarea{neuronid}; seq.general.finaldir = pinfo.general.finaldir{neuronid};
seq.general.datedir = pinfo.general.datedir{neuronid}; seq.general.animalname = pinfo.general.animalname{neuronid};
seq.general.sessionname = pinfo.general.sessionname{neuronid}; seq.general.sessionstartT = pinfo.general.sessionstartT{neuronid};
seq.general.sessionendT = pinfo.general.sessionendT{neuronid}; seq.general.sessionlength = pinfo.general.sessionlength{neuronid};
seq.general.eventname = evName; seq.general.genotype = pinfo.general.genotype{neuronid};
seq.general.sex = pinfo.general.sex{neuronid}; seq.general.age = pinfo.general.age{neuronid};
%%%%parameters
seq.parm = seqparm; seq.parm = rmfield(seq.parm, 'tmpfile'); 
seq.parm.sessType = pinfo.parm.sessType{neuronid}; seq.parm.eventType = evType;
%%%%%%%%%1.template properties - Not done: quality and similarities with other templates in the same finaldir
seq.tmp.tmpfile = seqparm.tmpfile; 
seq.tmp.tmpseq = tmpnow.tmp.tmpfinalseq{1}; 
seq.tmp.tmprank = []; %tmpnow.tmp.tmpfinalrank{1};
seq.tmp.cellnum = numel(tmpnow.tmp.tmpfinalseq{1}); seq.tmp.length = []; %max(seq.tmp.tmprank)-min(seq.tmp.tmprank);
seq.tmp.sessevname = tmpnow.evtfile;
seq.tmp.crrmode = tmpnow.parm.crrmode; 
if isfield(tmpnow.parm, 'sessev') 
    sessev = tmpnow.parm.sessev;
else
    sessev = 'unknown';
end
seq.tmp.sessev = sessev;
seq.tmp.filename = tmpnow.tmp.tmpfilename;  nlet = numel(tmpnow.tmp.tmpfileletter);
for (i = 1:nlet) seq.tmp.fileletter{i} = tmpnow.tmp.tmpfileletter(i); end
% [othertmpname,relationtype,matchscore, matchposprob, matchnegprob] = findrelationwithothertmps(tmpnow,alltmp, rrank);
% seq.tmp.othertempname = othertmpname; seq.tmp.relationtype = relationtype;
% seq.tmp.matchscore = matchscore; seq.tmp.posmatchprob = matchposprob; seq.tmp.negmatchprob = matchnegprob;
%%%%%%%%2.Decoded position properties: use shufposP/shufnegP to define significant events
nev = numel(evName); 
[seqnow, candevind, matchevind, posevind, negevind] = findquantifications(evseq, tottime, seqparm, nev); %%%This is the final output variable: P2D.(field)(1, nev)
%%%%% 3. shuffle significance
signow = findshufflesignicance(seqnow, shufevseq, tottime, seqparm, nev);
%%%%%%Output structure here
seq.seq = seqnow; seq.shufsig = signow; 
%%%%%%%% raw data output
seqdata.events.eventimes = evTimes; seqdata.events.candevind = candevind; 
seqdata.events.matchevind = matchevind; seqdata.events.posevind = posevind; seqdata.events.negevind = negevind;
seqdata.data.evseq = evseq; seqdata.data.shufevseq = shufevseq;

function sig1D = findshufflesignicance(P1D, shufevseq, tottime, seqparm, nev)
%%%%% evseq{nev}.(field){mm,nn}(1); shufevseq{nshuf, nev}.(field){mm,nn}(1)
sig1D = []; nshuf = seqparm.Nshuffle; S1D = cell(1, nshuf); evseqnow = cell(1, nev);
for (i = 1:nshuf) 
    for (k=1:nev) evseqnow{k} = shufevseq{i,k}; end
    [S1D{i}, ~, ~, ~, ~] = findquantifications(evseqnow, tottime, seqparm, nev);
end
%%%%variable Z scores and P values relative to shuffles
S1D = reformatcleanup(S1D); %%%%re-format S1D{nshuf}.(field)(1:nev) to S2D.(field){nev}(nshuf)
fname = fieldnames(P1D);
for (j = 1:numel(fname))% For each variable/field, compare P2D.(field)(1:nev) with S2D.(field)(nshuf:nev)
    f1name = strcat(fname{j}, 'Z'); f2name = strcat(fname{j}, 'P');
    sig1D.(f1name) = NaN*ones(1, nev); sig1D.(f2name) = NaN*ones(1, nev);
    if ~isempty(S1D)
      for (i=1:nev)
         vnow = S1D.(fname{j}){i}; 
         vnow = vnow(~isnan(vnow)); nn = numel(vnow); 
         pnow = P1D.(fname{j})(i);
         if nn >=3  
            shufmean = mean(vnow); shufstd = std(vnow); Z = (pnow-shufmean)/shufstd; 
            if Z>=0 
               P = numel(find(vnow>=pnow))/nn;
            else
               P = numel(find(vnow<=pnow))/nn;
            end
            sig1D.(f1name)(i) = Z; sig1D.(f2name)(i) = P;
         end
      end
    end
end
function P = reformatcleanup(S) %S{nshuf}.(field)(1:nev) to P.(field){nev}(1:nshuf)
nshuf = numel(S); P = []; 
if nshuf>0 
    fname = fieldnames(S{1});
    if ~isempty(fname)
        nev = numel(S{1}.(fname{1}));
        for (i=1:numel(fname))
             P.(fname{i}) = cell(1, nev); Tnow = zeros(nshuf,nev);
             for (j = 1:nshuf)
                 Tnow(j,:) = S{j}.(fname{i});
             end
             for (k = 1:nev)
                 P.(fname{i}){k} = Tnow(:,k);
             end
        end
    end
end

function [evseq, shufevseq, tmpnow, evName, evType, evT, tottime] = sequencedataandshuffle(pinfo, data, behav, bhdata, neuronid, tmpnow, seqparm)
%%%%%first get template and raw data to analyze
cellfilename = tmpnow.tmp.tmpfilename; tmpevname = tmpnow.evtfile;
nlet = numel(cellfilename); isok = zeros(1, nlet); indok = zeros(1, nlet);
for (i = 1:nlet)
    ti = find( strcmp(pinfo.general.parmfile(neuronid), cellfilename{i}) );
    if (numel(ti) == 1)
        isok(i) = 1; indok(i) = ti;
    end
end
iii = find(isok == 1); nlet = numel(iii); cellfilename = cellfilename(iii); cellrate = tmpnow.tmp.tmpspacerate(iii);
Xgrid = tmpnow.tmp.Xgrid; %%%%%%%%%%%%For now: only do it on a 1D space
shufcellrate = shuffltemplate(cellrate, seqparm); %%%This are the shuffled templates: shufcellrate{seqparm.Nshuffle}
indok = indok(iii); spikedata = cell(1, nlet);
for (i = 1:nlet) spikedata{i} = data.spike.spiketime{neuronid(indok(i))}; end %%%%This is the rawdata to analyze
%%%%%%% session sequences/decoding not decoded
sessName = pinfo.general.sessionname{neuronid(1)}; 
sessST = pinfo.general.sessionstartT{neuronid(1)}; sessET = pinfo.general.sessionendT{neuronid(1)}; 
%%%%%%% event sequences/decoding: only on selected events
%%%%%%%%%%%% !!!!!!!!!!!!temporary adjustment if 'modify' option is on: if theta, change to [max max] time
%[evName, evType, evT] = filterevents(pinfo, data, seqparm.evkeyword, seqparm.evkeytype, neuronid, 'modify', tmpevname);
[evName, evType, evT, tottime] = filterevents(pinfo, data, seqparm,  neuronid, 'none'); %%%neuroid contains cells within the same final dir
%%%%%%%%%%%%done temporary adjustment!!!!!!!!!
nev = numel(evName); 
%evsessid = zeros(1, nev); 
evsessName = cell(1, nev); %allevletter = cell(1, nev); allevtime = cell(1, nev); %%They latter two are dummy variables to be compatible with old sequencing methods
evseq = cell(1, nev); shufevseq = cell(seqparm.Nshuffle, nev);
for (i = 1:nev) %%%for each event file
    disp(['-------------------> event now: ' evName{i} '; number of cells found: ' num2str(numel(cellrate))]);
    %%%%work out which session the event is in spikedb
    %[evsessName{i}, evsessid(i)] = identifysession(evT{i}, sessName, sessST, sessET);
    [evsessName{i}, ~] = identifysession(evT{i}, sessName, sessST, sessET);
    evT{i}.start = evT{i}.start - seqparm.setbacktime; evT{i}.ent = evT{i}.ent + seqparm.setbacktime; 
    %%%%work out linearized positions for each lap in the event on the TEMPLATE trajectory
    finaldirnow = pinfo.general.finaldir{neuronid(1)};
    posid = []; rtime= []; rloc = [];
    joint = locatejoint(behav, bhdata, tmpevname, cellfilename); %%% the event (could be different from template - so may not be accurate) linearized on the template trajectory
    if ~isempty(joint)
        if (~isempty(evsessName{i}))
            posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evsessName{i}) );
        end
        if numel(posid) == 1
           postimestamp = bhdata.pos.postimestamp{posid}; 
           Pmarker = behav.parm.sessPmarker{posid}; allposmarker = behav.general.posMarker{posid}; 
           ik = find(strcmp(allposmarker, Pmarker)); 
           if (numel(ik) == 1) 
               posXX = bhdata.pos.XX{posid}{ik}*behav.parm.pixelXSize(posid); posYY = bhdata.pos.YY{posid}{ik}*behav.parm.pixelYSize(posid);
               spaceunit = (behav.parm.pixelXSize(posid) + behav.parm.pixelYSize(posid))/2;
               [rtime, rloc, ~] = linearizenow(joint, evT{i}, postimestamp*behav.parm.timeunit(posid), posXX, posYY, spaceunit); %%%actual running locations
           else
               disp(['-------------------> Warning: position diode for computing actual postions not found: ', Pmarker]);
           end
        end
    end
    %%%%work out the ep structure
    evTime = evT{i}; nep = numel(evTime.start); ep = cell(nep,1); realtime= cell(nep,1); realloc = cell(nep,1);
    for (ti = 1:nep)
         ep{ti,1}.start = evTime.start(ti); ep{ti, 1}.ent = evTime.ent(ti); ep{ti,1}.marker = evTime.marker{ti}; 
         realtime{ti,1} = rtime{ti}; realloc{ti,1} = rloc{ti}; 
    end
    %%%% spikedata within events
    datanow = cell(1, nlet); for (j = 1:nlet) [datanow{j}, ~] = SpikeEventFilter(spikedata{j}, evT{i}); end
    %%%%%%% Output are the raw decoded results
    evseq{i} = finddecoderesults(datanow, cellrate, Xgrid, seqparm, ep, realtime, realloc, evName{i}, 1, 1); %%saveprobdataflag = 1 to save raw prob data
    %%%%%%% compute shuffles
    saveshufrawdata = zeros(1, seqparm.Nshuffle); saveshufquandata = zeros(1, seqparm.Nshuffle);
    if (seqparm.saveshufrawdatarate<=seqparm.Nshuffle) ii = 1:seqparm.saveshufrawdatarate:seqparm.Nshuffle; saveshufrawdata(ii) = ones(size(ii)); end
    if (seqparm.saveshufquandatarate<=seqparm.Nshuffle) ii = 1:seqparm.saveshufquandatarate:seqparm.Nshuffle; saveshufquandata(ii) = ones(size(ii)); end
    for (j = 1:seqparm.Nshuffle)
        shufevseq{j, i} = finddecoderesults(datanow, shufcellrate{j}, Xgrid, seqparm, ep, realtime, realloc, evName{i}, saveshufrawdata(j), saveshufquandata(j)); %%save sample data to control final saved file size (<1GB)
    end
end
function evseq = finddecoderesults(datanow, cellrate, Xgrid, seqparm, ep, realtime, realloc, evName, saverawflag, savequantflag)
[mm,nn] = size(ep); prob = cell(mm,nn); timepoint = cell(mm,nn);
seqstart = cell(mm,nn); seqend = cell(mm,nn); seqmarker = cell(mm,nn); 
matchscore = cell(mm,nn); posmatchprob = cell(mm,nn); negmatchprob = cell(mm,nn); seqtime = cell(mm,nn);
shufZscore = cell(mm,nn); shufZposP = cell(mm,nn); shufZnegP = cell(mm,nn); shufposP = cell(mm,nn); 
shufnegP = cell(mm,nn); decodeerr = cell(mm,nn); nactivecell = cell(mm,nn); nactiveraw = cell(mm,nn); ntarget = cell(mm,nn); 
newquant = cell(mm,nn); quantdata = cell(mm,nn);
for (ii = 1:mm) %%%original lap index
     for (jj = 1:nn)
         [prob{ii,jj}, timepoint{ii,jj}, seqstart{ii,jj}, seqend{ii,jj}, seqmarker{ii,jj}, nactivecell{ii,jj}, nactiveraw{ii,jj}, ntarget{ii,jj}, newquant{ii,jj}, quantdata{ii,jj}, ...
                matchscore{ii,jj}, posmatchprob{ii,jj}, negmatchprob{ii,jj}, ...
                 shufZscore{ii,jj}, shufZposP{ii,jj}, shufZnegP{ii,jj},...
                 shufposP{ii,jj}, shufnegP{ii,jj}, decodeerr{ii,jj}] = ...
                findallpositionprob(datanow, cellrate, Xgrid, seqparm, ep{ii,jj}, realtime{ii,jj}, realloc{ii,jj}); %disp(decodeerr{ii,jj});
     end
end
evseq.seqdata = []; if saverawflag evseq.seqdata = prob; evseq.timepoint = timepoint; evseq.Xgrid = Xgrid; end
evseq.quantdata = []; if savequantflag evseq.quantdata = quantdata; end
%%%%%%%%%%%%%%% empty assignements
evseq.timepeakth = []; evseq.ratestd = []; evseq.spatratemean = []; evseq.spatpeakth = []; evseq.spatratestd = []; evseq.spatratemean = [];
%%%%%%%%%%%%%%%output
evseq.evname = evName; evseq.ep = ep; %evseq.allletter = allevletter; evseq.alltime = allevtime;
evseq.seq = []; evseq.seqstart = seqstart; evseq.seqend = seqend; evseq.seqmarker = seqmarker; evseq.seqtime = seqtime;
evseq.newquant = newquant;
evseq.matchscore = matchscore; evseq.posmatchprob = posmatchprob; evseq.negmatchprob = negmatchprob;
evseq.shufZscore = shufZscore; evseq.shufZposP = shufZposP; evseq.shufZnegP = shufZnegP;
evseq.shufposP = shufposP; evseq.shufnegP = shufnegP; evseq.decodeerr = decodeerr; 
evseq.nactivecell = nactivecell; evseq.nactiveraw = nactiveraw; evseq.ntarget = ntarget;

function [P, candevind, matchevind, posevind, negevind] = findquantifications(evseq, tottime, seqparm, nev)
%%%%%Old parameters: candidate events, significant positive/negative events
%nspl = 3; %%%calculate description statitstics only if this many samples
P.tottime = tottime; P.totalN = NaN*ones(1, nev); P.candN = NaN*ones(1, nev); P.posN = NaN*ones(1, nev); P.negN = NaN*ones(1, nev); P.matchN =NaN*ones(1,nev);
P.evTimeRate = NaN*ones(1, nev); P.candTimeRate = NaN*ones(1, nev); P.posTimeRate = NaN*ones(1, nev); P.negTimeRate = NaN*ones(1, nev); P.matchTimeRate = NaN*ones(1, nev);
P.candEvtRate = NaN*ones(1, nev); P.posEvtRate = NaN*ones(1, nev); P.negEvtRate = NaN*ones(1, nev);  P.matchEvtRate = NaN*ones(1, nev);
P.posRatio = NaN*ones(1, nev); P.negRatio = NaN*ones(1, nev); P.matchRatio = NaN*ones(1, nev);
candevind = cell(1,nev); matchevind=cell(1,nev); posevind=cell(1,nev); negevind=cell(1,nev);
%%%% Added new: mean/median matchscore, mean matchZscore
P.candmeanScore = NaN*ones(1, nev); P.candmedScore = NaN*ones(1, nev); P.candmeanZscore = NaN*ones(1, nev); P.candmedZscore = NaN*ones(1, nev);
P.candmeanAbsScore = NaN*ones(1, nev); P.candmedAbsScore = NaN*ones(1, nev); P.candmeanAbsZscore = NaN*ones(1, nev); P.candmedAbsZscore = NaN*ones(1, nev);
P.posmeanScore = NaN*ones(1, nev); P.posmedScore = NaN*ones(1, nev); P.posmeanZscore = NaN*ones(1, nev); P.posmedZscore = NaN*ones(1, nev);
P.negmeanScore = NaN*ones(1, nev); P.negmedScore = NaN*ones(1, nev); P.negmeanZscore = NaN*ones(1, nev); P.negmedZscore = NaN*ones(1, nev);
P.matchmeanScore = NaN*ones(1, nev); P.matchmedScore = NaN*ones(1, nev); P.matchmeanZscore = NaN*ones(1, nev);  P.matchmedZscore = NaN*ones(1, nev);
P.matchmeanAbsScore = NaN*ones(1, nev); P.matchmedAbsScore = NaN*ones(1, nev); P.matchmeanAbsZscore = NaN*ones(1, nev); P.matchmedAbsZscore = NaN*ones(1, nev);
%%%% Added new: mean/median decodeerr
P.candmeanErr = NaN*ones(1, nev); P.candmedErr = NaN*ones(1, nev); P.posmeanErr = NaN*ones(1, nev); P.posmedErr = NaN*ones(1, nev);
P.negmeanErr = NaN*ones(1, nev); P.negmedErr = NaN*ones(1, nev); P.matchmeanErr = NaN*ones(1, nev); P.matchmedErr = NaN*ones(1, nev);
for (i = 1:nev)
        [mm,nn] = size(evseq{i}.shufposP); shufposP = repack_single(evseq{i}.shufposP); shufnegP = repack_single(evseq{i}.shufnegP);
        nactraw = repack_single(evseq{i}.nactiveraw); nact = repack_single(evseq{i}.nactivecell); 
        candevind{i} = find(nact>=seqparm.minactivecell); 
        posevind{i} = find( (shufposP<seqparm.significancelevel) & (nact>=seqparm.minactivecell)); 
        negevind{i} = find( (shufnegP<seqparm.significancelevel) & (nact>=seqparm.minactivecell));
        matchevind{i} = sort(union(posevind{i}, negevind{i}));
    P.totalN(i) = mm*nn; P.candN(i) = numel(find(nactraw>=seqparm.minactivecell));
    P.posN(i) = numel(posevind{i}); P.negN(i) = numel(negevind{i}); P.matchN(i) = numel(matchevind{i});
    P.evTimeRate(i) = P.totalN(i)/tottime(i); P.candTimeRate(i) = P.candN(i)/tottime(i); 
    P.posTimeRate(i) = P.posN(i)/tottime(i); P.negTimeRate(i) = P.negN(i)/tottime(i); P.matchTimeRate(i) = P.matchN(i)/tottime(i);
    P.candEvtRate(i) = P.candN(i)/P.totalN(i); P.posEvtRate(i) = P.posN(i)/P.totalN(i); P.negEvtRate(i) = P.negN(i)/P.totalN(i); P.matchEvtRate(i) = P.matchN(i)/P.totalN(i);
    P.posRatio(i) = P.posN(i)/P.candN(i); P.negRatio(i) = P.negN(i)/P.candN(i);  P.matchRatio(i) = P.matchN(i)/P.candN(i);    
        matchscore = repack_single(evseq{i}.matchscore); shufZscore = repack_single(evseq{i}.shufZscore); err = repack_single(evseq{i}.decodeerr);       
    kk = candevind{i}; 
    if ~isempty(kk)
       T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.candmeanScore(i)=mean(T); P.candmedScore(i)=median(T); P.candmeanAbsScore(i)=mean(abs(T)); P.candmedAbsScore(i)=median(abs(T));  end
       T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.candmeanZscore(i)=mean(T); P.candmedZscore(i)=median(T);  P.candmeanAbsZscore(i)=mean(abs(T)); P.candmedAbsZscore(i)=median(abs(T)); end
       T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.candmeanErr(i)=mean(T); P.candmedErr(i)=median(T); end
    end
    kk = matchevind{i}; 
    if ~isempty(kk)
       T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.matchmeanScore(i)=mean(T); P.matchmedScore(i)=median(T); P.matchmeanAbsScore(i)=mean(abs(T)); P.matchmedAbsScore(i)=median(abs(T));  end
       T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.matchmeanZscore(i)=mean(T); P.matchmedZscore(i)=median(T);  P.matchmeanAbsZscore(i)=mean(abs(T)); P.matchmedAbsZscore(i)=median(abs(T)); end
       T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.matchmeanErr(i)=mean(T); P.matchmedErr(i)=median(T); end
    end
    kk = posevind{i}; 
    if ~isempty(kk)
       T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.posmeanScore(i)=mean(T); P.posmedScore(i)=median(T);   end
       T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.posmeanZscore(i)=mean(T); P.posmedZscore(i)=median(T); end   
       T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.posmeanErr(i)=mean(T); P.posmedErr(i)=median(T); end
    end
    kk = negevind{i}; 
    if ~isempty(kk)
       T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.negmeanScore(i)=mean(T); P.negmedScore(i)=median(T);   end
       T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.negmeanZscore(i)=mean(T); P.negmedZscore(i)=median(T); end    
       T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.negmeanErr(i)=mean(T); P.negmedErr(i)=median(T); end
    end
end
%%%% Added new: mean/median steps and otehrs
%%%%% All cand events stat: Nstep, Ngstep, PProb, EndDis, Length, MedianStep, MaxStep, MinStep, MeanStep, StepVar
P.candNstep = NaN*ones(1, nev); P.candNDstep = NaN*ones(1, nev); P.candNVstep = NaN*ones(1, nev); P.candmedPProb = NaN*ones(1, nev); P.candmeanPProb = NaN*ones(1, nev);
P.candmedLength = NaN*ones(1, nev); P.candmeanLength = NaN*ones(1, nev); P.candmedEndDis = NaN*ones(1, nev); P.candmeanEndDis = NaN*ones(1, nev); 
P.candmeanStep = NaN*ones(1, nev); P.candstdStep = NaN*ones(1, nev); P.candmedStep = NaN*ones(1, nev);
P.candmaxStep = NaN*ones(1, nev); P.candminStep = NaN*ones(1, nev); 
%%%%% All match events stat: Nstep, Ngstep, PProb, EndDis, Length, MedianStep, MaxStep, MinStep, MeanStep, StepVar
P.matchNstep = NaN*ones(1, nev); P.matchNDstep = NaN*ones(1, nev); P.matchNVstep = NaN*ones(1, nev); P.matchmedPProb = NaN*ones(1, nev); P.matchmeanPProb = NaN*ones(1, nev); 
P.matchmedLength = NaN*ones(1, nev); P.matchmeanLength = NaN*ones(1, nev); P.matchmedEndDis = NaN*ones(1, nev); P.matchmeanEndDis = NaN*ones(1, nev); 
P.matchmeanStep = NaN*ones(1, nev); P.matchstdStep = NaN*ones(1, nev); P.matchmedStep = NaN*ones(1, nev);
P.matchmaxStep = NaN*ones(1, nev); P.matchminStep = NaN*ones(1, nev);
for (i=1:nev) 
    %%%candidate events
    kk = candevind{i}; prop = repackout(evseq{i}.newquant); 
    if ~isempty(kk)
       NN = prop.Nstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candNstep(i) = median(NN(ii)); end
       NN = prop.NVstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candNVstep(i) = median(NN(ii)); end %%%number of Valid steps
       NN = prop.NDstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candNDstep(i) = median(NN(ii)); end %%%number of decoded positions
       NN = prop.MedianPProb(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmedPProb(i) = median(NN(ii)); end
       NN = prop.MeanPProb(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmeanPProb(i) = mean(NN(ii)); end
       %disp(['Prob: ' num2str([numel(ii) median(NN(ii)) mean(NN(ii)) std(NN(ii))])]);
       NN = prop.Length(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmedLength(i) = median(NN(ii)); P.candmeanLength(i) = mean(NN(ii)); end 
       %disp(['Length: ' num2str([numel(ii) median(NN(ii)) mean(NN(ii)) std(NN(ii))])]);
       NN = prop.EndDis(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmedEndDis(i) = median(NN(ii)); P.candmeanEndDis(i) = mean(NN(ii)); end
       NN = prop.MeanStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmeanStep(i) = mean(NN(ii)); end
       NN = prop.StepStd(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candstdStep(i) = mean(NN(ii)); end
       NN = prop.MedianStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmedStep(i) = median(NN(ii)); end
       NN = prop.MaxStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmaxStep(i) = median(NN(ii)); end
       NN = prop.MinStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candminStep(i) = median(NN(ii)); end
    end
    %%%match events
    kk = matchevind{i};
    if ~isempty(kk)
       NN = prop.Nstep(kk); P.matchNstep(i) = median(NN); 
       NN = prop.NVstep(kk); P.matchNVstep(i) = median(NN); %%%number of Valid steps
       NN = prop.NDstep(kk); P.matchNDstep(i) = median(NN); %%%number of decoded positions      
       NN = prop.MedianPProb(kk); P.matchmedPProb(i) = median(NN); NN = prop.MeanPProb(kk); P.matchmeanPProb(i) = mean(NN);
       NN = prop.Length(kk); P.matchmedLength(i) = median(NN); P.matchmeanLength(i) = mean(NN);
       NN = prop.EndDis(kk); P.matchmedEndDis(i) = median(NN); P.matchmeanEndDis(i) = mean(NN);
       NN = prop.MeanStep(kk); P.matchmeanStep(i) = mean(NN); NN = prop.StepStd(kk); P.matchstdStep(i) = mean(NN);
       NN = prop.MedianStep(kk); P.matchmedStep(i) = median(NN);
       NN = prop.MaxStep(kk); P.matchmaxStep(i) = median(NN); NN = prop.MinStep(kk); P.matchminStep(i) = median(NN); 
    end
end

function [prob, timepoint, seqstart, seqend, seqmarker, nactivecell, nactiveraw, ntarget, newquant, quantdata, ...
    matchscore, posmatchprob, negmatchprob, shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decodeerr] = findallpositionprob(spiketime, cellrate, Xgrid, seqparm, ep, realtime, realloc)
decparm.winsize = seqparm.windowsize; decparm.shifttime = seqparm.shifttime;
[timepoint, prob] = DE_Seq_DecodePositionProbability(spiketime, Xgrid, cellrate, ep, decparm);
%seqstart = ep.start; seqend = ep.ent; seqmarker = ep.marker;
%%%%%%%% quantification and shuffle significance for individual events
[newquant, quantdata, ntarget, seqstart, seqend, seqmarker, ...
    matchscore, posmatchprob, negmatchprob, shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decodeerr] = finddecodingquantifications(timepoint, prob, Xgrid, seqparm, realtime, realloc);
nactiveraw = findactivecell(spiketime, ep); nactivecell = mean(nactiveraw)*ones(size(matchscore)); %%%%Be careful here: assuming nep = 1, but matchscore could be 1 (WithinEvents)or multiple (FreeSequences)
function [prop, propdata, ntarget, seqstart, seqend, seqmarker, ...
    matchscore, posmatchprob, negmatchprob, shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decodeerr] = finddecodingquantifications(timepoint, prob, Xgrid, seqparm, alltime, allloc)
%%%%%timepoint{nep}(ntimepoint); prob{nep}{ntimepoint}(nXgrid); Xgrid(nXgrid) 
%%%%%%%%%%%determine the decoded position, then compute: decode error; pearson R, p; shuffle Z, p 
nep = numel(timepoint); %%%% this is for single lap at this point
Xgrid = Xgrid{1}; %%%%for 1D space
%%%%%%%%%%%%%%%determine decoded peak positions: apply thresholding if necessary
%propdata.timepoint = cell(1, nep); propdata.peakloc = cell(1, nep); propdata.peakprob = cell(1, nep); propdata.stepdis = cell(1, nep);
peakloc = cell(1, nep); peakprob = cell(1, nep); stepdis = cell(1, nep);
for (i = 1:nep) %%%be careful about peakprob and maxstep filtering - sequence is important: should calculate first then filter -- the following is OK!
    timenow = timepoint{i}; ntime = numel(timenow); peakx = NaN*ones(size(timenow)); peakp = NaN*ones(size(timenow)); dis = NaN*ones(size(timenow));
    for (j = 1:ntime) 
        %disp( [size(prob{i}{j}) size(Xgrid)] );
        [mm, iii] = max(prob{i}{j}); 
        if mm > seqparm.probthres/numel(find(~isnan(prob{i}{j})))  %%%if propthres exists (default 0 = no exclusion)
            peakx(j) = Xgrid(iii); peakp(j) = mm; 
        end
    end
    for (j = 2:ntime)
        if (~isnan(peakx(j-1))) && (~isnan(peakx(j)))
            dis(j) = abs(peakx(j)-peakx(j-1));
        end
    end
    %%%%% The following is not implemented, beacause it is hard to decide which of the two neighbors to exclude if distance is abnornmally high 
%     iii = find(dis>seqparm.stepdisthres); %%% if any max jump distance to select peak locations: 
%     if ~isempty(iii)
%         dis(iii) = NaN*ones(size(iii)); peakx(iii) = NaN*ones(size(iii)); peakp(iii) = NaN*ones(size(iii));
%     end
    peakloc{i} = peakx; peakprob{i} = peakp; stepdis{i} = dis; %%%% Currently this is for a single lap (event)
end
if strcmp(seqparm.eventoption, 'FreeSequences') %%%%%%If FreeSequences - decide sequences within the event - multiple outputs
   [seqstartind, seqendind, Pout] = findfreematchingsequences(peakloc{nep}, seqparm); %%%%Be careful here: assuming nep = 1;
   ntarget = numel(seqstartind); disnow = cell(1, ntarget); peaklocnow = cell(1, ntarget);
   peakprobnow = cell(1, ntarget); timepointnow = cell(1, ntarget); 
   for (ttk = 1:ntarget)
       disnow{ttk} = stepdis{nep}(seqstartind(ttk):seqendind(ttk)); 
       peaklocnow{ttk} = peakloc{nep}(seqstartind(ttk):seqendind(ttk));
       peakprobnow{ttk} = peakprob{nep}(seqstartind(ttk):seqendind(ttk)); 
       timepointnow{ttk} = timepoint{nep}(seqstartind(ttk):seqendind(ttk));
   end
elseif strcmp(seqparm.eventoption, 'WithinEvents')%%%%%%If WithinEvent - take the whole event as a single output; 
   ntarget = 1; disnow = stepdis; peaklocnow = peakloc; peakprobnow = peakprob; timepointnow = timepoint; %%%already in (single) cell array
end
%%%%%%% new version properties
prop = findnewproperties(disnow, peaklocnow, peakprobnow, timepointnow);
%%%old variables
[matchscore, posmatchprob, negmatchprob, shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decodeerr, seqstart, seqend, seqmarker] =...
        findoldproperties(peaklocnow, timepointnow, alltime, allloc, seqparm);
%%%%%% if free sequencing, need to re-asign the scores and probabilities 
if strcmp(seqparm.eventoption, 'FreeSequences')
   matchscore = Pout.matchscore; posmatchprob = Pout.posmatchprob; negmatchprob = Pout.negmatchprob;
   shufZscore = NaN*ones(1, ntarget); shufZposP = NaN*ones(1, ntarget); shufZnegP = NaN*ones(1, ntarget);
   shufposP = posmatchprob; shufnegP = negmatchprob; 
end
propdata.timepoint = timepoint; propdata.peakloc = peakloc; propdata.peakprob = peakprob; propdata.stepdis = stepdis; 

function P = findnewproperties(dis, peakx, peakprob, timepoint)
nep = numel(timepoint); P.Nstep = NaN*ones(1, nep); P.NDstep = NaN*ones(1, nep);  P.NVstep = NaN*ones(1, nep); 
P.MedianPProb = NaN*ones(1, nep); P.MeanPProb = NaN*ones(1, nep);
P.EndDis = NaN*ones(1, nep); P.Length = NaN*ones(1, nep);
P.MedianStep = NaN*ones(1, nep); P.MaxStep = NaN*ones(1, nep); P.MinStep = NaN*ones(1, nep); 
P.MeanStep = NaN*ones(1, nep); P.StepStd = NaN*ones(1, nep);  
for (i = 1:nep)
    P.Nstep(i) = numel(timepoint{i}); 
    iii = find(~isnan(peakx{i})); if ~isempty(iii) peaknow = peakx{i}(iii); P.EndDis(i) = peaknow(numel(peaknow)) - peaknow(1); end 
    iii = find(~isnan(dis{i})); P.NVstep(i) = numel(iii);
    if (~isempty(iii))
        disnow = dis{i}(iii); P.MedianStep(i) = median(disnow); P.MinStep(i) = min(disnow); P.MaxStep(i) = max(disnow);
        P.MeanStep(i) = mean(disnow); P.StepStd(i) = std(disnow); P.Length(i) = sum(disnow);
    end
    iii = find(~isnan(peakprob{i})); 
    if (~isempty(iii))
         P.NDstep(i) = numel(iii); pnow = peakprob{i}(iii); P.MedianPProb(i) = median(pnow); P.MeanPProb(i) = mean(pnow);
    end
end
function [matchscore, posmatchprob, negmatchprob, shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decodeerr, seqstart, seqend, seqmarker] =...
    findoldproperties(peakloc, timepoint, alltime, allloc, seqparm)
nshuffle = seqparm.Nwithineventshuffle; windowsize = seqparm.windowsize;
nep = numel(timepoint); matchscore = NaN*ones(1, nep); posmatchprob = NaN*ones(1, nep); negmatchprob = NaN*ones(1, nep); decodeerr = NaN*ones(1, nep);
shufZscore = NaN*ones(1, nep); shufZposP = NaN*ones(1, nep); shufZnegP = NaN*ones(1, nep); shufposP = NaN*ones(1, nep); shufnegP = NaN*ones(1, nep);
seqstart = NaN*ones(1, nep); seqend = NaN*ones(1, nep); seqmarker = cell(1, nep);
for (i = 1:nep)
    iii = find(~isnan(peakloc{i}));
    if ~isempty(iii)
       timenow = timepoint{i}(iii); locnow = peakloc{i}(iii); ntime = numel(timenow);
       [timenow, iii] = sort(timenow); locnow = locnow(iii); 
       seqstart(i) = timenow(1); seqend(i) = timenow(numel(timenow)); %%%seqmarker{i} = 'decoding1D';
       %%%Pearlson R, Ps
       [rr, pp] = corrcoef(timenow, locnow); %disp(rr)
       if numel(rr) > 1 %%%if timenow or locnow is empty, rr is 1, not a 2x2 matrix
          matchscore(i) = rr(1,2);
          if matchscore(i) >= 0
             posmatchprob(i) = pp(1,2); 
          else
             negmatchprob(i) = pp(1,2); 
          end
       end
       %%%errors
       if ~isempty(alltime)
          realloc = NaN*ones(size(locnow)); %disp(timegrid'); disp(alltime);
          for (ii = 1:ntime)
              jj = find( (alltime>=timenow(ii)-windowsize/2) & (alltime<timenow(ii)+windowsize/2) );
              if (numel(jj) >= 1) 
                  aaa = allloc(jj); aaa = aaa(~isnan(aaa));
                  realloc(ii) = mean(aaa); 
              end
          end
          dis = abs(locnow-realloc); 
          decodeerr(i) = median(dis(~isnan(dis)));
       end
       %%%shuffle Z (Zp), Sp
       if nshuffle > 1
          shuffleR = zeros(1, nshuffle);
          for (ii = 1:nshuffle)
            %C = corrcoef(timenow,locnow(randperm(ntime))); shuffleR(i) = C(1,2);
            shuffleR(ii) = Utilities_FindCorrCoef_self(timenow,locnow(randperm(ntime))); %%%both inputs are row vectors
          end
          ss = std(shuffleR);
          if (ss ~= 0)
           shufZscore(i) = (matchscore(i) - mean(shuffleR))/ss; 
          else
           shufZscore(i) = NaN;
          end
          ipsi = 0.00000001;
          if (shufZscore(i) >= 0)
           shufZposP(i) = 1- normcdf(shufZscore(i), 0, 1); 
           shufposP(i) = numel(find(shuffleR>= (matchscore(i)-ipsi) ))/nshuffle;
          elseif (shufZscore(i) < 0)
           shufZnegP(i) = normcdf(shufZscore(i), 0, 1); 
           shufnegP(i) = numel(find(shuffleR<= (matchscore(i)+ipsi) ))/nshuffle;
          end 
       end
    end
end
function [startindout, endindout, Pout] = findfreematchingsequences(peakloc, seqparm)
%%%%Identify all possible significant sequences (rank order determined), then resolve the overlap (take the most signidicant ones)
peakind = find(~isnan(peakloc)); peaknow = peakloc(peakind); nlet = numel(peakind); 
maxlet = min([nlet seqparm.freeMaxTimepoints]); minlet = max([4 seqparm.freeMinTimepoints]);
seqleng = [maxlet:-1:minlet];
istart = cell(1, numel(seqleng)); iend = cell(1, numel(seqleng)); 
score = cell(1, numel(seqleng)); mscore = cell(1, numel(seqleng));
for (i=1:numel(seqleng))
    lengnow = seqleng(i); startind = [1:1:nlet-lengnow+1]; endind = [lengnow:1:nlet]; %nseqnow = nlet-lengnow + 1; 
    [pscore, matchscore] = findseqscoreprob_multiple(peaknow, startind, endind); %, seqparm.Nwithineventshuffle);
    tt = find(pscore < seqparm.significancelevel); score{i} = pscore(tt); %%%score is significance: the lower the better
    istart{i} = startind(tt); iend{i} = endind(tt); mscore{i} = matchscore(tt);
end
istart = cell2mat(istart); iend = cell2mat(iend); score = cell2mat(score); mscore = cell2mat(mscore); %%%already sorted from long to short
%[startind, endind, pscore, ind] = resolveoverlap_pvaltakeover(score, istart, iend); %%%% -------------------------This or below needs to settle for specific problems--------------------------
[startind, endind, pscore, ind] = resolveoverlap_lengtakeover(score, istart, iend);
Pout.matchscore = mscore(ind); Pout.posmatchprob = NaN*ones(size(ind)); Pout.negmatchprob = NaN*ones(size(ind));
for (i = 1:numel(startind))
    if Pout.matchscore(i) >= 0
       Pout.posmatchprob(i) = pscore(i); 
    else
       Pout.negmatchprob(i) = pscore(i); 
    end
end
startindout = peakind(startind); endindout = peakind(endind); 
function [pscore, matchscore] = findseqscoreprob_multiple(peaknow, startind, endind) % nshuffle) 
%randorder or linear correlation?
pscore = NaN*ones(1, numel(startind)); matchscore = NaN*ones(1, numel(startind));
for (i = 1:numel(startind))
    timerank = 1:(endind(i)-startind(i)+1); 
    %peakrank=peaknow(startind(i):endind(i));   %%%linear correlaton -------------------------This or below need to settle for specific problems--------------------------
    [~, peakrank]=sort(peaknow(startind(i):endind(i)));   %Rankorder is too liberal for sequence detection:  
    [rr, pp] = corrcoef(timerank, peakrank); %%%both timerank and peakrank need to have at least two elements; otherwise numel(rr) and numel(pp) =1
    matchscore(i) = rr(1,2); pscore(i) = pp(1,2); %%%this seems significantly faster than the option below
    %[matchscore(i), pscore(i)] = corr(timerank', peakrank', 'Type', 'Spearman'); %%%use spearman correlation - shuffle p still will not solve the issue of sig corr between [1 2 3 4...] and [0 0 1 1...]
end
function [newstart, newend, newscore, ind] = resolveoverlap_lengtakeover(score, sstart, send)
ist = false(1, numel(score)); %%%set as if a target sequence
leng = send-sstart+1;
for (i = 1:numel(score))
    iii = findoverlap(sstart(ist), send(ist), sstart(i), send(i)); 
    if isempty(iii)
        ist(i) = 1;
    else
        kkk = find(ist);
        if leng(i)>=max(leng(ist(kkk(iii))))&& score(i)<min(score(ist(kkk(iii))))
            ist(i) = 1; ist(kkk(iii)) = zeros(1, numel(iii));
        end
    end
end
newstart = sstart(ist); newend = send(ist); newscore = score(ist); ind = find(ist);
[newstart, kkk] = sort(newstart); newend = newend(kkk); newscore = newscore(kkk); ind = ind(kkk);
function [newstart, newend, newscore, ind] = resolveoverlap_pvaltakeover(score, sstart, send)
ist = false(1, numel(score)); %%%set as if a target sequence
for (i = 1:numel(score))
    iii = findoverlap(sstart(ist), send(ist), sstart(i), send(i)); 
    if isempty(iii)
        ist(i) = 1;
    else
        kkk = find(ist);
        if score(i)<min(score(ist(kkk(iii))))
            ist(i) = 1; ist(kkk(iii)) = zeros(1, numel(iii));
        end
    end
end
newstart = sstart(ist); newend = send(ist); newscore = score(ist); ind = find(ist);
[newstart, kkk] = sort(newstart); newend = newend(kkk); newscore = newscore(kkk); ind = ind(kkk);
function ind = findoverlap(newstart, newend, sigstart, sigend)
ind = [];
iii = find( (newstart<sigstart) & (newend>sigstart) ); ind = union(ind, iii);
iii = find( (newstart<sigend) & (newend>sigend) ); ind = union(ind, iii);
iii = find( (newstart>=sigstart) & (newend<=sigend) ); ind = union(ind, iii);

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

function joint = locatejoint(behav, bhdata, evname, tmpfilename)
joint = [];
%%%get the final directory
tind = strfind(tmpfilename{1}, 'final'); tmpfdir = tmpfilename{1}(1:tind+4);
sessid = []; evid = []; nsess = numel(behav.general.finaldir);
for (i = 1:nsess)
    evnamenow = behav.general.eventname{i}; finaldirnow = behav.general.finaldir{i};
    jjj = find(strcmp(evnamenow, evname));
    if (numel(jjj)==1) && strcmp(finaldirnow, tmpfdir)
        sessid = i; evid = jjj; break
    end
end
if ~isempty(sessid)
    posltrfile = bhdata.pos.ltrfilename{sessid}; posltr = bhdata.pos.posltr{sessid};
    evPosltr = behav.parm.eventPosltr{sessid};
    for (tt = 1:numel(posltrfile))
        if (strcmp(posltrfile{tt}, evPosltr{evid}))
           joint = posltr{tt}; 
           joint(:,1) = joint(:,1)*behav.parm.pixelXSize(sessid); joint(:,2) = joint(:,2)*behav.parm.pixelYSize(sessid);
           break
        end
    end
end

function [realtime, reallocx, pout] = linearizenow(joint, evT, postimestamp, xpos, ypos, spaceunit)
nep = numel(evT.start); reallocx = cell(1, nep); realtime = cell(1, nep); pout = [];
for (tt = 1:nep)
    evposind = find( (postimestamp>=evT.start(tt)) & (postimestamp<=evT.ent(tt)) );
    realtime{tt} = postimestamp(evposind); reallocx{tt} = NaN*ones(size(realtime{tt}));
    if ~isempty(evposind)
          [xout, ~, pout] = LinearizePathXY_new(joint, xpos(evposind), ypos(evposind), 1, 0, spaceunit);
          reallocx{tt} = xout; 
    end
end
function shufcellrate = shuffltemplate(cellrate, seqparm)
nshuffle = seqparm.Nshuffle; shufcellrate = cell(1,nshuffle);
ncell = numel(cellrate); 
if ncell>0
   nbin = numel(cellrate{1});
   if strcmp(seqparm.shufflemode, 'GlobalID')
      for (i = 1:nshuffle) shufcellrate{i} = cellrate(randperm(ncell)); end
   elseif strcmp(seqparm.shufflemode, 'CircularSlide')
      for (i = 1:nshuffle)
          shufcellrate{i} = cellrate{i};
          for (j = 1:ncell)
              IX = mod((1:nbin)+ ceil(rand*nbin), nbin) + 1; shufcellrate{i}{j} = cellrate{i}{j}(IX);
          end
      end
   else
      disp('----------> can not generate shuffled data or shuffle mode not available; shuffle significance not computed');
   end
end

function pp = repackout(prop)  %%prop{nep, 1}.(field)(1 or ntarget) to pp.(field)(1:nep*ntarget) 
[mm,nn]=size(prop); pp = [];
if (mm*nn>=1)
   if (nn~=1)
       disp(['--------------------warning: unexpected dimensions encountered!']);
   end
   allf = fieldnames(prop{1,1});
   for (k = 1:numel(allf))
       ppnow = cell(1, mm*nn); %ppnow = NaN*ones(1, mm*nn);
       for (i = 1:mm)
       for (j = 1:nn)
           %ppnow((i-1)*nn+j) = prop{i,j}.(allf{k})(1);
           ppnow{(i-1)*nn+j} = prop{i,j}.(allf{k});
       end
       end
       pp.(allf{k}) = cell2mat(ppnow); %pp.(allf{k}) = ppnow;
   end
end
function pp = repack_single(Fnow) %%[prop.(singlefield)]{nlap, 1}.(1 or ntarget) to pp(nep*ntarget) 
[mm,nn]=size(Fnow); ppnow = cell(1, mm*nn);
if (nn~=1)
    disp(['--------------------warning: unexpected dimensions encountered!']);
end
for (i = 1:mm)
    for j = 1:nn
        ppnow{(i-1)*nn+j} = Fnow{i,j};
    end
end
pp = cell2mat(ppnow);
% function pp = repack_oldvar(prop) %%DOES NOT WORK YET -since seqdata contains more complex arrays- prop.(field){nlap, 1}.(1 or ntarget) to pp.(field)(nep*ntarget) 
% allf = fieldnames(prop);
% for (k = 1:numel(allf))
%     disp(allf{k});
%     Fnow = prop.(allf{k}); [mm,nn]=size(Fnow); ppnow = cell(1, mm*nn);
%     if (nn~=1)
%        disp(['--------------------warning: unexpected dimensions encountered!']);
%     end
%     for (i = 1:mm)
%     for j = 1:nn
%         ppnow{(i-1)*nn+j} = Fnow{mm,nn};
%     end
%     end
%     pp.(allf{k}) = cell2mat(ppnow);
% end

function nact = findactivecell(spiketime, ep)
nep = numel(ep.start); ncell = numel(spiketime); cellact = zeros(ncell, nep); 
for (i = 1:nep)
    for (j = 1:ncell)
        if numel(find( (spiketime{j}>=ep.start(i)) & (spiketime{j}<=ep.ent(i)) )) > 0
            cellact(j,i) = 1;
        end
    end
end
nact = sum(cellact);
         
function tottime = findtottime(evNames, evTimes, pinfo, data, cellid, seqparm)
%%%How to compute total time is tricky: user selected for now
%%%options: session - use entire session
nev = numel(evNames); tottime = zeros(1, nev); %%%%%for each event file
for (i = 1:nev)
    %%%find event total time: below is to set the event times back the original value for session identification
    evTimes{i}.start = evTimes{i}.start - seqparm.setbacktime; evTimes{i}.ent = evTimes{i}.ent + seqparm.setbacktime;
    [~, evsessid] = identifysession(evTimes{i}, pinfo.general.sessionname{cellid}, pinfo.general.sessionstartT{cellid}, pinfo.general.sessionendT{cellid});%%%%work 
    if isempty(evsessid)
       disp(['-----------------Warning: corresponding session not found for the event: ', evNames{i}]); 
    else
        sessnow = pinfo.general.sessionname{cellid}{evsessid};
        switch seqparm.tottimeoption
            case 'session'
                tottime(i) = pinfo.general.sessionlength{cellid}(evsessid);
            case 'events'
                evkeyW = seqparm.tottimekeyword; evkeyT = seqparm.tottimekeytype;
                allevnames = pinfo.general.eventname{cellid}; allevtypes = pinfo.parm.eventtype{cellid}; 
                evsel = checkifcorrectevents(allevnames, allevtypes, evkeyW, evkeyT, [], []); %%%if fullfil tottime event criteria
                evpos = find(evsel == 1);
                for (j=1:numel(evpos))
                    if contains(lower(allevnames{evpos(j)}), lower(sessnow)) %%%if belong to the session
                       evtimenow =  data.events.eventtimes{cellid}{evpos(j)};
                       tottime(i) = tottime(i) + sum(evtimenow.ent - evtimenow.start);
                    end
                end
        end
    end
end
function [evSess, evsessid] = identifysession(evTime, sessionname, startT, endT)
evSess = []; evsessid = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii}; evsessid = iii;
end

function [nall, npos, nneg, ncand, posevind, negevind, candevind] = findmatchnum(nactivecell, posP, negP, seqparm)
%function [npos, posP, posmm, nneg, negP, negmm, posSig, negSig] = findtheorysig(evseq, slevel, rrank, eventoption)
%%%%%%here seq{i,j} corresponds to ep{i,j}
%npos = 0; nneg = 0; ncand = 0;
[mm,nn] = size(posP); shufposP = reshape(cell2mat(posP), mm*nn, 1); nall = numel(shufposP);
shufnegP = reshape(cell2mat(negP), mm*nn, 1); nact = reshape(cell2mat(nactivecell), mm*nn, 1);
candevind = find(nact>=seqparm.minactivecell); ncand = numel(candevind);
posevind = numel(find( (shufposP<seqparm.significancelevel) & (nact>=seqparm.minactivecell)) ); npos = numel(posevind);
negevind = numel(find( (shufnegP<seqparm.significancelevel) & (nact>=seqparm.minactivecell)) ); nneg = numel(negevind);

function [evName, evType, evT, tottime] = filterevents(pinfo, data, seqparm, neuronid, ifmod)
evkeyword = seqparm.evkeyword; evkeytype = seqparm.evkeytype; evkeynoword = seqparm.evkeynoword; evkeynotype = seqparm.evkeynotype;
evName = pinfo.general.eventname{neuronid(1)}; %%%neuronid here contains only cells in the same final dir
evType = pinfo.parm.eventtype{neuronid(1)};
evT = data.events.eventtimes{neuronid(1)}; 
%%%%%%%%%% !!!!modify theta events to [max max] 
nev = numel(evName);
if strcmp(ifmod, 'modify')
  for (i = 1:nev)
    if ~isempty(strfind(evName{i}, 'theta'))
        nlap = numel(evT{i}.start); 
        %startT = evT{i}.ent(1:nlap-1); endT = evT{i}.ent(2:nlap); %%Max to Max
        startT = evT{i}.start(1:nlap-1); endT = evT{i}.start(2:nlap); %%Min to Min
        iii = find(endT-startT<0.2); startT = startT(iii); endT = endT(iii);
        evT{i}.start = startT(iii); evT{i}.ent = endT(iii);
        %disp(['--------------> Warning: theta event times adjusted [Maxend(N) Maxend(N+1)]']);
        disp(['--------------> Warning: theta event times adjusted [Minend(N) Minend(N+1)]']);
    end
  end
end
%%%%%%%%%%!!!!done modification
evsel = checkifcorrectevents(evName, evType, evkeyword, evkeytype, evkeynoword, evkeynotype);
evpos = find(evsel == 1);
evName = evName(evpos); evType = evType(evpos); evT = evT(evpos);
tottime = findtottime(evName, evT, pinfo, data, neuronid(1), seqparm);
function evsel = checkifcorrectevents(evName, evType, evkeyword, evkeytype, evkeynoword, evkeynotype)
evsel = ones(1,numel(evName));
for (i = 1:numel(evName))
    if (~isempty(evkeytype))||(~isempty(evkeyword)) %%%for inclusion
      if isempty(evkeytype)
         if ~contains(lower(evName{i}), lower(evkeyword)) evsel(i) = 0; end 
      elseif isempty(evkeyword)
         if ~strncmpi(evType{i}, evkeytype, 3) evsel(i) = 0; end 
      else
         if ~(contains(lower(evName{i}), lower(evkeyword)) && strncmpi(evType{i}, evkeytype, 3))
             evsel(i) = 0; 
         end 
      end
    end
    if (~isempty(evkeynotype))||(~isempty(evkeynoword)) %%%for exclusion
       if isempty(evkeynotype)
          if contains(lower(evName{i}), lower(evkeynoword)) evsel(i) = 0; end 
       elseif isempty(evkeynoword)
          if strncmpi(evType{i}, evkeynotype, 3) evsel(i) = 0; end
       else
          if contains(lower(evName{i}), lower(evkeynoword)) || strncmpi(evType{i}, evkeynotype, 3)
              evsel(i) = 0; 
          end   
       end
    end
end
