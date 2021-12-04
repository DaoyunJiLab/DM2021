function DM_Seq_2DBayesianDecoding_events_Callback
%%Do Bayesian decoding analysis on a database (.spikedb), require a matching template.tmp file
%%This one DECODEs A 2D SPACE!!
%%%%%%%%%%%%%%%ONLY DO IT WITHIN EVENTS SPECIFIED BY KEYWORD/KEYTYPE%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%work on selected cells (groups), group cells into days, and work on each day
%%%%%%%% find if templates exist and then work on each template (if no template found for a date, warn and skip)
%%%%%%%%%%%%%%%%%%%for each event: output decoded prob and replay quantifications
%%%%event parsing options: within events
hf = gcbf; hgroup = getappdata(hf, 'hgroup'); ok = 1; disp('-----> 2D Bayesian decoding: using templates as items');
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
    seqparm.seqtype = 'Bayesian2D'; seqparm.eventoption = 'WithinEvents';
    input = inputdlg({'Enter event keyword'; 'Enter event type'; 'Enter event keyNOword'; 'Enter event NOtype'}, 'Event selection', 4, {'track'; 'run'; 'first'; ''}); 
    if (~isempty(input))
        seqparm.evkeyword = input{1}; seqparm.evkeytype = input{2}; seqparm.evkeynoword = input{3}; seqparm.evkeynotype = input{4};
    else
        ok = 0;
    end
end
if ok
    %input = inputdlg({'Min peakprob allowed(x)'; 'Max stepdis allowed (cm)'; 'Target: min median peakprob'; 'Target: [min max] continuous steps'; 'Target: [min max] stepdis (cm)'},...
    %    'Decoding quality criteria', 5, {'0'; '1000000'; '0'; '5 1000'; '0 50'}); 
    input = inputdlg({'Min peakprob allowed(x avg)'; 'Target: min median peakprob'; 'Target: [min max] continuous steps'; 'Target: [min max] stepdis (cm)'; 'Target: [min max] length (cm)'},...
        'Quality/target criteria', 5, {'1'; '0'; '5 1000'; '0 50'; '0 5000'}); 
    if (~isempty(input))
       seqparm.probthres = str2num(input{1}); seqparm.minmedpprob = str2num(input{2}); %seqparm.stepdisthres = str2num(input{2}); 
       aa = str2num(input{3}); if numel(aa)==2 seqparm.minNstep = aa(1); seqparm.maxNstep = aa(2); else ok = 0; disp('--------> invalid paramer values; aborted.'); end
       aa = str2num(input{4}); if numel(aa)==2 seqparm.minstepdis = aa(1); seqparm.maxstepdis = aa(2); else ok = 0; disp('--------> invalid paramer values; aborted.'); end
       aa = str2num(input{5}); if numel(aa)==2 seqparm.minlength = aa(1); seqparm.maxlength = aa(2); else ok = 0; disp('--------> invalid paramer values; aborted.'); end
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
    input = inputdlg({'Min active cell #'; 'Window size (s)'; 'Window shift time (s)'; 'Setback time (s)'},...
        'Decoding parameters', 4,  {'4'; '0.3'; '0.1'; '0.02'}); % {'4'; '0.3'; '0.03'; '0'}); 
    if (~isempty(input))
       seqparm.minactivecell = str2num(input{1});
       seqparm.windowsize = str2num(input{2}); seqparm.shifttime = str2num(input{3}); seqparm.setbacktime = str2num(input{4}); 
    else
       ok = 0;
    end
end
if ok
    input = inputdlg({'Significance level'; 'Number of local shuffles'; 'Number of global shuffles'; 'Gloabal shuffle mode (GlobalID/CircularSlide)'; 'Shuffle data save (raw/quant) rates (1/#)'},...
        'Significance parameters', 5, {'0.05'; '50'; '5'; 'GlobalID'; '5 5'}); %{'0.05'; '1000'; '100'; 'GlobalID'; '100 40'}); 
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
end
if ok
    cpathname = fullfile(cd, '*.tmp');
    [fname, pname] = uigetfile(cpathname, 'Select a template file to open:');
    if fname ~= 0
       tmpfile = fullfile(pname, fname);
       disp(['--------> template file: ', tmpfile]); seqparm.tmpfile = tmpfile;
       S = load(tmpfile, '-mat'); tmp = S.pinfo; S = [];
       if strcmp(tmp(1).tmp.tmptype, '2Dratemap')
          disp(['--------> number of templates found: ', num2str(numel(tmp))]);
       else
          ok = 0; disp('--------> templates not generated in 2D spaces; aborted.');
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
       seqdata.parentfile = currentfilename; seqdata.seqmethod = 'Bayesian2D';
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
    disp(['--------> final dir now: ', fdirnow]);
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
        [evseq, shufevseq, tmpnow(j), evName, evType, evT, tottime] = sequencedataandshuffle2D(pinfo, data, behav, bhdata, neuronid, tmpnow(j), seqparm); %evseq{nevent file}.field{nevents, 1}(1); 
        %%%%%%%%%%%%
        %t1 = cputime; disp(['-------------> decoding cpu time: ', num2str(t1 - t0)]);
        [seqnow, seqdatanow] = computeassignseqproperties2D(evseq, shufevseq, seqparm, pinfo, neuronid, tmpnow(j), evName, evType, evT, tottime); %shufevseq{nshuf, nevent file}.field.(morefield){nevents, 1}(1)
        t2 = cputime; disp(['----------------> 2d decoding/shuffling/quantification cpu time: ', num2str(t2 - t0)]);
        %t2 = cputime; disp(['-------------> quantification cpu time: ', num2str(t2 - t1)]);
        [seq, seqdata] = assignseqoutput(seq, seqdata, seqnow, seqdatanow, nt);
    end
   end
end

function [evseq, shufevseq, tmpnow, evName, evType, evT, tottime] = sequencedataandshuffle2D(pinfo, data, behav, bhdata, neuronid, tmpnow, seqparm)
%%%%%first get template and raw spike data: match tmp and spike cell names
cellfilename = tmpnow.tmp.tmpfilename; tmpevname = tmpnow.evtfile;
nlet = numel(cellfilename); isok = zeros(1, nlet); indok = zeros(1, nlet);
for (i = 1:nlet)
    ti = find( strcmp(pinfo.general.parmfile(neuronid), cellfilename{i}) );
    if (numel(ti) == 1)
        isok(i) = 1; indok(i) = ti;
    end
end
iii = find(isok == 1); nlet = numel(iii); cellfilename = cellfilename(iii); 
Xgrid{1} = tmpnow.tmp.xgrid; Xgrid{2} = tmpnow.tmp.ygrid; cellrate = tmpnow.tmp.ratemap(iii); %%%ratemap (1:ygrid, 1:xgrid) could be multidimensional: size = (m,n,...)
shufcellrate = shuffltemplate(cellrate, seqparm); %%%This are the shuffled templates: shufcellrate{seqparm.Nshuffle}
%%%%%%%%%%%get raw spikedata to decode
indok = indok(iii); spikedata = cell(1, nlet);
for (i = 1:nlet) spikedata{i} = data.spike.spiketime{neuronid(indok(i))}; end %%%%This is the matching rawdata to analyze
disp(['----------------> number of cells found: ', num2str(numel(cellrate))]);
%%%%%%%transfer basic session variables - no decoding done on overall sessions
sessName = pinfo.general.sessionname{neuronid(1)}; 
sessST = pinfo.general.sessionstartT{neuronid(1)}; sessET = pinfo.general.sessionendT{neuronid(1)}; nsess = numel(sessName);
%%%%%%%Decoding/sequences only done on selected events
           %%%%%%%%%%%% !!!!!!!!!!!!temporary adjustment if 'modify' option is on: if theta, change to [max max] time
           %[evName, evType, evT] = filterevents(pinfo, data, seqparm.evkeyword, seqparm.evkeytype, neuronid, 'modify', tmpevname);
[evName, evType, evT, tottime] = filterevents(pinfo, data, seqparm,  neuronid, 'none');
%%%%%%%%%%%%done temporary adjustment!!!!!!!!!
nev = numel(evName); 
evsessid = zeros(1, nev); evsessName = cell(1, nev); allevletter = cell(1, nev); allevtime = cell(1, nev);
evseq = cell(1, nev); shufevseq = cell(seqparm.Nshuffle, nev);
for (i = 1:nev) %%%for each event file
    %%%%work out which session the event is in spikedb
    [evsessName{i}, evsessid(i)] = identifysession(evT{i}, sessName, sessST, sessET);
    evT{i}.start = evT{i}.start - seqparm.setbacktime; evT{i}.ent = evT{i}.ent + seqparm.setbacktime; 
    %%%%find actual 2D (and linearized if possible) positions for each lap in a event
    finaldirnow = pinfo.general.finaldir{neuronid(1)}; nep = numel(evT{i}.start);
    posid = []; evid = []; realtime= cell(nep,1); realx = cell(nep,1); realy = cell(nep,1); real1dx = cell(nep,1);
    [joint, spaceunit] = locatejoint(behav, bhdata, tmpevname, cellfilename); %%% the event (could be different from template - so may not be accurate) linearized on the template trajectory
    if (~isempty(evsessName{i}))
        posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evsessName{i}) );
    end
    if numel(posid) == 1
        postimestamp = bhdata.pos.postimestamp{posid}; 
        Pmarker = behav.parm.sessPmarker{posid}; allposmarker = behav.general.posMarker{posid}; 
        ik = find(strcmp(allposmarker, Pmarker)); rtime = []; rx = []; ry=[]; r1dx=[];
        if (numel(ik) == 1) 
            posXX = bhdata.pos.XX{posid}{ik}*behav.parm.pixelXSize(posid); posYY = bhdata.pos.YY{posid}{ik}*behav.parm.pixelYSize(posid);
            [rtime, rx, ry, r1dx] = findrealposnow(joint, evT{i}, postimestamp*behav.parm.timeunit(posid), posXX, posYY, spaceunit);%%% 1D real locs: same events lineared on all possible linear trajecs 
        else
            disp(['--------------> Warning: position diode for computing actual postions not found: ', Pmarker]);
        end
    else
        isp(['--------------> Warning: position data cannot be identified in the behavioral database: ', evName{1}]);
    end
    %%%%work out the ep structure - this is for the old version that events can be treated with different options
    %%% This is quite annoying, but have to work around it. Now event.start(1:n) -> ep{n,1}.start(1); each ep{ii,jj} below contains only one event
    evTime = evT{i}; ep = []; 
    for (ti = 1:nep)
         ep{ti,1}.start = evTime.start(ti); ep{ti, 1}.ent = evTime.ent(ti); ep{ti,1}.marker = evTime.marker{ti};
         if ~isempty(rtime) realtime{ti,1} = rtime{ti}; realx{ti,1} = rx{ti}; realy{ti,1} = ry{ti}; real1dx{ti,1} = r1dx{ti}; end
    end
    datanow = cell(1, nlet); for (j = 1:nlet) [datanow{j}, ~] = SpikeEventFilter(spikedata{j}, evT{i}); end
    %%%%%%% Output are the raw decoded results
    evseq{i} = finddecoderesults(datanow, cellrate, Xgrid, seqparm, ep, realtime, realx, realy, real1dx, evName{i}, allevletter, allevtime, joint, spaceunit, 1, 1);
    %%%%%%% compute shuffles
    saveshufrawdata = zeros(1, seqparm.Nshuffle); saveshufquandata = zeros(1, seqparm.Nshuffle);
    if (seqparm.saveshufrawdatarate<=seqparm.Nshuffle) ii = 1:seqparm.saveshufrawdatarate:seqparm.Nshuffle; saveshufrawdata(ii) = ones(size(ii)); end
    if (seqparm.saveshufquandatarate<=seqparm.Nshuffle) ii = 1:seqparm.saveshufquandatarate:seqparm.Nshuffle; saveshufquandata(ii) = ones(size(ii)); end
    for (j = 1:seqparm.Nshuffle)
        shufevseq{j, i} = finddecoderesults(datanow, shufcellrate{j}, Xgrid, seqparm, ep, realtime, realx, realy, real1dx, evName{i}, allevletter, allevtime, joint, spaceunit, saveshufrawdata(j), saveshufquandata(j));
    end
end
function evseq = finddecoderesults(datanow, cellrate, Xgrid, seqparm, ep, realtime, realx, realy, real1dx, evName, allevletter, allevtime, joint, spaceunit, saverawflag, savequantflag)
[mm,nn] = size(ep); prob = cell(mm,nn); timepoint = cell(mm,nn);
seqstart = cell(mm,nn); seqend = cell(mm,nn); seqmarker = cell(mm,nn); 
matchscore = cell(mm,nn); posmatchprob = cell(mm,nn); negmatchprob = cell(mm,nn); seqtime = cell(mm,nn);
shufZscore = cell(mm,nn); shufZposP = cell(mm,nn); shufZnegP = cell(mm,nn); shufposP = cell(mm,nn); 
shufnegP = cell(mm,nn); decodeerr = cell(mm,nn); nactivecell = cell(mm,nn);  newquant = cell(mm,nn); quantdata = cell(mm,nn);
for (ii = 1:mm) %%%original event index
     for (jj = 1:nn)
         [prob{ii,jj}, timepoint{ii,jj}, seqstart{ii,jj}, seqend{ii,jj}, seqmarker{ii,jj}, nactivecell{ii,jj}, newquant{ii,jj}, quantdata{ii,jj},  ...
                matchscore{ii,jj}, posmatchprob{ii,jj}, negmatchprob{ii,jj}, ...
                 shufZscore{ii,jj}, shufZposP{ii,jj}, shufZnegP{ii,jj},...
                 shufposP{ii,jj}, shufnegP{ii,jj}, decodeerr{ii,jj}] = ...
                findallpositionprob(datanow, cellrate, Xgrid, seqparm, ep{ii,jj}, realtime{ii,jj}, realx{ii,jj}, realy{ii,jj}, real1dx{ii,jj}, joint, spaceunit);
     end
end
evseq.seqdata = []; if saverawflag evseq.seqdata = prob; evseq.timepoint = timepoint; evseq.Xgrid = Xgrid;  end
evseq.quantdata = []; if savequantflag evseq.quantdata = quantdata; end
%%%%%%%%%%%%%%% empty assignements
evseq.timepeakth = []; evseq.ratestd = []; evseq.spatratemean = []; evseq.spatpeakth = []; evseq.spatratestd = []; evseq.spatratemean = [];
%%%%%%%%%%%%%%%output
evseq.evname = evName; evseq.ep = ep; evseq.allletter = allevletter; evseq.alltime = allevtime;
evseq.seq = []; evseq.seqstart = seqstart; evseq.seqend = seqend; evseq.seqmarker = seqmarker; evseq.seqtime = seqtime;
evseq.newquant = newquant;
evseq.matchscore = matchscore; evseq.posmatchprob = posmatchprob; evseq.negmatchprob = negmatchprob;
evseq.shufZscore = shufZscore; evseq.shufZposP = shufZposP; evseq.shufZnegP = shufZnegP;
evseq.shufposP = shufposP; evseq.shufnegP = shufnegP; evseq.decodeerr = decodeerr; evseq.nactivecell = nactivecell;

% %%%%%%%%%%%%%%% old output
% evseq.evname = evName; evseq.ep = ep; evseq.allletter = allevletter; evseq.alltime = allevtime;
% evseq.seq = []; evseq.seqstart = seqstart; evseq.seqend = seqend; evseq.seqmarker = seqmarker; evseq.seqtime = seqtime;
% evseq.trajevname = jointevname; evseq.trajjoint = joint; 
% evseq.matchscore = matchscore; evseq.posmatchprob = posmatchprob; evseq.negmatchprob = negmatchprob;
% evseq.shufZscore = shufZscore; evseq.shufZposP = shufZposP; evseq.shufZnegP = shufZnegP;
% evseq.shufposP = shufposP; evseq.shufnegP = shufnegP; evseq.decodeerr = decodeerr; evseq.nactivecell = nactivecell;
% evseq.prop2D = prop2D;

function [seq, seqdata] = computeassignseqproperties2D(evseq, shufevseq, seqparm, pinfo, allneuronid, tmpnow, evName, evType, evTimes, tottime)
%%%%Properties to compute: 1-template propert; 2-decoded position properties; 3-global significance
%%%%%general variables
neuronid = allneuronid(1); %seq.general.seqtype = 'Bayesian2D'; %%%already in seq.parm.seqtype
seq.general.tmpID = strcat(tmpnow.parm.animaldate, '_', tmpnow.evtfile); 
seq.general.recarea = pinfo.general.recarea{neuronid}; seq.general.finaldir = pinfo.general.finaldir{neuronid};
seq.general.datedir = pinfo.general.datedir{neuronid}; seq.general.animalname = pinfo.general.animalname{neuronid};
seq.general.sessionname = pinfo.general.sessionname{neuronid}; seq.general.sessionstartT = pinfo.general.sessionstartT{neuronid};
seq.general.sessionendT = pinfo.general.sessionendT{neuronid}; seq.general.sessionlength = pinfo.general.sessionlength{neuronid};
seq.general.eventname = evName; seq.general.genotype = pinfo.general.genotype{neuronid}; %seq.general.lineartrajname = jointevname;
seq.general.sex = pinfo.general.sex{neuronid}; seq.general.age = pinfo.general.age{neuronid};
%%%%parameters
seq.parm = seqparm; seq.parm = rmfield(seq.parm, 'tmpfile'); 
seq.parm.sessType = pinfo.parm.sessType{neuronid}; seq.parm.eventType = evType;
%%%%%%%%%1.template properties - quality and similarities with other templates in the same finaldir
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
%%%%%%%%2. 2D decoded position properties: 
nev = numel(evName);
[P2D, candevind, targevind] = find2Dquantifications(evseq, tottime, seqparm, nev); %%%This is the final output variable: P2D.(field)(1, nev)
%%%%%%%%%%%%%%% 2.1D decoded position properties if available: use shufposP/shufnegP to define significant events
%%% 1D properties are computed against each linear trajectories (on the day: joint, jointevname)
%%%% nj = 1;          %%%now decide to quantify 1D only on the template traj
P1D = find1Dquantifications(evseq, tottime, candevind, targevind, P2D.totalN, P2D.candN, P2D.targN, seqparm, nev); %%%P1D.(field)(nev)
%%%%% 3. shuffle significance for both 2D and 1D variables
[sig2D, sig1D] = findshufflesignicance(P2D, P1D, shufevseq, tottime, seqparm, nev);
%%%%%%Output structure here
seq.seq = P1D; seq.shufsig = sig1D; %%%%this is the same as in 1D decoding
seq.seq2D = P2D; seq.shufsig2D = sig2D; 
%%%%%%Raw data output first
seqdata.events.eventimes = evTimes; seqdata.events.candevind = candevind; seqdata.events.targevind = targevind;
%seqdata.lineartraj.joint = joint; seqdata.lineartraj.jointevname = jointevname; %%%not sure what to do with joint/jointrevname
seqdata.data.evseq = evseq; seqdata.data.shufevseq = shufevseq;

function [sig2D, sig1D] = findshufflesignicance(P2D, P1D, shufevseq, tottime, seqparm, nev)
%%%%% evseq{nev}.(field){mm,nn}(1); shufevseq{nshuf, nev}.(field){mm,nn}(1)
sig2D = []; sig1D = []; nshuf = seqparm.Nshuffle; S2D = cell(1, nshuf); S1D = cell(1, nshuf); evseqnow = cell(1, nev);
for (i = 1:nshuf) 
    [S2D{i}, candevindnow, targevindnow] = find2Dquantifications(shufevseq(i,:), tottime, seqparm, nev);
    S1D{i} = find1Dquantifications(shufevseq(i,:), tottime, candevindnow, targevindnow, S2D{i}.totalN, S2D{i}.candN, S2D{i}.targN, seqparm, nev);
end
%%%%2D variable Z scores and P values relative to shuffles
S2D = reformatcleanup(S2D); %%%%re-format S2D{nshuf}.(field)(1:nev) to S2D.(field){nev}(nshuf)
fname = fieldnames(P2D);
for (j = 1:numel(fname))% For each variable/field, compare P2D.(field)(1:nev) with S2D.(field)(nshuf:nev)
    f1name = strcat(fname{j}, 'Z'); f2name = strcat(fname{j}, 'P');
    sig2D.(f1name) = NaN*ones(1, nev); sig2D.(f2name) = NaN*ones(1, nev);
    for (i=1:nev)
         vnow = S2D.(fname{j}){i}; vnow = vnow(~isnan(vnow)); nn = numel(vnow); 
         %disp(fname{j}); disp(P2D.(fname{j})); disp('****');
         pnow = P2D.(fname{j})(i);
         if nn >=3  
            shufmean = mean(vnow); shufstd = std(vnow); Z = (pnow-shufmean)/shufstd; 
            if Z>=0 
               P = numel(find(vnow>=pnow))/nn;
            else
               P = numel(find(vnow<=pnow))/nn;
            end
         end
         sig2D.(f1name)(i) = Z; sig2D.(f2name)(i) = P;
    end
end
%%%%1D variable Z scores and P values relative to shuffles
S1D = reformatcleanup(S1D); %%%%re-format S2D{nshuf}.(field)(1:nev) to S2D.(field){nev}(nshuf)
fname = fieldnames(P1D);
for (j = 1:numel(fname))% For each variable/field, compare P2D.(field)(1:nev) with S2D.(field)(nshuf:nev)
    f1name = strcat(fname{j}, 'Z'); f2name = strcat(fname{j}, 'P');
    sig1D.(f1name) = NaN*ones(1, nev); sig1D.(f2name) = NaN*ones(1, nev);
    for (i=1:nev)
         vnow = S1D.(fname{j}){i}; vnow = vnow(~isnan(vnow)); nn = numel(vnow); 
         pnow = P1D.(fname{j})(i);
         if nn >=3  
            shufmean = mean(vnow); shufstd = std(vnow); Z = (pnow-shufmean)/shufstd; 
            if Z>=0 
               P = numel(find(vnow>=pnow))/nn;
            else
               P = numel(find(vnow<=pnow))/nn;
            end
         end
         sig1D.(f1name)(i) = Z; sig1D.(f2name)(i) = P;
    end
end
% %%%%%%%%%%%%%this needs to dedone
% S1D = reformatcleanup1D(S1D); %%%%re-format S1D{nshuf}.(field){1:nev}(1:nj) to S1D/S2D.(field){nev}{nj}(nshuf)
% fname = fieldnames(P1D);
% for (j = 1:numel(fname))% For each variable/field, compare P1D.(field){nev}{nj} with S2D.(field){nev}{nj}(nshuf)
%     %%%%all field{nev}{nj} contains vector numbers ( no cells unlide P2D)
%     f1name = strcat(fname{j}, 'Z'); f2name = strcat(fname{j}, 'P');
%     sig1D.(f1name) = cell(1, nev); sig1D.(f2name) = cell(1, nev);
%     for (i=1:nev)
%         P1D.(f1name){i} = NaN*ones(1, nj); P1D.(f2name){i} = NaN*ones(1, nj);
%         for (k= 1:nj)
%             vnow = S1D.(fname{j}){i}{k}; vnow = vnow(~isnan(vnow)); nn = numel(vnow); 
%             pnow = P2D.(fname{j}){i}(k);
%             if nn >=3  
%                shufmean = mean(vnow); shufstd = std(vnow); Z = (pnow-shufmean)/shufstd; 
%                if Z>=0 
%                    P = numel(find(vnow>=pnow))/nn;
%                else
%                    P = numel(find(vnow<=pnow))/nn;
%                end
%             end
%             sig1D.(f1name){i}(k) = Z; sig1D.(f2name){i}(k) = P;
%         end
%     end
% end
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
function P = reformatcleanup1D(S) %%%%re-format S{nshuf}.(field){1:nev}(1:nj) to S.(field){nev}{nj}(nshuf)
nshuf = numel(S); P = [];
if nshuf>0 
    fname = fieldnames(S{1});
    if ~isempty(fname)
        nev = numel(S{1}.(fname{1}));
        for (i=1:numel(fname)) %%%all fields in S1D are cells with numbers(1:nj)
              P.(fname{i}) = cell(1, nev);
        end
        if nev >0
            nj = numel(S{1}.(fname{1}){1});
            for (i=1:numel(fname))
                for (k = 1:nev) 
                    P.(fname{i}){k} = cell(1,nj);
                    Tnow = zeros(nshuf,nj);
                    for (j = 1:nshuf)
                        Tnow(j,:) = S{j}.(fname{i}){k};
                    end
                    for (m=1:nj)
                        P.(fname{i}){k}{m} = Tnow(:,m);
                    end
                end
            end
        end
    end
end
    
function P = find1Dquantifications(evseq, tottime, candevind, targevind,totalN, candN, targN, seqparm, nj)
%%%%%%%%%canevind not used here: assuming all computed events are candidate events
%%%%%candidate/target events: matchscore(time-pos correlation), matchZ (corrrelative to shuffle), 
%%%%%                  significant positive (forward)/negative(reverse) matche N and ratio
P.meanScoreCand = NaN*ones(1, nj); P.medScoreCand = NaN*ones(1, nj); P.meanScoreZCand = NaN*ones(1, nj); %%% mean/med score values, Z means
P.med1dErrCand = NaN*ones(1, nj); P.mean1dErrCand = NaN*ones(1, nj); %%%nj here is the same as nev below
P.posNCand = NaN*ones(1, nj); P.posTimeRateCand = NaN*ones(1, nj); P.posEvtRateCand = NaN*ones(1, nj); P.posRatioCand = NaN*ones(1, nj);
P.negNCand = NaN*ones(1, nj); P.negTimeRateCand = NaN*ones(1, nj); P.negEvtRateCand = NaN*ones(1, nj); P.negRatioCand = NaN*ones(1, nj);
P.matchNCand = NaN*ones(1, nj); P.matchTimeRateCand = NaN*ones(1, nj); P.matchEvtRateCand = NaN*ones(1, nj); P.matchRatioCand = NaN*ones(1, nj);
P.meanScoreTarg = NaN*ones(1, nj); P.medScoreTarg = NaN*ones(1, nj); P.meanScoreZTarg = NaN*ones(1, nj); P.med1dErrTarg = NaN*ones(1, nj); P.mean1dErrTarg = NaN*ones(1, nj);
P.posNTarg = NaN*ones(1, nj); P.posTimeRateTarg = NaN*ones(1, nj); P.posEvtRateTarg = NaN*ones(1, nj); P.posRatioTarg = NaN*ones(1, nj);
P.negNTarg = NaN*ones(1, nj); P.negTimeRateTarg = NaN*ones(1, nj); P.negEvtRateTarg = NaN*ones(1, nj); P.negRatioTarg = NaN*ones(1, nj);
P.matchNTarg = NaN*ones(1, nj); P.matchTimeRateTarg = NaN*ones(1, nj); P.matchEvtRateTarg = NaN*ones(1, nj); P.matchRatioTarg = NaN*ones(1, nj);

for (j = 1:nj)    
    %%% format of decoding result: evseq{i}.shufnegP{mm, nn}(1): mm = nep, nn = 1 scaling factor
    %%%%%%% candidate events
        shufposall = repackoutsinglefield(evseq{j}.shufposP); shufnegall = repackoutsinglefield(evseq{j}.shufnegP); %%prop{nep, 1}(1) to pp(1:nep)
        nactall = repackoutsinglefield(evseq{j}.nactivecell); 
        scoreall = repackoutsinglefield(evseq{j}.matchscore); scoreZall = repackoutsinglefield(evseq{j}.shufZscore);
        errall = repackoutsinglefield(evseq{j}.decodeerr); 
        cind = candevind{j}; shufpos = shufposall(cind); shufneg = shufnegall(cind); nact = nactall(cind); score = scoreall(cind); scoreZ = scoreZall(cind); err = errall(cind);
    P.posNCand(j) = numel(find( (shufpos<seqparm.significancelevel) & (nact>=seqparm.minactivecell)) );
    P.negNCand(j) = numel(find( (shufneg<seqparm.significancelevel) & (nact>=seqparm.minactivecell)) );
    P.matchNCand(j) = P.posNCand(j) + P.negNCand(j);
    P.meanScoreCand(j) = mean(score(~isnan(score))); P.medScoreCand(j) = median(score(~isnan(score))); P.meanScoreZCand(j) = mean(scoreZ(~isnan(scoreZ)));
    P.med1dErrCand(j) = median(err(~isnan(err))); P.mean1dErrCand(j) = mean(err(~isnan(err)));
    P.posEvtRateCand(j) = P.posNCand(j)/totalN(j); P.posTimeRateCand(j) = P.posNCand(j)/tottime(j); P.posRatioCand(j) = P.posNCand(j)/candN(j);
    P.negEvtRateCand(j) = P.negNCand(j)/totalN(j); P.negTimeRateCand(j) = P.negNCand(j)/tottime(j); P.negRatioCand(j) = P.negNCand(j)/candN(j);
    P.matchEvtRateCand(j) = P.matchNCand(j)/totalN(j); P.matchTimeRateCand(j) = P.matchNCand(j)/tottime(j); P.matchRatioCand(j) = P.matchNCand(j)/candN(j);
    %%%%%%% target events    
        ind = targevind{j}; shufpos = shufposall(ind); shufneg = shufnegall(ind); nact = nactall(ind); score = scoreall(ind); scoreZ = scoreZall(ind); err = errall(ind);
    P.posNTarg(j) = numel(find( (shufpos<seqparm.significancelevel) & (nact>=seqparm.minactivecell)) );
    P.negNTarg(j) = numel(find( (shufneg<seqparm.significancelevel) & (nact>=seqparm.minactivecell)) ); 
    P.matchNTarg(j) = P.posNTarg(j) + P.negNTarg(j);
    P.meanScoreTarg(j) = mean(score(~isnan(score))); P.medScoreTarg(j) = median(score(~isnan(score))); P.meanScoreZTarg(j) = mean(scoreZ(~isnan(scoreZ)));
    P.med1dErrTarg(j) = median(err(~isnan(err))); P.mean1dErrTarg(j) = mean(err(~isnan(err)));
    P.posEvtRateTarg(j) = P.posNTarg(j)/totalN(j); P.posTimeRateTarg(j) = P.posNTarg(j)/tottime(j); P.posRatioTarg(j) = P.posNTarg(j)/targN(j);
    P.negEvtRateTarg(j) = P.negNTarg(j)/totalN(j); P.negTimeRateTarg(j) = P.negNTarg(j)/tottime(j); P.negRatioTarg(j) = P.negNTarg(j)/targN(j);
    P.matchEvtRateTarg(j) = P.matchNTarg(j)/totalN(j); P.matchTimeRateTarg(j) = P.matchNTarg(j)/tottime(j); P.matchRatioTarg(j) = P.matchNTarg(j)/targN(j);    
end

function [P, candevind, targevind] = find2Dquantifications(evseq, tottime, seqparm, nev)
%%%%% All cand events stat: Nstep, Ngstep, PProb, EndDis, Length, MedianStep, MaxStep, MinStep, MeanStep, StepVar
P.tottime = tottime; P.totalN = NaN*ones(1, nev); P.evTimeRate = NaN*ones(1, nev);%%%evseq{nev}.(field){mm,nn}(1)
P.candN = NaN*ones(1, nev); P.candTimeRate = NaN*ones(1, nev); P.candEvtRate = NaN*ones(1, nev);
P.targN = NaN*ones(1, nev); P.targTimeRate = NaN*ones(1, nev); P.targEvtRate = NaN*ones(1, nev); P.targRatio = NaN*ones(1, nev); 
P.NstepCand = NaN*ones(1, nev); P.NDstepCand = NaN*ones(1, nev); P.NVstepCand = NaN*ones(1, nev); P.NCstepCand = NaN*ones(1, nev);
P.medPProbCand = NaN*ones(1, nev); P.meanPProbCand = NaN*ones(1, nev);
P.medErrCand = NaN*ones(1, nev); P.meanErrCand = NaN*ones(1, nev); 
P.medLengthCand = NaN*ones(1, nev); P.meanLengthCand = NaN*ones(1, nev); P.medCLengthCand = NaN*ones(1, nev); P.meanCLengthCand = NaN*ones(1, nev); 
P.medEndDisCand = NaN*ones(1, nev); P.meanEndDisCand = NaN*ones(1, nev); 
P.meanStepCand = NaN*ones(1, nev); P.varStepCand = NaN*ones(1, nev); P.medStepCand = NaN*ones(1, nev);
P.maxStepCand = NaN*ones(1, nev); P.minStepCand = NaN*ones(1, nev); 
%%%%% All target events stat: Nstep, Ngstep, PProb, EndDis, Length, MedianStep, MaxStep, MinStep, MeanStep, StepVar
P.NstepTarg = NaN*ones(1, nev); P.NDstepTarg = NaN*ones(1, nev); P.NVstepTarg = NaN*ones(1, nev); P.NCstepTarg = NaN*ones(1, nev);
P.medPProbTarg = NaN*ones(1, nev); P.meanPProbTarg = NaN*ones(1, nev); 
P.medErrTarg = NaN*ones(1, nev); P.meanErrTarg = NaN*ones(1, nev); 
P.medLengthTarg = NaN*ones(1, nev); P.meanLengthTarg = NaN*ones(1, nev); P.medCLengthTarg = NaN*ones(1, nev); P.meanCLengthTarg = NaN*ones(1, nev);
P.medEndDisTarg = NaN*ones(1, nev); P.meanEndDisTarg = NaN*ones(1, nev); 
P.meanStepTarg = NaN*ones(1, nev); P.varStepTarg = NaN*ones(1, nev); P.medStepTarg = NaN*ones(1, nev);
P.maxStepTarg = NaN*ones(1, nev); P.minStepTarg = NaN*ones(1, nev);
targevind = cell(1, nev); candevind = cell(1, nev);
for (i = 1:nev) %for each event file
    prop2D = repackout(evseq{i}.newquant);  %%evseq{i}.newquant.(field){nep, 1}(1) to prop2D.(field)(1:nep) 
    ncell = repackoutsinglefield(evseq{i}.nactivecell); %%evseq{i}.nactivecell{nep, 1}(1) to ncell(1:nep) 
    P.totalN(i) = numel(prop2D.Nstep); 
    kk = find(ncell>=seqparm.minactivecell); P.candN(i) = numel(kk); candevind{i} = kk; %%%candidate events
    kk = find(prop2D.IsTarget); P.targN(i) = numel(kk); targevind{i} = kk; %%%target events
end
%%%%%% basic event rates
P.evTimeRate = P.totalN./tottime; P.candTimeRate = P.candN./tottime; P.candEvtRate = P.candN./P.totalN; 
P.targTimeRate = P.targN./tottime; P.targEvtRate = P.targN./P.totalN; P.targRatio = P.targN./P.candN;
for (i=1:nev) %%%%decoding properties
    %%%candidate events
    kk = candevind{i}; prop2D = repackout(evseq{i}.newquant); %%evseq{i}.newquant.(field){nep, 1}(1)
    NN = prop2D.Nstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.NstepCand(i) = median(NN(ii)); end
    NN = prop2D.NDstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.NDstepCand(i) = median(NN(ii)); end %%%number of Decoded steps
    NN = prop2D.NVstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.NVstepCand(i) = median(NN(ii)); end %%%number of Valid steps
    NN = prop2D.NCstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.NCstepCand(i) = median(NN(ii)); end %%%number of Continuous target steps
    NN = prop2D.MedianPProb(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medPProbCand(i) = median(NN(ii)); end
    NN = prop2D.MeanPProb(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.meanPProbCand(i) = mean(NN(ii)); end
    NN = prop2D.MedDecode2dErr(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medErrCand(i) = median(NN(ii)); end
    NN = prop2D.MeanDecode2dErr(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.meanErrCand(i) = mean(NN(ii)); end
    NN = prop2D.Length(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medLengthCand(i) = median(NN(ii)); P.meanLengthCand(i) = mean(NN(ii)); end
    NN = prop2D.CLength(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medCLengthCand(i) = median(NN(ii)); P.meanCLengthCand(i) = mean(NN(ii)); end
    NN = prop2D.EndDis(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medEndDisCand(i) = median(NN(ii)); P.meanEndDisCand(i) = mean(NN(ii)); end
    NN = prop2D.MeanStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.meanStepCand(i) = mean(NN(ii)); end
    NN = prop2D.StepVar(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.varStepCand(i) = mean(NN(ii)); end
    NN = prop2D.MedianStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medStepCand(i) = median(NN(ii)); end
    NN = prop2D.MaxStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.maxStepCand(i) = median(NN(ii)); end
    NN = prop2D.MinStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.minStepCand(i) = median(NN(ii)); end
    %%%target events
    kk = targevind{i};
    if ~isempty(kk)
       NN = prop2D.Nstep(kk); P.NstepTarg(i) = median(NN); NN = prop2D.NDstep(kk); P.NDstepTarg(i) = median(NN);
       NN = prop2D.NVstep(kk); P.NVstepTarg(i) = median(NN); NN = prop2D.NCstep(kk); P.NCstepTarg(i) = median(NN);
       NN = prop2D.MedianPProb(kk); P.medPProbTarg(i) = median(NN); NN = prop2D.MeanPProb(kk); P.meanPProbTarg(i) = mean(NN);
       NN = prop2D.MedDecode2dErr(kk); P.medErrTarg(i) = median(NN); NN = prop2D.MeanDecode2dErr(kk); P.meanErrTarg(i) = mean(NN);
       NN = prop2D.Length(kk); P.medLengthTarg(i) = median(NN); P.meanLengthTarg(i) = mean(NN);
       NN = prop2D.CLength(kk); P.medCLengthTarg(i) = median(NN); P.meanCLengthTarg(i) = mean(NN);   
       NN = prop2D.EndDis(kk); P.medEndDisTarg(i) = median(NN); P.meanEndDisTarg(i) = mean(NN);
       NN = prop2D.MeanStep(kk); P.meanStepTarg(i) = mean(NN); NN = prop2D.StepVar(kk); P.varStepTarg(i) = mean(NN);
       NN = prop2D.MedianStep(kk); P.medStepTarg(i) = median(NN);
       NN = prop2D.MaxStep(kk); P.maxStepTarg(i) = median(NN); NN = prop2D.MinStep(kk); P.minStepTarg(i) = median(NN); 
    end
end

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

function [prob, timepoint, seqstart, seqend, seqmarker, nactivecell, prop2D, prop2Ddata, ... %%% correspodong to newquant, quantdata: decoding a single event file: [start end] times
    matchscore, posmatchprob, negmatchprob, shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decode1derr] = ...
               findallpositionprob(spiketime, cellrate, Xgrid, seqparm, ep, realtime, realx, realy, real1dx, joint, spaceunit)
decparm.winsize = seqparm.windowsize; decparm.shifttime = seqparm.shifttime;
[timepoint, prob] = DE_Seq_DecodePositionProbability(spiketime, Xgrid, cellrate, ep, decparm); %since ep only contains 1 event, timepoint{1}
seqstart = ep.start; seqend = ep.ent; seqmarker = ep.marker;
nactivecell = findactivecell(spiketime, ep);
%%%% Quantify 2D decoded properties
[prop2D, prop2Ddata] = find2Ddecodingquantifications(timepoint, prob, Xgrid, seqparm, realtime, realx, realy);
%%%%%%%%1D is done separately for each available linear traj - NOT DONE!!!!
%%%%%%%%this - only for the template traj - since too many layers of hiararchy - no way to display in DataManager
% nj = numel(joint); matchscore = cell(1, nj); posmatchprob = cell(1, nj); negmatchprob = cell(1, nj); shufZscore = cell(1, nj); shufZposP = cell(1, nj);
% shufZnegP = cell(1, nj); shufposP = cell(1, nj); shufnegP = cell(1,nj); decode1derr = cell(1,nj);
% for (i = 1:nj)
%      [matchscore{i}, posmatchprob{i}, negmatchprob{i}, ...
%                  shufZscore{i}, shufZposP{i}, shufZnegP{i}, shufposP{i}, shufnegP{i}, decode1derr{i}] = finddecodingquantifications(timepoint, prop2D.peakx, prop2D.peaky, seqparm, realtime, real1dx{i}, joint{i});
% end
[matchscore, posmatchprob, negmatchprob, ...
                  shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decode1derr] = finddecodingquantifications(timepoint, prop2Ddata.peakx, prop2Ddata.peaky, seqparm, realtime, real1dx, joint, spaceunit);
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
                allevnames = pinfo.general.eventname; allevtypes = pinfo.parm.eventType; allevtimes = data.event.eventtimes{cellid};
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

function [evName, evType, evT, tottime] = filterevents(pinfo, data, seqparm, neuronid, ifmod)
evkeyword = seqparm.evkeyword; evkeytype = seqparm.evkeytype; evkeynoword = seqparm.evkeynoword; evkeynotype = seqparm.evkeynotype;
evName = pinfo.general.eventname{neuronid(1)};
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

function [prop2D, prop2Ddata] = find2Ddecodingquantifications(timepoint, prob, Xgrid, seqparm, alltime, allx, ally)
%%%%%timepoint{nep}(ntimepoint); prob{nep}{ntimepoint}(nY, nX); Xgrid{1} = X, Xgrid{2} = Y. 
nep = numel(timepoint); windowsize = seqparm.windowsize;
%%%%%determine decoded 2D positions and step distances
peakx = cell(1, nep); peaky = cell(1, nep); peakprob = cell(1, nep); dis = cell(1, nep);
for (i = 1:nep) %%%due to event re-formatiting, here nep = 1
    timenow = timepoint{i}; ntime = numel(timenow);
    peakx{i} = NaN*ones(size(timenow)); peaky{i} = NaN*ones(size(timenow)); peakprob{i} = NaN*ones(size(timenow)); dis{i} = NaN*ones(size(timenow));
    %nonselect = NaN*ones(size(timenow)); %%%careful about selection sequence: calculate first then select out -- This is OK
    for (j = 1:ntime)
        ttt = find(~isnan(prob{i}{j})); 
        if (numel(ttt)>0)
           [xind, yind, pprob] = findpeakxyind(prob{i}{j}); 
           if pprob > seqparm.probthres/numel(ttt) %%% if any threshold to select peak locations: higher than 0 is imposed
                peakx{i}(j) = Xgrid{1}(xind); peaky{i}(j) = Xgrid{2}(yind); peakprob{i}(j) = pprob;
           end
        end
    end
    for (j = 2:ntime)
        if isempty(find(isnan([peakx{i}(j-1) peaky{i}(j-1) peakx{i}(j) peaky{i}(j)])))
            dis{i}(j) = sqrt((peakx{i}(j)-peakx{i}(j-1))^2 + (peaky{i}(j)-peaky{i}(j-1))^2);
        end
    end
  %%%%% The following is not implemented, beacause it is hard to decide which of the two neighbors to exclude if distance is abnornmally high 
%     for (j = 2:ntime)
%         iii = find(dis{i}>seqparm.stepdisthres); %%% if any max jump distance to select peak locations: exclude if max jump is too high
%         if ~isempty(iii)
%             dis{i}(iii) = NaN*ones(size(iii)); peakx{i}(iii) = NaN*ones(size(iii)); peaky{i}(iii) = NaN*ones(size(iii)); peakprob{i}(iii) = NaN*ones(size(iii));
%         end
%     end
end
prop2Ddata.peakx = peakx; prop2Ddata.peaky = peaky; prop2Ddata.peakprob = peakprob; prop2Ddata.stepdis = dis;
%%%%% 2D quantification: meanpprob, median distance, overall distance, distance variance. distance = speed if timebin is constant 
prop2D.Nstep = NaN*ones(1, nep); prop2D.NDstep = NaN*ones(1, nep); prop2D.NVstep = NaN*ones(1, nep); prop2D.NCstep = NaN*ones(1, nep); 
prop2D.MedianPProb = NaN*ones(1, nep); prop2D.MeanPProb = NaN*ones(1, nep);
prop2D.EndDis = NaN*ones(1, nep); prop2D.Length = NaN*ones(1, nep); prop2D.CLength = NaN*ones(1, nep);
prop2D.MedianStep = NaN*ones(1, nep); prop2D.MaxStep = NaN*ones(1, nep); prop2D.MinStep = NaN*ones(1, nep); 
prop2D.MeanStep = NaN*ones(1, nep); prop2D.StepVar = NaN*ones(1, nep); 
for (i = 1:nep)
    prop2D.Nstep(i) = numel(timepoint{i})-1;
    iii = find(~isnan(dis{i})); prop2D.NVstep(i) = numel(iii); 
    if (~isempty(iii))
        disnow = dis{i}(iii); prop2D.MedianStep(i) = median(disnow); prop2D.MinStep(i) = min(disnow); prop2D.MaxStep(i) = max(disnow);
        prop2D.MeanStep(i) = mean(disnow); prop2D.StepVar(i) = var(disnow);
        prop2D.EndDis(i) = disnow(numel(disnow)) - disnow(1); prop2D.Length(i) = sum(disnow);
        [prop2D.NCstep(i), prop2D.CLength(i)] = findNcontinuoussteps(dis{i}, seqparm.minstepdis, seqparm.maxstepdis);
    end
    iii = find(~isnan(peakprob{i})); prop2D.NDstep(i) = numel(iii); 
    if (~isempty(iii))
        pnow = peakprob{i}(iii); prop2D.MedianPProb(i) = median(pnow); prop2D.MeanPProb(i) = mean(pnow);
    end
end
%%%%% 2D event quantifications: 2D (replay/trajectory) event detection numbers 
%%%%%    thresholds: seqparm.minmedpprob, seqparm.minNstep, seqparm.maxNstep, seqparm.minlength, seqparm.maxlength, 
prop2D.IsTarget = zeros(1, nep); %prop2D.NTarget = 0;
ind1 = find(prop2D.MedianPProb >= seqparm.minmedpprob); 
ind2 = find( (prop2D.NCstep>=seqparm.minNstep) & (prop2D.NCstep<=seqparm.maxNstep) );
ind3 = find( (prop2D.CLength>=seqparm.minlength) & (prop2D.CLength<=seqparm.maxlength) );
iii = intersect(intersect(ind1, ind2), ind3);
if ~isempty(iii)
    prop2D.IsTarget(iii) = ones(1, numel(iii)); %prop2D.NTarget = numel(iii);
end
%%%%% 2D quantification: decoding error if available
MedDecode2dErr = NaN*ones(1, nep); MeanDecode2dErr = NaN*ones(1, nep);
for (i = 1:nep)
    timenow = timepoint{i}; ntime = numel(timenow); realx = NaN*ones(size(timenow)); realy = NaN*ones(size(timenow));
    %%%%real 2D positions
    for (j = 1:ntime)
        jj = find( (alltime>=timenow(j)-windowsize/2) & (alltime<timenow(j)+windowsize/2) );
        if (numel(jj) >= 1) 
            xx = allx(jj); yy = ally(jj);
            xx = xx( (~isnan(xx)) & (~isnan(yy)) ); yy = yy( (~isnan(xx)) & (~isnan(yy)) );
            realx(j) = median(xx); realy(j) = median(yy);
        end
    end
    ind = find( (~isnan(realx)) & (~isnan(realy)) & (~isnan(peakx{i})) & (~isnan(peaky{i})) ); 
    if ~isempty(ind)
       aa = zeros(1, numel(ind));
       for kk = 1:numel(ind)
           aa(kk) = sqrt((realx(ind(kk))-peakx{i}(ind(kk)))^2 + (realy(ind(kk))-peaky{i}(ind(kk)))^2);
       end
       MedDecode2dErr(i) = median(aa); MeanDecode2dErr(i) = mean(aa); %%%takes the median/mean value as the error for the entire event
    end
end
prop2D.MedDecode2dErr = MedDecode2dErr; prop2D.MeanDecode2dErr = MeanDecode2dErr;
function [peakx, peaky, peakv] = findpeakxyind(ratemap)
[pprate, IY] = max(ratemap); [peakv, IX] = max(pprate); 
peakx = IX; peaky = IY(IX); %%%For ratemaps, x = column, y = row 
function [Cstep, CLeng] = findNcontinuoussteps(dnow, minstepdis, maxstepdis)
nc =0; maxc = 0; countD = 0; maxD = 0;
for (i=1:numel(dnow))
    if (dnow(i)<=maxstepdis) && (dnow(i)>=minstepdis)
       nc=nc+1; countD=countD+dnow(i); 
       if (nc>maxc) maxc = nc; maxD = countD; end
    else
       nc=0; countD=0; 
    end
end
Cstep = maxc; CLeng = maxD;

function [matchscore, posmatchprob, negmatchprob, ...
                 shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decode1derr] = finddecodingquantifications(timepoint, peakx, peaky, seqparm, alltime, all1dx, joint, spaceunit)
%%%%%timepoint{nep}(ntimepoint); prob{nep}{ntimepoint}(nY, nX); Xgrid{1} = X, Xgrid{2} = Y. 
%%%%%%%%%%% For 1D info if avaialble: determine the decoded position, then compute: decode error; pearson R, p; shuffle Z, p 
nshuffle = seqparm.Nwithineventshuffle; nep = numel(timepoint); windowsize = seqparm.windowsize;
matchscore = NaN*ones(1, nep); posmatchprob = NaN*ones(1, nep); negmatchprob = NaN*ones(1, nep); decode1derr = NaN*ones(1, nep);
shufZscore = NaN*ones(1, nep); shufZposP = NaN*ones(1, nep); shufZnegP = NaN*ones(1, nep); shufposP = NaN*ones(1, nep); shufnegP = NaN*ones(1, nep);
if (~isempty(joint)) && (~isempty(all1dx))
for (i = 1:nep) %%% due to event reformatting, here nep = 1
    %%%%%determine decoded positions
    timenow = timepoint{i}; ntime = numel(timenow); locnow = NaN*ones(size(timenow));
    for (j = 1:ntime) 
        if (~isnan(peakx{i}(j))) && (~isnan(peaky{i}(j)))
            [locnow(j), ~, ~] = LinearizePathXY_new(joint, peakx{i}(j), peaky{i}(j), 1, 0, spaceunit); %%% set to "1" meaning no continuity check, which only works for continuous moving; 
        end
    end
    iii = find(~isnan(locnow));
    if ~isempty(iii)
       timenow = timenow(iii); locnow = locnow(iii); ntime = numel(timenow);
       [timenow, iii] = sort(timenow); locnow = locnow(iii);
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
                  aaa = all1dx(jj); aaa = aaa(~isnan(aaa));
                  realloc(ii) = median(aaa); 
              end
          end
          dis = abs(locnow-realloc); 
          decode1derr(i) = median(dis(~isnan(dis)));
       end
       %%%shuffle Z (Zp), Sp
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
% 
% function [joint, jointevname] = locatejoint(behav, bhdata, tmpfilename)
% joint = []; nj = 0;
% %%%get the final directory
% tind = strfind(tmpfilename{1}, 'final'); tmpfdir = tmpfilename{1}(1:tind+4);
% sessid = find(strcmp(behav.general.finaldir, tmpfdir)); %%%behavioral data items are sessions
% for (i = 1:numel(sessid))
%     posltrfile = bhdata.pos.ltrfilename{sessid(i)}; jjdat = bhdata.pos.posltr{sessid(i)};
%     for (tt = 1:numel(posltrfile))
%         if ~isempty(jjdat{tt})
%            nj = nj +1; joint{nj} = jjdat{tt}; jointevname{nj} = posltrfile{tt};
%            joint{nj}(:,1) = jjdat{tt}(:,1)*behav.parm.pixelXSize(sessid); joint{nj}(:,2) = jjdat{tt}(:,2)*behav.parm.pixelYSize(sessid);
%         end
%     end
% end
function [joint, spaceunit] = locatejoint(behav, bhdata, evname, tmpfilename)
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
    evPosltr = behav.parm.eventPosltr{sessid}; spaceunit = (behav.parm.pixelXSize(sessid)+ behav.parm.pixelYSize(sessid))/2;
    for (tt = 1:numel(posltrfile))
        if (strcmp(posltrfile{tt}, evPosltr{evid}))
           joint = posltr{tt}; 
           joint(:,1) = joint(:,1)*behav.parm.pixelXSize(sessid); joint(:,2) = joint(:,2)*behav.parm.pixelYSize(sessid);
           break
        end
    end
end
function [realtime, realx, realy, real1dx] = findrealposnow(joint, evT, postimestamp, xpos, ypos, spaceunit)
nep = numel(evT.start); realx = cell(1, nep); realy = cell(1, nep); realtime = cell(1, nep); %disp('&&&&'); disp(joint); disp('%%%%');
real1dx = cell(1, nep); %%%compute real locations on only one (template) linear trajectories
for (tt = 1:nep)
    evposind = find( (postimestamp>=evT.start(tt)) & (postimestamp<=evT.ent(tt)) );
    realtime{tt} = postimestamp(evposind); realx{tt} =  xpos(evposind); realy{tt} =  ypos(evposind);
    real1dx{tt} = NaN*ones(size(realtime{tt}));
    if (~isempty(evposind)) && (~isempty(joint))
        [xout, ~, ~] = LinearizePathXY_new(joint, realx{tt}, realy{tt}, 0.95, 0, spaceunit);
        real1dx{tt} = xout; 
    end
end

function pp = repackout(prop)  %%prop{nep, 1}.(field)(1) to pp.(field)(1:nep) 
[mm,nn]=size(prop); pp = [];
if (mm*nn>=1)
   if (nn~=1)
       disp(['--------------------warning: unexpected dimensions encountered!']);
   end
   allf = fieldnames(prop{1,1});
   for (k = 1:numel(allf))
       ppnow = NaN*ones(1, mm*nn);
       for (i = 1:mm)
       for j = 1:nn
           ppnow((i-1)*nn+j) = prop{i,j}.(allf{k})(1);
       end
       end
       pp.(allf{k}) = ppnow;
   end
end
function pp = repackoutsinglefield(prop)  %%prop{nep, 1}(1) to pp(1:nep) 
[mm,nn]=size(prop); pp = NaN*ones(1, mm*nn);
if (nn~=1)
    disp(['--------------------warning: unexpected dimensions encountered!']);
end
for (i = 1:mm)
for (j = 1:nn)
    pp((i-1)*nn+j) = prop{i,j}(1);
end
end

function shufcellrate = shuffltemplate(cellrate, seqparm)
nshuffle = seqparm.Nshuffle; shufcellrate = cell(1,nshuffle);
ncell = numel(cellrate); 
if ncell>0
   if strcmp(seqparm.shufflemode, 'GlobalID')
      for (i = 1:nshuffle) shufcellrate{i} = cellrate(randperm(ncell)); end
   elseif strcmp(seqparm.shufflemode, 'CircularSlide') %%%%Cannot slide for 2D templates
      [nx, ny] = nsize(cellrate{1}); %%%this is for a 2D array
      for (i = 1:nshuffle)
          shufcellrate{i} = cellrate{i};
          for (j = 1:ncell)
              IX = mod((1:nx)+ ceil(rand*nx), nx) + 1; IY = mod((1:ny)+ ceil(rand*ny), ny) + 1;
              shufcellrate{i}{j} = cellrate{i}{j}(IX, IY);
          end
      end
   else
      disp('----------> can not generate shuffled data or shuffle mode not available; shuffle significance not computed');
   end
end



