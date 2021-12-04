function DM_Seq_EventItemized_Callback
%%%transform sequencing/decoding results in a .seqdb file from template items as event items 

hf = gcbf; pinfonow = getappdata(hf, 'pinfo'); datanow = getappdata(hf, 'data');
plotparm = getappdata(hf, 'plotparm');
%%%%data will not be changed;
%%%recomupte everything based on events, may need to redo the shuffling
disp('-----> Transforming a .seqdb file as event itemized......');
ok = 1;
%%%%general/parm variables: work out event items fdir by fdir
seqtype = 'Bayesian';
if isfield(pinfonow.parm, 'seqtype') seqtype = pinfonow.parm.seqtype{1}; end
pinfo = []; pinfo.general = [];
fdirall = pinfonow.general.finaldir;
fdir = unique(fdirall);
nitem = 0; tmpid = [];
allvar = fieldnames(pinfonow.general);
for (i = 1:numel(fdir))
    allevnow = cell(1,0);
    tmpnow = find( strcmp(fdirall, fdir{i}) );
    for (j = 1:numel(tmpnow))
        allevnow = [allevnow pinfonow.general.eventname{tmpnow(j)}];
    end
    evnow = unique(allevnow);
    for (k = 1:numel(evnow))
        nitem = nitem + 1; tmpid{nitem} = [];
        for (j = 1:numel(tmpnow))
            if ~isempty(find(strcmp(pinfonow.general.eventname{tmpnow(j)}, evnow{k})))
                tmpid{nitem} = [tmpid{nitem} tmpnow(j)];
            end
        end
        for (m = 1:numel(allvar))
            if (~strncmpi(allvar{m}, 'session', 4)) && (~strncmpi(allvar{m}, 'tmpID', 5))
                pinfo.general.(allvar{m}){nitem} = pinfonow.general.(allvar{m}){tmpid{nitem}(1)};    
            end
        end
        %%%need to reset eventnames, tmpID (below)
        eventname{nitem} = evnow{k};
        %pinfo.parm.eventType{nitem} = pinfonow.parm.eventType{tmpid{nitem}(1)};
    end
end
%%%reset tmpID
for (i = 1:nitem)
    pinfo.general.tmpname{i} = pinfonow.general.tmpID(tmpid{i});
end
%%%eventType
for (n = 1:nitem)
     evidnow = find(strcmp(pinfonow.general.eventname{tmpid{n}(1)}, eventname{n}));
     if numel(evidnow) == 1
        pinfo.parm.eventType{n} = pinfonow.parm.eventType{tmpid{n}(1)}{evidnow};
        evID = strcat(pinfonow.general.animalname{tmpid{n}(1)}, '_', pinfonow.general.datedir{tmpid{n}(1)}, '_',...
            eventname{n});
        pinfo.general.eventname{n} = evID;
     else
        disp(['---------> warning: event id not found: ', pinfo.general.eventname{n}, ' for template: ', pinfonow.general.tmpID{tmpid{n}(1)}]);
     end
end

% %%%%re-arrange the parm/tmp variables
% allvar = fieldnames(pinfonow.parm);
% for (m = 1:numel(allvar))
%     if ~strcmp(allvar{m}, 'eventType')
%         for (n = 1:nitem)
%              pinfo.parm.(allvar{m}){n} = pinfonow.parm.(allvar{m})(tmpid{n});
%         end
%     end
% end
% allvar = fieldnames(pinfonow.tmp);
% for (m = 1:numel(allvar))
%     if ~strcmp(allvar{m}, 'eventType')
%         for (n = 1:nitem)
%              pinfo.tmp.(allvar{m}){n} = pinfonow.tmp.(allvar{m})(tmpid{n});
%         end
%     end
% end
% 
% %%%%re-arrange the seq variables
% allvar = fieldnames(pinfonow.seq);
% for (m = 1:numel(allvar))
%     if (~strcmp(allvar{m}, 'evTotalN'))&(~strcmp(allvar{m}, 'evCandN'))&(~strcmp(allvar{m}, 'evCandLapRate'))...
%             &(~strcmp(allvar{m}, 'evCandTimeRate'))&(~strcmp(allvar{m}, 'evTotalTimeRate'))
%         newvarnow = strrep(allvar{m}, 'ev', 'tmp');
%         for (n = 1:nitem)
%             pinfo.seq.(allvar{m}){n} = NaN*ones(size(tmpid{n}));
%             for (k = 1:numel(tmpid{n}))
%                 %disp(pinfo.general.eventname{n})
%                 %disp(pinfonow.general.eventname{tmpid{n}(k)})
%                 evidnow = find(strcmp(pinfonow.general.eventname{tmpid{n}(k)}, pinfo.general.eventname{n}));
%                 if numel(evidnow) == 1
%                    pinfo.seq.(allvar{m}){n}(k) = pinfonow.seq.(allvar{m}){tmpid{n}(k)}(evidnow);
%                 else
%                    disp(['---------> warning: event id not found: ', pinfo.general.eventname{n}, ' for template: ', pinfonow.general.tmpID{tmpid{n}(k)}]);
%                 end
%             end
%         end
%     end
% end

%%%%compute event variables
pinfo.event = [];
pinfo.event.totalN = cell(1, nitem); pinfo.event.totalTimeRate = cell(1, nitem);
pinfo.event.candN = cell(1, nitem); pinfo.event.candTimeRate = cell(1, nitem); pinfo.event.candLapRate = cell(1, nitem);
pinfo.event.posN = cell(1, nitem); pinfo.event.posTimeRate = cell(1, nitem); pinfo.event.posLapRate = cell(1, nitem); 
pinfo.event.posRatio = cell(1, nitem); pinfo.event.posShuffleZ = cell(1, nitem); pinfo.event.posShuffleP = cell(1, nitem);
pinfo.event.negN = cell(1, nitem); pinfo.event.negTimeRate = cell(1, nitem); pinfo.event.negLapRate = cell(1, nitem); 
pinfo.event.negRatio = cell(1, nitem); pinfo.event.negShuffleZ = cell(1, nitem); pinfo.event.negShuffleP = cell(1, nitem);
pinfo.event.matchN = cell(1, nitem); pinfo.event.matchTimeRate = cell(1, nitem); pinfo.event.matchLapRate = cell(1, nitem);
pinfo.event.matchRatio = cell(1, nitem); pinfo.event.matchShuffleZ = cell(1, nitem); pinfo.event.matchShuffleP = cell(1, nitem);
datanow.event.shufmatchN = cell(1, nitem);
for (n = 1:nitem)
    %%%%%%update event numbers
    evidnow = find(strcmp(pinfonow.general.eventname{tmpid{n}(1)}, eventname{n}));
    if numel(evidnow) ~= 1
        disp(['---------> warning: event id not found: ', pinfo.general.eventname{n}, ' for template: ', pinfonow.general.tmpID{tmpid{n}(1)}]);
    else
        pinfo.event.totalN{n} = pinfonow.seq.evTotalN{tmpid{n}(1)}(evidnow); pinfo.event.totalTimeRate{n} = pinfonow.seq.evTotalTimeRate{tmpid{n}(1)}(evidnow);
        pinfo.event.candN{n} = pinfonow.seq.evCandN{tmpid{n}(1)}(evidnow); pinfo.event.candTimeRate{n} = pinfonow.seq.evCandTimeRate{tmpid{n}(1)}(evidnow);
        pinfo.event.candLapRate{n} = pinfonow.seq.evCandLapRate{tmpid{n}(1)}(evidnow);
        %%%%%%count sig /total/match/pos/neg matched numbers
        [posN, negN, matchN] = findmatchnumber(pinfonow, datanow, tmpid{n}, eventname{n}, pinfo.event.totalN{n});
        totaltime = pinfo.event.totalN{n}/pinfo.event.totalTimeRate{n};
        pinfo.event.posN{n} = posN; pinfo.event.posTimeRate{n} = posN/totaltime; 
        pinfo.event.posLapRate{n} = posN/pinfo.event.totalN{n}; pinfo.event.posRatio{n} = posN/pinfo.event.candN{n};
        pinfo.event.negN{n} = negN; pinfo.event.negTimeRate{n} = negN/totaltime; 
        pinfo.event.negLapRate{n} = negN/pinfo.event.totalN{n}; pinfo.event.negRatio{n} = negN/pinfo.event.candN{n};
        pinfo.event.matchN{n} = matchN; pinfo.event.matchTimeRate{n} = matchN/totaltime; 
        pinfo.event.matchLapRate{n} = matchN/pinfo.event.totalN{n}; pinfo.event.matchRatio{n} = matchN/pinfo.event.candN{n};
        %%%%%%compute match/pos/neg matched number Z, P
        evTimes = datanow.events.eventimes{tmpid{n}(1)}{evidnow};
        if strcmp(seqtype, 'Bayesian')
           [posZ, negZ, matchZ, posP, negP, matchP, shufmatchN] = findshufflesignificance_Bayesian(posN, negN, matchN, pinfonow, datanow, tmpid{n}, eventname{n}, evTimes);
        else
           [posZ, negZ, matchZ, posP, negP, matchP, shufmatchN] = findshufflesignificance_SeqMatching(posN, negN, matchN, pinfonow, datanow, tmpid{n}, eventname{n}, evTimes); 
        end
        pinfo.event.posShuffleZ{n} = posZ; pinfo.event.posShuffleP{n} = posP; 
        pinfo.event.negShuffleZ{n} = negZ; pinfo.event.negShuffleP{n} = negP;
        pinfo.event.matchShuffleZ{n} = matchZ; pinfo.event.matchShuffleP{n} = matchP;
        datanow.event.shufmatchN{n} = shufmatchN;
    end
end
%%%%%group data
pinfo.work = struct([]); data = datanow;
data.grouplist = []; 
data.grouplist.groupname{1} = 'List0'; 
data.grouplist.groupindex{1} = 1:numel(pinfo.general.eventname);
data.grouplist.grouptype{1} = 'Manual'; data.grouplist.groupcrit{1} = []; data.grouplist.groupparents{1} = []; 
fname = get(hf, 'Name'); ftt = strfind(fname, '__'); currentfilename = fname(ftt+2:numel(fname));
[pp, nn, ee] = fileparts(currentfilename);
data.parentfile = currentfilename; data.seqmethod = 'Bayesian';
filename = fullfile(pp, strcat(nn, '_eventitemized', '.seqdb'));
hmain = DataManager_DataManager_Callback;
setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); S = [];
DataManager_PlotSpikeDatabase(hmain, pinfo, data, 'Event', 'eventname', '.seqdb'); 
set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', filename));
save(filename, 'pinfo', 'data', '-mat');
disp('************************');

function [posZ, negZ, matchZ, posP, negP, matchP, shufmatchN] = findshufflesignificance_SeqMatching(posN, negN, matchN, pinfonow, datanow, tmpid, evName, evTimes)
posP = NaN; posZ = NaN; negP = NaN; negZ = NaN; matchP = NaN; matchZ = NaN;
evseq = datanow.data.evseq(tmpid); %%%all templates containing this event
allparmvar = fieldnames(pinfonow.parm);
for (i = 1:numel(allparmvar))
    seqparm.(allparmvar{i}) = pinfonow.parm.(allparmvar{i}){tmpid(1)};
end
% %%%%before the shuffle computation, need to get the source data: read at the spike trains and template data
% S = load(datanow.parentfile, '-mat'); %load the parent file
% pinfo = S.pinfo; data = S.data; S = [];%original spike data
% S = load(pinfonow.tmp.tmpfile{tmpid(1)}, '-mat'); tmp = S.pinfo; S = [];

%%%compute the shuffled numbers for each template
nshuffle = seqparm.Nshuffle; ntmp = numel(tmpid); ok = 1;
shufposN = NaN*ones(1, nshuffle); shufnegN = NaN*ones(1, nshuffle); shufmatchN = NaN*ones(1, nshuffle);
shufposP = cell(nshuffle, ntmp); shufnegP = cell(nshuffle, ntmp); nact = cell(nshuffle, ntmp); mm= 0; nn = 0;
for (tti = 1:ntmp) 
    evidnow = find(strcmp(pinfonow.general.eventname{tmpid(tti)}, evName));
    if numel(evidnow) ~= 1
        disp(['---------> warning: event id not found: ', evName, ' for template: ', pinfonow.general.tmpID{tmpid(tti)}]);
    else
        %%%%first: get template and raw data f
        ep = evseq{tti}{evidnow}.ep; [mm,nn] = size(ep);
        [MCroot, MCname, DAname, DEname, ANname, MAname] = CurrentVersion;
        rfile = fullfile(MCroot, 'DataExplorer', 'Sequence', 'rrankprob.mat'); 
        rrank = [];
        if (exist(rfile, 'file') == 2)
             S = load(rfile); rrank = S.rrank; S = [];
        else
             ok = 0; disp('-----> ranking probability file for selected rank mode not found; aborted');
        end
        if ok
             tmpseq = pinfonow.tmp.tmpseq{tmpid(tti)}; fileletter = cell2mat(pinfonow.tmp.fileletter{tmpid(tti)});
             tmprank = pinfonow.tmp.tmprank{tmpid(tti)};
             if strcmp(seqparm.shufflemode, 'GlobalID')
                 [seq, ~, posmatchprob, negmatchprob, ~, ~] = generateshuffle_gloabalID_seqmatch(tmpseq, fileletter, evseq{tti}{evidnow}, seqparm, rrank, tmprank, seqparm.eventoption);%, alltime);
             elseif strcmp(seqparm.shufflemode, 'LocalID')
                 [seq, ~, posmatchprob, negmatchprob, ~, ~] = generateshuffle_localID_seqmatch(tmpseq, fileletter, evseq, seqparm, rrank, tmprank, seqparm.eventoption);%, alltime);
             elseif strcmp(seqparm.shufflemode, 'CircleSlide')
                 sessST = pinfonow.general.sessionstartT{tmpid(tti)}; sessET = pinfonow.general.sessionendT{tmpid(tti)}; 
                 [seq, ~, posmatchprob, negmatchprob, ~, ~] = generateshuffle_circleslide_seqmatch(tmpseq, fileletter, evseq, seqparm, ...
                       sessST, sessET, rrank, tmprank, seqparm.eventoption, seqparm.eventtimingoption);
             else
                 ok = 0; disp('----------> can not generate shuffled data or shuffle mode not available; shuffle significance not computed');
             end
             if ok
                for (i = 1:nshuffle) %%%forget the ep{mm,nn} structure here: orginal event times are in evTimes
                     nactnow = cell(mm,nn);
                     for (ii = 1:mm)
                     for (jj = 1:nn)
                         for (kkti = 1:numel(seq{i}{mm,nn})) nactnow{mm,nn} = numel(seq{i}{mm,nn}{kkti});  end
                     end
                     end
                     shufposP{tti,i} = cell2mat(posmatchprob{i}); shufnegP{tti,i} = cell2mat(negmatchprob{i});
                     nact{tti,i} = cell2mat(nactnow);
                end
             end
        end
     end
end
if ok
for i = 1:nshuffle
    shufposPnow = NaN*ones(ntmp,mm*nn); shufnegPnow = NaN*ones(ntmp,mm*nn); nactnow = NaN*ones(ntmp,mm*nn);
    for (tti = 1:ntmp)
        if (numel(shufposP{tti,i}) == mm*nn)
            shufposPnow(tti,:) = shufposP{tti,i}; shufnegPnow(tti,:) = shufnegP{tti,i}; nactnow(tti,:) = nact{tti,i};
        else
            disp('----------> event lap number not match; results from a template may be ignored');
        end
    end
    [shufposN(i), shufnegN(i), shufmatchN(i)] = findmatchnumbernow(nactnow, shufposPnow, shufnegPnow, seqparm.significancelevel, seqparm.minactivecell);
    %disp([shufposN(i) shufnegN(i) shufmatchN(i)]);
end
posP = numel(find(shufposN>=posN))/nshuffle; negP = numel(find(shufnegN>=negN))/nshuffle; matchP = numel(find(shufmatchN>=matchN))/nshuffle; 
shufposN = shufposN(~isnan(shufposN)); shufnegN = shufnegN(~isnan(shufnegN)); shufmatchN = shufmatchN(~isnan(shufmatchN));
if ~isempty(shufposN) posZ = (posN-mean(shufposN))/std(shufposN); end
if ~isempty(shufnegN) negZ = (negN-mean(shufnegN))/std(shufnegN); end
if ~isempty(shufmatchN) matchZ = (matchN-mean(shufmatchN))/std(shufmatchN); end
end

function [posN, negN, matchN] = findmatchnumber(pinfonow, datanow, tmpid, eventname, evN)
ntmp = numel(tmpid); slevel = pinfonow.parm.significancelevel{1}; candN = pinfonow.parm.minactivecell{1};
shufposP = NaN*ones(ntmp, evN); shufnegP = NaN*ones(ntmp, evN); nactnum = NaN*ones(evN,1);
posPname = 'shufposP'; negPname = 'shufnegP';
for (k = 1:ntmp)
    if isfield(pinfonow.parm, 'seqtype')
        if strcmp(pinfonow.parm.seqtype{k}, 'SeqMatching')
            posPname = 'posmatchprob'; negPname = 'negmatchprob';
        end
    end
    evidnow = find(strcmp(pinfonow.general.eventname{tmpid(k)}, eventname));
    if numel(evidnow) == 1
        shufposP(k,:) = cell2mat(datanow.data.evseq{tmpid(k)}{evidnow}.(posPname));
        shufnegP(k,:) = cell2mat(datanow.data.evseq{tmpid(k)}{evidnow}.(negPname));
        nactnum = cell2mat(datanow.data.evseq{tmpid(k)}{evidnow}.nactivecell);
    else
        disp(['---------> warning: event id not found: ', eventname, ' for template: ', pinfonow.general.tmpID{tmpid(k)}]);
    end
end
[posN, negN, matchN] = findmatchnumbernow(nactnum, shufposP, shufnegP, slevel, candN);

function [posN, negN, matchN] = findmatchnumbernow(nactnum, shufposP, shufnegP, slevel, candN)
%%%%shufposP[ntmp, nlap]
[mm,nn] = size(nactnum);
if nn == 1
   nactnum = nactnum';
elseif (mm>1) && (nn>1)
   nactnum = min(nactnum);
end
shufposP = min(shufposP); shufnegP = min(shufnegP);
posN = numel(find( (shufposP<slevel) & (nactnum>candN) )); 
negN = numel(find( (shufnegP<slevel) & (nactnum>candN) )); 
matchN = numel(find( ((shufposP<slevel) | (shufnegP<slevel)) & (nactnum>candN) )); 

function [posZ, negZ, matchZ, posP, negP, matchP, shufmatchN] = findshufflesignificance_Bayesian(posN, negN, matchN, pinfonow, datanow, tmpid, evName, evTimes)
%%%%%3. shuffle significance
posP = NaN; posZ = NaN; negP = NaN; negZ = NaN; matchP = NaN; matchZ = NaN;
evseq = datanow.data.evseq(tmpid); %%%all templates containing this event
allparmvar = fieldnames(pinfonow.parm);
for (i = 1:numel(allparmvar))
    seqparm.(allparmvar{i}) = pinfonow.parm.(allparmvar{i}){tmpid(1)};
end
%%%%before the shuffle computation, need to get the source data: read at the spike trains and template data
S = load(datanow.parentfile, '-mat'); %load the parent file
pinfo = S.pinfo; data = S.data; S = [];%original spike data
S = load(pinfonow.tmp.tmpfile{tmpid(1)}, '-mat'); tmp = S.pinfo; S = [];
%%%compute the shuffled numbers for each template
nshuffle = seqparm.Nshuffle; 
ntmp = numel(tmpid);
shufposN = NaN*ones(1, nshuffle); shufnegN = NaN*ones(1, nshuffle); shufmatchN = NaN*ones(1, nshuffle);
shufposP = cell(nshuffle, ntmp); shufnegP = cell(nshuffle, ntmp); nact = cell(nshuffle, ntmp); mm= 0; nn = 0;
for (tti = 1:ntmp) 
    evidnow = find(strcmp(pinfonow.general.eventname{tmpid(tti)}, evName));
    if numel(evidnow) ~= 1
        disp(['---------> warning: event id not found: ', evName, ' for template: ', pinfonow.general.tmpID{tmpid(tti)}]);
    else
        %%%%first: get template and raw data for decoding
        Xgrid = evseq{tti}{evidnow}.Xgrid;
        cellfilename = pinfonow.tmp.filename{tmpid(tti)};
        ok = 1;
        [tmpnow, ok] = findtemplates(tmp, cellfilename);
      if ~ok
        disp(['---------> warning: template not found; template may be ignored: ', pinfonow.general.tmpID{tmpid(tti)}]);
      else
        nlet = numel(cellfilename); isok = zeros(1, nlet); indok = zeros(1, nlet);
        for (i = 1:nlet)
            ti = find( strcmp(pinfo.general.parmfile, cellfilename{i}) );
            if (numel(ti) == 1) isok(i) = 1; indok(i) =ti; end
        end
        iii = find(isok == 1); nlet = numel(iii); cellfilename = cellfilename(iii); cellrate = tmpnow.tmp.tmpspacerate(iii);
        indok = indok(iii); spikedata = cell(1, nlet);
        spikedata = data.spike.spiketime(indok); %%%%This is the rawdata to analyze
        dataok = cell(1, nlet);
        for (j = 1:nlet) [dataok{j}, ~] = SpikeEventFilter(spikedata{j}, evTimes); end
        ep = evseq{tti}{evidnow}.ep; [mm,nn] = size(ep);
        %%%%second: generate shuffles - all done on templates
        ncell = numel(cellrate); nbin = numel(Xgrid); shufcellrate = cell(1,nshuffle);
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
           ok = 0; disp('----------> can not generate shuffled data or shuffle mode not available; shuffle significance not computed');
        end
        %%%%third: compute decoding score (matchscore) for each shuffle
        if ok 
           for (i = 1:nshuffle) %%%forget the ep{mm,nn} structure here: orginal event times are in evTimes
                shufposPnow = cell(mm,nn); shufnegPnow = cell(mm,nn); nactnow = cell(mm,nn);
                for (ii = 1:mm)
                for (jj = 1:nn)
                     [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, shufposPnow{ii,jj}, shufnegPnow{ii,jj}, ~, nactnow{ii,jj}] = ...
                        findallpositionprob(dataok, shufcellrate{i}, Xgrid, seqparm, ep{ii,jj}, [], []);
                end
                end
                shufposP{tti,i} = cell2mat(shufposPnow); shufnegP{tti,i} = cell2mat(shufnegPnow);
                nact{tti,i} = cell2mat(nactnow);
           end
        end
      end
    end
end
if ok
for i = 1:nshuffle
    shufposPnow = NaN*ones(ntmp,mm*nn); shufnegPnow = NaN*ones(ntmp,mm*nn); nactnow = NaN*ones(ntmp,mm*nn);
    for (tti = 1:ntmp)
        if (numel(shufposP{tti,i}) == mm*nn)
            shufposPnow(tti,:) = shufposP{tti,i}; shufnegPnow(tti,:) = shufnegP{tti,i}; nactnow(tti,:) = nact{tti,i};
        else
            disp('----------> event lap number not match; results from a template may be ignored');
        end
    end
    [shufposN(i), shufnegN(i), shufmatchN(i)] = findmatchnumbernow(nactnow, shufposPnow, shufnegPnow, seqparm.significancelevel, seqparm.minactivecell);
    %disp([shufposN(i) shufnegN(i) shufmatchN(i)]);
end
posP = numel(find(shufposN>=posN))/nshuffle; negP = numel(find(shufnegN>=negN))/nshuffle; matchP = numel(find(shufmatchN>=matchN))/nshuffle; 
shufposN = shufposN(~isnan(shufposN)); shufnegN = shufnegN(~isnan(shufnegN)); shufmatchN = shufmatchN(~isnan(shufmatchN));
if ~isempty(shufposN) posZ = (posN-mean(shufposN))/std(shufposN); end
if ~isempty(shufnegN) negZ = (negN-mean(shufnegN))/std(shufnegN); end
if ~isempty(shufmatchN) matchZ = (matchN-mean(shufmatchN))/std(shufmatchN); end
end

function [prob, timepoint, seqstart, seqend, seqmarker, matchscore, posmatchprob, negmatchprob, ...
                 shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decodeerr, nactivecell] = findallpositionprob(spiketime, cellrate, Xgrid, seqparm, ep, realtime, realloc)
decparm.winsize = seqparm.windowsize; decparm.shifttime = seqparm.shifttime;
[timepoint, prob] = DE_Seq_DecodePositionProbability(spiketime, Xgrid, cellrate, ep, decparm);
seqstart = ep.start; seqend = ep.ent; seqmarker = ep.marker;
%%%%%%%%Need to work out the rest - quantification and significance
nactivecell = findactivecell(spiketime, ep);
[matchscore, posmatchprob, negmatchprob, ...
                 shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decodeerr] = finddecodingquantifications(timepoint, prob, Xgrid, seqparm, realtime, realloc);

function nact = findactivecell(spiketime, ep)
nep = numel(ep.start); ncell = numel(spiketime); cellact = zeros(ncell, nep); 
nact = zeros(1, nep);
for (i = 1:nep)
    for (j = 1:ncell)
        if numel(find( (spiketime{j}>=ep.start(i)) & (spiketime{j}<=ep.ent(i)) )) > 0
            cellact(j,i) = 1;
        end
    end
end
nact = sum(cellact);
      
function [tmpnow, ok] = findtemplates(tmp, cellfilename)
tmpnow = []; ok = 0;
for (i = 1:numel(tmp))
    tmpfilename = tmp(i).tmp.tmpfilename;
    if (numel(tmpfilename) == numel(cellfilename))
        allok = zeros(1, numel(cellfilename));
        for (j = 1:numel(cellfilename))
            if ~isempty(find(strcmp(tmpfilename, cellfilename{j}))) allok(j) = 1; end
        end
        if isempty(find(allok==0))
            tmpnow = tmp(i); ok =1; break
        end
    end
end


function [evSess, evsessid] = identifysession(evTime, sessionname, startT, endT)
evSess = []; evsessid = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii}; evsessid = iii;
end

function [matchscore, posmatchprob, negmatchprob, ...
                 shufZscore, shufZposP, shufZnegP, shufposP, shufnegP, decodeerr] = finddecodingquantifications(timepoint, prob, Xgrid, seqparm, alltime, allloc)
%%%%%timepoint{nep}(ntimepoint); prob{nep}{ntimepoint}(nXgrid); Xgrid(nXgrid) 
%%%%%%%%%%%determine the decoded position, then compute: decode error; pearson R, p; shuffle Z, p 
%%%%%%%%%%%now only for 1D position decoding
nshuffle = seqparm.Nwithineventshuffle; nep = numel(timepoint); timebinsize = seqparm.shifttime;
matchscore = NaN*ones(1, nep); posmatchprob = NaN*ones(1, nep); negmatchprob = NaN*ones(1, nep); decodeerr = NaN*ones(1, nep);
shufZscore = NaN*ones(1, nep); shufZposP = NaN*ones(1, nep); shufZnegP = NaN*ones(1, nep); shufposP = NaN*ones(1, nep); shufnegP = NaN*ones(1, nep);
Xgrid = Xgrid{1}; %%%%for now: only do it for 1D space
for (i = 1:nep)
    %%%%%determine decoded positions
    timenow = timepoint{i}; ntime = numel(timenow); locnow = NaN*ones(size(timenow));
    for (j = 1:ntime) 
        %disp( [size(prob{i}{j}) size(Xgrid)] );
        [mm, iii] = max(prob{i}{j}); 
        if ~isnan(mm) locnow(j) = Xgrid(iii); end
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
              jj = find( (alltime{i}>=timenow(ii)) & (alltime{i}<timenow(ii)+timebinsize) );
              if (numel(jj) >= 1) 
                  aaa = allloc{i}(jj); aaa = aaa(~isnan(aaa));
                  realloc(ii) = mean(aaa); 
              end
          end
          dis = abs(locnow-realloc); 
          decodeerr(i) = mean(dis(~isnan(dis)));
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

%%%%%%%%%%%%%%%%%%%Below is for computing shuffle significance of sequence matching results%%%%%%%%%%%%%%%
function [posP, posE, negP, negE, posSig, negSig] = findseqmatchingshufflesig(tmpseq, npos, nneg, sessST, sessET, fileletter, evseq, rrank, seqparm, sessev, tmprank, eventoption) %, alltime)
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
   [~, ~, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_gloabalID_seqmatch(tmpseq, fileletter, evseq, seqparm, rrank, tmprank, eventoption);%, alltime);
elseif strcmp(seqparm.shufflemode, 'LocalID')
   [~, ~, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_localID_seqmatch(tmpseq, fileletter, evseq, seqparm, rrank, tmprank, eventoption);%, alltime);
elseif strcmp(seqparm.shufflemode, 'CircleSlide')
   [~, ~, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_circleslide_seqmatch(tmpseq, fileletter, evseq, seqparm, ...
       sessST, sessET, rrank, tmprank, eventoption, timingoption);
else
    ok = 0; disp('----------> can not generate shuffled data or shuffle mode not available; shuffle significance not computed');
end
if ok
   shufnpos = zeros(1,nshuffle); shufnneg = zeros(1,nshuffle);
   for (k = 1:nshuffle)
     [shufnpos(k), shufnneg(k)] = findmatchnum_seqmatch(posmatchprob{k}, negmatchprob{k}, seqstart{k}, seqend{k}, slevel, eventoption); 
   end
   posP = numel(find(shufnpos>=npos))/nshuffle; negP = numel(find(shufnneg>=nneg))/nshuffle;
   posE = mean(shufnpos); negE = mean(shufnneg); posSig = std(shufnpos); negSig = std(shufnneg);
end
function [seq, matchscore, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_gloabalID_seqmatch(tmpseq, fileletter, evseq, seqparm, rrank, tmprank, eventoption) %alltime
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
function [seq, matchscore, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_localID_seqmatch(tmpseq, fileletter, evseq, seqparm, rrank, tmprank, eventoption) %, alltime)
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
function [seq, matchscore, posmatchprob, negmatchprob, seqstart, seqend] = generateshuffle_circleslide_seqmatch(tmpseq, fileletter, evseq, seqparm, sessST, sessET, rrank, tmprank, eventoption, timingoption)
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
function [npos, nneg] = findmatchnum_seqmatch(evposprob, evnegprob, evseqstart, evseqend, slevel, eventoption)
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
