function DM_Seq_ShowSequenceLocations_Callback
%%% Display where (current locations) and what (decoded locations) sequences detected in a .seqdb
%%%    only show sequences associated with one template (may be multiple events)

hf = gcbf; pinfonow = getappdata(hf, 'pinfo'); datanow = getappdata(hf, 'data'); tagnow = get(gcbo, 'Tag');
plotparm = getappdata(hf, 'plotparm');
hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); cellind = find(spikeselection==1);
%hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
%for (kk = 1:numel(grpind)) cellind = union(cellind, datanow.grouplist.groupindex{grpind(kk)}); end
ok = 1;
disp('-----> Displaying target locations ......');
if (plotparm.evselect == 0) %%%select sessions
      disp('-----------> shuffle significance cannot be computed: need to select the event option'); ok = 0;
elseif ~isfield(datanow, 'data')
      ok = 0; disp('-----------> aborted: data not available for stripped databases');
end
if ok
   if (numel(cellind)>1) disp('-----------> more than 1 items selected; only display the first selection'); end
if (plotparm.linkbehav == 0)
       disp(['-----------> no behav data linked']); ok = 0;
else
       behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
end
cellind = cellind(1); %%%%for sequence disiplay: only display the  first selected template
if (plotparm.evselect == 0) %%%select sessions
    seqdata = datanow.data.sessseq{cellind}; 
    evnames = pinfonow.general.sessionname{cellind}; 
    eventoption = pinfonow.parm.sessionoption{cellind}; 
    timingoption = pinfonow.parm.sessiontimingoption{cellind};
else
    seqdata = datanow.data.evseq{cellind}; someparm.eventoption = pinfonow.parm.eventoption{cellind};
    someparm.pthre = []; someparm.sthre = []; someparm.maxstep = [];
    if isfield(pinfonow.parm, 'probthres') someparm.pthre = pinfonow.parm.probthres{cellind}; end
    if isfield(pinfonow.parm, 'speedthres') someparm.sthre = pinfonow.parm.speedthres{cellind}; end
    if isfield(pinfonow.parm, 'maxstepdis') someparm.maxstep = pinfonow.parm.maxstepdis{cellind}; end
    evnames = pinfonow.general.eventname{cellind}; evType = pinfonow.parm.eventType{cellind}; 
    tmpID = pinfonow.general.tmpID{cellind}; animaldate = strcat(pinfonow.general.animalname{cellind}, '_', pinfonow.general.datedir{cellind});
end
end
if ok
    seqtype = [];
    if (isfield(pinfonow.parm, 'seqtype')) someparm.seqtype = pinfonow.parm.seqtype{cellind}; end
    if contains(someparm.seqtype, 'evtitemized') 
        nev = numel(tmpID); evnames = erase(evnames, [animaldate '_']); 
        title = strcat('Event-', evnames, '___Tmp-', tmpID);
    else
        nev = numel(evnames); title = strcat('Tmp-', tmpID, '___Event-', evnames);
    end
end
if ok 
    evtind = [];
    if contains(someparm.seqtype, 'Bayesian') && (~contains(someparm.seqtype, '2D')) %%%choices provided only for 1D Bayesian
       varlist = {'Positive (forward)'; 'Negative (reverse)'; 'All targets (1D forw/rev)'};
       [ind, ok] = listdlg('ListString', varlist, 'Name', 'What sequences');
       if ~ok
           disp('-----------> no plot type selected or cancelled');
       else
           someparm.sstype = varlist{ind(1)};
           evtind = findidentifiedevents(ind(1), datanow, cellind);
       end
    elseif strcmp(seqtype, 'Bayesian2D')
       someparm.sstype = 'All 2D targets';
       evtind = findidentifiedevents(4, datanow, cellind);
    end
    if isempty(evtind)
        ok=0; disp('-----------> no target sequences are found; aborted');
    end
end
if ok
    rrmode = []; if isfield(pinfonow.parm, 'rankmode') rrmode = pinfonow.parm.rankmode{cellind}; end 
    parentfilename = []; if isfield(datanow, 'parentfile') parentfilename = datanow.parentfile; end
    tmpevname = pinfonow.tmp.sessevname{cellind}; tmpfilename = pinfonow.tmp.filename{cellind};
    rrank = []; rfile = [];
    [MCroot, ~, ~, ~, ~, ~] = CurrentVersion;
    if strcmp(rrmode, 'PairRank')
        rfile = fullfile(MCroot, 'DataExplorer', 'Sequence', 'rrankprob.mat'); 
    elseif strcmp(rrmode, 'OrderRank')
        rfile = fullfile(MCroot, 'DataExplorer', 'Sequence', 'orderrankprob.mat');
    end
    if (exist(rfile, 'file') == 2)
       S = load(rfile); rrank = S.rrank; S = [];
    elseif strcmp(rrmode, 'PairRank')
       ok = 0; disp('-----> pair ranking probability file for selected rank mode not found; aborted');
    end
end
if ok
   for (i = 1:nev)
       evseqnow = seqdata{i}; %%%for the current event 
       tagstr{1} = strcat('Identified sequence type: ', someparm.sstype);
       evT = datanow.events.eventimes{cellind}{i};
       if contains(someparm.seqtype, 'evtitemized')
           plotcurrentdecodingevents(evseqnow, seqtype, behav, bhdata, evtind{i}, evnames, evT, cellind, pinfonow, title{i}, tagstr, someparm);
       else
           plotcurrentdecodingevents(evseqnow, seqtype, behav, bhdata, evtind{i}, evnames{i}, evT, cellind, pinfonow, title{i}, tagstr, someparm);
       end
   end
end
disp('*******************************');

function plotcurrentdecodingevents(evseq, seqtype, behav, bhdata, evtind, evName, evT, cellid, pinfonow, title, tagstr, someparm)
%%%% evseq.(fieldnames){1:nlap,1}{1 or ntarget} for each current event file
%%%% evtind[targetlaps indices]
ok = 1; finaldirnow = pinfonow.general.finaldir{cellid};
%%%% The following are what need to be assigned
Xgrid = evseq.Xgrid;
ntarg = numel(evtind); %%%[mm,nn] = size(evseq.timepoint); for now: mm = nlap, nn = 1
% disp(evtind)
% disp(numel(evseq.timepoint))
% disp('*****')
%%%%%%%need to rearrange evseq.timepoint and evseq.quantdata for 1D Bayesian free sequencing
timepointnow = evseq.timepoint; quantdatanow = evseq.quantdata;
if strcmp(someparm.eventoption, 'FreeSequences')
    [timepointnow, quantdatanow] = rearange(evseq.quantdata); %{evtind(ii),1}{1}; evseq.quantdata{evtind(ii),1}
end
timepoints = cell(ntarg, 1); decodelocx = cell(ntarg, 1); decodelocy = cell(ntarg, 1); decode1dloc = cell(ntarg, 1);
realtime = cell(ntarg, 1); reallocx = cell(ntarg, 1); reallocy = cell(ntarg, 1); real1dloc = cell(ntarg, 1); 
%%%%%%%decoded timepoints and decoded locations
for (ii = 1:ntarg)
    timepoints{ii} = timepointnow{evtind(ii),1}{1}; decode1dloc{ii} = NaN*ones(size(timepoints{ii}));
    decodelocx{ii} = NaN*ones(size(timepoints{ii})); decodelocy{ii} = NaN*ones(size(timepoints{ii})); 
end
if contains(someparm.seqtype, '2D') %%% if 2D sequences
    for (ii = 1:ntarg)
         decodelocx{ii} = quantdatanow{evtind(ii),1}.peakx{1}; decodelocy{ii,1} = quantdatanow{evtind(ii),1}.peaky{1}; 
         %%%decode1dloc{ii} = evseq.quantdata{evtind(ii),1}.peak1dx{1}; this needs to be decoded out below from peakx, peaky
    end
else %%%%if 1D decoding
    %if strcmp(someparm.eventoption, 'FreeSequences')
       
    %else
       for (ii = 1:ntarg)
            decode1dloc{ii} = quantdatanow{evtind(ii),1}.peakloc{1};  %disp([numel(timepoints{ii}) numel(decode1dloc{ii})]);
       end
    %end
end
%%% Real 2D and 1D (if possible) positions for the decoded timepoints
evid = []; pout = [];
[evsessname, ~] = identifysession(evT, pinfonow.general.sessionname{cellid}, pinfonow.general.sessionstartT{cellid}, pinfonow.general.sessionendT{cellid});
posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evsessname) );
if numel(posid) == 1
    %%%%2D position data 
    postimestamp = bhdata.pos.postimestamp{posid}*behav.parm.timeunit(posid);  
    allposmarker = behav.general.posMarker{posid}; Pmarker = behav.parm.sessPmarker{posid}; 
    posXX = []; posYY = []; ik = find(strcmp(allposmarker, Pmarker)); 
    if (numel(ik) == 1) posXX = bhdata.pos.XX{posid}{ik}*behav.parm.pixelXSize(posid); posYY = bhdata.pos.YY{posid}{ik}*behav.parm.pixelYSize(posid); end
    [realtime, reallocx, reallocy] = locatepositionnow(timepoints, postimestamp, posXX, posYY);
    %%%for 1D position data
    if (isfield(behav.general, 'eventname'))
        evid = find(strcmp(behav.general.eventname{posid}, evName));
    else
        evid = find(strcmp(behav.behavior.eventname{posid}, evName));
    end
else 
    disp(['-------------> warning: position data for the event not found: ', finaldirnow, '___', evName]);
end
if (numel(evid)==1)
    %%%%have to linearize the current postion on the spot 
    joint = locatejoint(behav, bhdata, posid, evid);
    if ~isempty(joint)
       joint(:,1) = joint(:,1)*behav.parm.pixelXSize(posid); joint(:,2) = joint(:,2)*behav.parm.pixelYSize(posid);
       spaceunit = (behav.parm.pixelXSize(posid)+behav.parm.pixelYSize(posid))/2;
       [real1dloc, pout] = linearizenow(joint, realtime, reallocx, reallocy, spaceunit);
       %%%%linearize decoded 2D locations
       if strcmp(someparm.seqtype, 'Bayesian2D')
           for (ii = 1:ntarg)
               [xout, ~, pout] = LinearizePathXY_new(joint, decodelocx{ii}, decodelocy{ii}, 0.8, 0, spaceunit);
               decode1dloc{ii} = xout;
           end
       end
    else
       disp(['-------------> warning: linearization file not found in behavior database: ', finaldirnow, '___', evName]); 
    end
else
    disp(['-------------> warning: event name not found in behavior database: ', finaldirnow, '___', evName]);
end
if ok 
   %%%% plot out data: (i) 2D real, (ii) 2D decoded locations, (iii) 1D real, (iv) 1D decoded locations
   if ~isempty(joint)
      plot1Dpositions(realtime, real1dloc, Xgrid, pout, title, tagstr, 'Real 1D location (cm)', [0 0 0]); 
      plot1Dpositions(timepoints, decode1dloc, Xgrid, pout, title, tagstr, 'Decoded 1D location (cm)', [0.8 0 0.6]);
   end
   plot2Dpositions(realtime, reallocx, reallocy, Xgrid, title, tagstr, 'Real 2D location (cm)', [0 0 0]); 
   if strcmp(seqtype, 'Bayesian2D')
      plot2Dpositions(timepoints, decodelocx, decodelocy, Xgrid, title, tagstr, 'Decoded 2D location (cm)', [0.8 0 0.6]);
   end
end

function plot2Dpositions(timepoint, locx, locy, Xgrid, title, tagstr, labelstr, col)
if iscell(Xgrid) && numel(Xgrid)==2
    minx = 0.9*min(Xgrid{1}); maxx = 1.1*max(Xgrid{1}); 
    miny = 0.9*min(Xgrid{2}); maxy = 1.1*max(Xgrid{2});
else
    [minx, maxx] = determinelimits(locx); [miny, maxy] = determinelimits(locy);
end
hf = figure('Name', title); hax = axes('Parent', hf, 'XLim', [minx maxx], 'YLim', [miny maxy], 'NextPlot', 'add'); 
xlabel (labelstr); ylabel(labelstr);
for (i = 1:numel(timepoint))
%     size(loc{i})
%     size(timepoint{i})
%     disp('***********');
    line(hax, locx{i}, locy{i}, 'Color', col, 'Linewidth', 0.2, 'Marker', '.', 'MarkerSize', 5, 'MarkerFaceColor', col, 'MarkerEdgeColor', col); 
end
%%%% start and end positions
for (i = 1:numel(timepoint))
    line(hax, locx{i}(1), locy{i}(1), 'LineStyle', 'none', 'Marker', '>', 'MarkerSize', 10, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]); 
    line(hax, locx{i}(numel(locx{i})), locy{i}(numel(locy{i})), 'LineStyle', 'none', 'Marker', '<',...
        'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]); 
end
% for (j = 1:numel(pout))
%     line(hax, [pout mint], [pout maxt], 'Color', [0.5 0.5 0.5], 'Linewidth', 0.2);
% end
ind = strfind(title, '___'); str1 = title(ind+3:numel(title));
text('Parent', hax, 'String', str1, 'Interpreter', 'none', 'Units', 'normalized', 'Position', [0.05 0.96]);
for (i = 1:numel(tagstr))
    text('Parent', hax, 'String', tagstr{i}, 'Interpreter', 'none', 'Units', 'normalized', 'Position', [0.05 0.92-(i-1)*0.04]);
end

function plot1Dpositions(timepoint, loc, Xgrid, pout, title, tagstr, labelstr, col)
%xlimnow = [0 1.5*max(Xgrid{1})]; %xlimnow = [0.9*min(Xgrid{1}) 1.1*max(Xgrid{1})]; 
%----Be careful here: Xgrid for 1D decoding is the linearized path, but for 2D is just the real data of first dimension
if ~isempty(pout) 
    xlimnow = [-5 1.05*max(pout)];
else
    xlimnow = [-5 1.5*max(Xgrid{1})];
end
[mint, maxt] = determinelimits(timepoint);
hf = figure('Name', title); hax = axes('Parent', hf, 'XLim', xlimnow, 'YLim', [mint maxt], 'NextPlot', 'add'); 
xlabel (labelstr); ylabel('Time (s)');
for (i = 1:numel(timepoint))
    %disp([numel(timepoint{i}) numel(loc{i})]);
    line(hax, loc{i}, timepoint{i}, 'Color', col, 'Linewidth', 0.2, 'Marker', '.', 'MarkerSize', 10, 'MarkerFaceColor', col, 'MarkerEdgeColor', col); 
end
%%%% start and end positions
for (i = 1:numel(timepoint))
    line(hax, loc{i}(1), timepoint{i}(1), 'LineStyle', 'none', 'Marker', '>', 'MarkerSize', 10, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]); 
    line(hax, loc{i}(numel(loc{i})), timepoint{i}(numel(loc{i})), 'LineStyle', 'none', 'Marker', '<',...
        'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]); 
end
for (j = 1:numel(pout))
    line(hax, [pout(j) pout(j)], [mint maxt], 'Color', [0.8 0.8 0.8], 'Linewidth', 0.5);
end
ind = strfind(title, '___'); str1 = title(ind+3:numel(title));
text('Parent', hax, 'String', str1, 'Interpreter', 'none', 'Units', 'normalized', 'Position', [0.05 0.96]);
for (i = 1:numel(tagstr))
    text('Parent', hax, 'String', tagstr{i}, 'Interpreter', 'none', 'Units', 'normalized', 'Position', [0.05 0.92-(i-1)*0.04]);
end

function [evSess, evsessid] = identifysession(evTime, sessionname, startT, endT)
evSess = []; evsessid = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii}; evsessid = iii;
end

function joint = locatejoint(behav, bhdata, sessid, evid)  
joint = [];
posltrfile = bhdata.pos.ltrfilename{sessid}; posltr = bhdata.pos.posltr{sessid};
evPosltr = behav.parm.eventPosltr{sessid};
for (tt = 1:numel(posltrfile))
    if (strcmp(posltrfile{tt}, evPosltr{evid}))
       joint = posltr{tt}; break
    end
end

function [mint, maxt] = determinelimits(timepoint)
minnow = zeros(1, numel(timepoint)); maxnow = zeros(1, numel(timepoint));
for (i = 1:numel(timepoint))
   minnow(i) = min(timepoint{i}); maxnow(i) = max(timepoint{i}); 
end
mint = 0.98*min(minnow); maxt = 1.02*max(maxnow); 
if isempty(mint) mint = 0; maxt = 1; end
if ~(maxt>mint) maxt = mint +1 ; end

function evtind = findidentifiedevents(choiceind, seqdata, cellind)
%%% varlist = {'Positive (forward)'; 'Negative (reverse)'; 'All matching (forward & reverse)'};
%%% evtind{nev}
switch choiceind
    case 1
          evtind = seqdata.events.posevind{cellind}; 
    case 2
          evtind = seqdata.events.negevind{cellind};
    case 3
          evtind = seqdata.events.matchevind{cellind};
    case 4
          evtind = seqdata.events.targevind{cellind};
end

function [realtime, reallocx, reallocy] = locatepositionnow(timepoints, postimestamp, xpos, ypos)
nlap= numel(timepoints); realtime = cell(nlap, 1); reallocx = cell(nlap, 1); reallocy = cell(nlap, 1);
for (tt = 1:nlap)
    if ~isempty(timepoints{tt})
       evposind = find( (postimestamp>=min(timepoints{tt})) & (postimestamp<=max(timepoints{tt})) );
       realtime{tt} = postimestamp(evposind); reallocx{tt} = xpos(evposind); reallocy{tt} = ypos(evposind);
    end
end
function [real1dloc, pout] = linearizenow(joint, realtime, xpos, ypos, spaceunit)
real1dloc = cell(size(realtime)); pout = [];
for (tt = 1:numel(realtime))
    if ~isempty(realtime{tt})
       [xout, ~, pout] = LinearizePathXY_new(joint, xpos{tt}, ypos{tt}, 0.8, 0, spaceunit);
       real1dloc{tt} = xout; 
    end
end

function [timepointout, Pout] = rearange(quantdata) %For free sequencing only
%%% features already rearranged: evseeq.shufposP{lapind,1}(ntarget) -> shufposP(evind) = cell2mat(shufposP{ntarget}) 
%%% timepoint needs to be re-arranged back to the old format: from evseq.quantdata{nlap, 1}.(field){ntargetnow} to evseq.quantdata{ntotaltarget, 1}.(field){1}
ntarget = 0; [mm,nn] = size(quantdata); fnames = fieldnames(quantdata{1,1}); fnames = setdiff(fnames, {'timepoint'});
for (i = 1:mm)
    for (j = 1:nn)
        ntarget = ntarget + numel(quantdata{mm,nn}.timepoint);
    end
end
timepointout = cell(ntarget, 1); Pout = cell(ntarget, 1); n = 0;
% for (i = 1:ntarget)
%     timepointout{i,1} = cell(1,1);
%     for (j = 1:numel(fnames)) Pout{i,j}.(fnames) = cell(1,1); end
% end
for (i=1:mm)
    for (j=1:nn)
        for (k = 1:numel(quantdata{i,j}.timepoint))
             n = n+1; timepointout{n,1}{1} = quantdata{i,j}.timepoint{k};
             for (tj = 1:numel(fnames))
                 Pout{n,1}.(fnames{tj}){1} = quantdata{i,j}.(fnames{tj}){k};
             end
        end
    end
end




    