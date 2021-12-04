function DM_Seq_PlotSequenceVelocity
%%%based on sequencing results of a .seqdb (for the Tau mouse project)
%%%      plot how animal velocity varies within the sequence and compare with the whole session velocity
%%%          -- do this for track running events and significant open/sleep sessions
%%%      For each event/session of a group:
%%%          -- compute a zscore of the mean sequence speed against the speed distrution of the entire event/session
%%%          -- plot the distribution of the sequence speed against the entire speed distribution
ok = 1;
hf = gcbf; pinfonow = getappdata(hf, 'pinfo'); datanow = getappdata(hf, 'data');
plotparm = getappdata(hf, 'plotparm');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); grpind = find(groupselection == 1); 
ngroup = numel(grpind); cellind = cell(1, ngroup); grpnames = cell(1, ngroup);
for (kk = 1:numel(grpind)) 
    cellind{kk} = datanow.grouplist.groupindex{grpind(kk)}; grpnames{kk} = datanow.grouplist.groupname{grpind(kk)}; 
end
if (plotparm.linkbehav == 0);
    disp(['--------> no behavioral data linked']); ok = 0;
else
    behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
end

if ok
plotdistflag = 0;
ss = questdlg('Plot velocity distribution for all the sessions?');
if (strcmp(ss, 'Yes'))
    plotdistflag = 1;
end
evtraj = {'selfTrack'; 'open'; 'sleep'}; ntraj = numel(evtraj);
evsigval = 1; %1.645;
plotval = cell(ngroup, ntraj); mm = NaN*ones(ngroup, ntraj); se = NaN*ones(ngroup, ntraj); 
md = NaN*ones(ngroup, ntraj);  nn = NaN*ones(ngroup, ntraj); 
ppgroup = NaN*ones(ngroup-1, ntraj); pptraj = NaN*ones(ngroup, ntraj-1);
for (i = 1:ngroup) %%%%self track running events
    tmpname = pinfonow.tmp.sessevname(cellind{i}); score = pinfonow.seq.evPosMatchNshufsig(cellind{i});
    evname = pinfonow.general.eventname(cellind{i}); evtype = pinfonow.parm.eventType(cellind{i});
    for (k = 1:numel(cellind{i}))
        [S, T] = strtok(tmpname{k}, '_'); S = lower(S); T = lower(T);
        evnow = evname{k}; evtypenow = evtype{k}; scorenow = score{k};
        sessyes = zeros(1, numel(evnow)); trajyes = zeros(1, numel(evnow)); 
        typeyes = zeros(1, numel(evnow)); scoreyes = zeros(1, numel(evnow));
        for (j = 1:numel(evnow))
                if (~isempty(strfind(lower(evnow{j}), S))) sessyes(j) = 1; end
                if (~isempty(strfind(lower(evnow{j}), T))) trajyes(j) = 1; end
                if strcmp(evtypenow{j}, 'run') typeyes(j) = 1; end
                if (scorenow(j) > evsigval) scoreyes(j) = 1; end
        end
        ii = find(sessyes & trajyes & typeyes & scoreyes); 
        if (numel(ii) == 1) 
            aa = findeventmeansequencevelocity(datanow, cellind{i}(k), ii, plotparm.significancelevel, ...
                    pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, evname{k}{ii}, plotdistflag);  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
            if (~isnan(aa))
                plotval{i,1} = [plotval{i,1} aa];
            end
        end
    end
end
for (i = 1:ngroup) %%%%open/sleep sessions
    tmpname = pinfonow.tmp.sessevname(cellind{i}); score = pinfonow.seq.sessPosMatchNshufsig(cellind{i});
    sessname = pinfonow.general.sessionname(cellind{i}); sesstype = pinfonow.parm.sessType(cellind{i});
    for (k = 1:numel(cellind{i}))
        sessnow = sessname{k}; sesstypenow = sesstype{k}; scorenow = score{k};
        for (j = 1:numel(sessnow))
             if (strcmp(sesstypenow{j}, 'open')) && (scorenow(j) > evsigval) 
                aa = findsessmeansequencevelocity(datanow, cellind{i}(k), j, plotparm.significancelevel, ...
                    pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, sessname{k}{j}, plotdistflag);  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
                if (~isnan(aa))
                    plotval{i,2} = [plotval{i,2} aa];
                end
             elseif (strcmp(sesstypenow{j}, 'sleep')) && (scorenow(j) > evsigval) 
                aa = findsessmeansequencevelocity(datanow, cellind{i}(k), j, plotparm.significancelevel, ...
                    pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, sessname{k}{j}, plotdistflag);  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
                if (~isnan(aa))
                    plotval{i,3} = [plotval{i,3} aa];
                end
             end
        end
    end
end
for (i = 1:ngroup) %%%%basic statistics
        for (j = 1:numel(evtraj))
            aa = plotval{i,j}; aa = aa(~isnan(aa)); plotval{i,j} = aa;
            nn(i,j) = numel(aa);
            if (nn(i,j) > 1)
                mm(i,j) = mean(aa); se(i,j) = std(aa)/sqrt(nn(i,j)); md(i,j) = median(aa);
            end
        end
end
for (i = 1:ngroup-1) %%%%group comparasons: all compare with the first one
        for (j = 1:ntraj)
            if (~isempty(plotval{1,j})) && (~isempty(plotval{i+1,j}))
               [~, ppgroup(i,j)] = ttest2(plotval{1,j}, plotval{i+1,j});
            end
        end
end
for (j = 1:ntraj-1) %%%%traj comparasons: all compare with the first one
        for (i = 1:ngroup)
            if (~isempty(plotval{i,1})) && (~isempty(plotval{i,j+1}))
               [~, pptraj(i,j)] = ttest2(plotval{i,1}, plotval{i,j+1});
            end
        end
end

hf = figure('Name', ['SignificantSequences - velocity']); wid = 1/ngroup; gap = 0.05;
for (i = 1:ngroup)
        hax(i) = axes('Parent', hf, 'NextPlot', 'add', 'Units', 'normalized', 'Position', [(i-1)*wid+gap 0.05 wid-gap 0.95]);
        for (j = 1:ntraj)
            YY = plotval{i,j}; XX = j*ones(1, nn(i,j)); 
            line(XX, YY, 'Parent', hax(i), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30);
            line([j-0.2 j+0.2], [mm(i,j) mm(i,j)], 'Parent', hax(i), 'LineWidth', 2); %%%mean value line 
        end
        line([0 ntraj+1], [evsigval evsigval], 'Parent', hax(i)); %%%significance line
        str = [grpnames{i}, ': selftrack, open, sleep'];
        text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.96]);
        str = ['nn=[', num2str(nn(i,:)), ']; mean=[', num2str(mm(i,:)), ']']; 
        text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.92]);
        str = ['se=[', num2str(se(i,:)), ']; median=[', num2str(md(i,:)), ']'];
        text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.88]);
        str = ['sess/ev comparisons (t-test against the 1st) p=', num2str(pptraj(i,:))]; %, '%10.5e')];
        text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.84]);
        thnow = 0.84;
        for (k = 1:ntraj) %%%%group comparisons
            thnow = thnow - 0.04;
            str = [evtraj{k}, ' group comparisons (t-test against the 1st) p=', num2str(ppgroup(:,k)')]; %, '%10.5e')];
            text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 thnow]);
        end
end
end
disp('*******************************');

function val = findsessmeansequencevelocity(datanow, cellind, sessind, siglevel,...
    pinfo, data, behav, bhdata, finaldirnow, sessnamenow, plotdistflag)
val = NaN;
seqnow = datanow.data.sessseq{cellind}{sessind};
seqstart = cell2mat(seqnow.seqstart); seqend = cell2mat(seqnow.seqend); 
posmatchprob = cell2mat(seqnow.posmatchprob); %matchscore = seqnow.matchscore;  negmatchprob = seqnow.negmatchprob;
iii = find( posmatchprob < siglevel ); 
posmatchprob = posmatchprob(iii); seqstart = seqstart(iii); seqend = seqend(iii); %%%only the significant sequences
[newstart, newend, ~, ~] = resolveoverlap(posmatchprob, seqstart, seqend); %%resolve the overlaps
 %%%this contains all the center times of the significant sequences
%%%%locate event position data
posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, sessnamenow) );
if numel(posid == 1)
    postimestamp = bhdata.pos.postimestamp{posid}*behav.parm.timeunit(posid);
    allvel = bhdata.sess.AllVel{posid};
    %%%%%%%%%sequence velocity
    nseq = numel(newstart); seqvel = NaN*ones(1, nseq);
    for (i = 1:nseq)
        iii = find( (postimestamp>=newstart(i)) & (postimestamp<=newend(i)) );
        seqvel(i) = mean(allvel(iii));
    end
    seqvel = seqvel(~isnan(seqvel));
    %%%%%%%%%sess velocity
    allvel = allvel(~isnan(allvel));
    if (~isempty(seqvel)) && (~isempty(allvel)) 
       %val = (mean(seqvel)-mean(allvel))/std(allvel);
       val = mean(seqvel)/mean(allvel);
       %%%%plot velocity distribution
       if plotdistflag %%%if plot and compare the distribution
          plotveldistribution(allvel, seqvel, finaldirnow, sessnamenow);
       end
    end
end

function val = findeventmeansequencevelocity(datanow, cellind, evind, siglevel,...
    pinfo, data, behav, bhdata, finaldirnow, evnamenow, plotdistflag)
val = NaN;
seqnow = datanow.data.evseq{cellind}{evind};
seqstart = cell2mat(seqnow.seqstart); seqend = cell2mat(seqnow.seqend); 
posmatchprob = cell2mat(seqnow.posmatchprob); %matchscore = seqnow.matchscore;  negmatchprob = seqnow.negmatchprob;
iii = find( posmatchprob < siglevel ); 
posmatchprob = posmatchprob(iii); seqstart = seqstart(iii); seqend = seqend(iii); %%%only the significant sequences
[newstart, newend, ~, ~] = resolveoverlap(posmatchprob, seqstart, seqend); %%resolve the overlaps
 %%%this contains all the center times of the significant sequences
%%%%locate event position data
evTime = data.events.eventimes{cellind}{evind};
evSess = identifysession(evTime, pinfo.general.sessionname{cellind}, pinfo.general.sessionstartT{cellind}, pinfo.general.sessionendT{cellind});
posid = []; 
if (~isempty(evSess))
    posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess) );
end
if numel(posid == 1)
    postimestamp = bhdata.pos.postimestamp{posid}*behav.parm.timeunit(posid);
    allvel = bhdata.sess.AllVel{posid};
    %%%%%%%%%sequence velocity
    nseq = numel(newstart); seqvel = NaN*ones(1, nseq);
    for (i = 1:nseq)
        iii = find( (postimestamp>=newstart(i)) & (postimestamp<=newend(i)) );
        seqvel(i) = mean(allvel(iii));
    end
    seqvel = seqvel(~isnan(seqvel));
    %%%%%%%%%lap velocity
    lapvel = [];
    for (i = 1:numel(evTime.start))
        iii = find( (postimestamp>=evTime.start(i)) & (postimestamp<=evTime.ent(i)) );
        lapvel = [lapvel allvel(iii)];
    end
    lapvel = lapvel(~isnan(lapvel));
    if (~isempty(seqvel)) && (~isempty(lapvel)) 
       %val = (mean(seqvel)-mean(lapvel))/std(lapvel);
       %val = mean(seqvel)/mean(lapvel);
       val = mean(seqvel)/median(lapvel);
       %%%%plot velocity distribution
       if plotdistflag %%%if plot and compare the distribution
          plotveldistribution(lapvel, seqvel, finaldirnow, evnamenow);
       end
    end
end

function plotveldistribution(allvel, seqvel, fdir, evname)
nbin = 100; maxV = max(allvel); binsize = 20; %maxV/nbin;
X = 0:binsize:maxV; Y = histc(allvel, X); 
hf = figure('Name', 'Sequence velocity'); hax = axes('Parent', hf, 'NextPlot', 'add');
bar(hax, X, Y, 'EdgeColor', [0 0 1], 'FaceColor', [0 0 1]);
mm = mean(Y(~isnan(Y)));
line(seqvel, mm*ones(1,numel(seqvel)), 'Parent', hax, 'LineStyle', 'none', ...
    'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0]);
str{1} = [fdir, '__', evname];
str{2} = ['sess/ev mean=', num2str(mean(allvel)), ';se=', num2str(std(allvel)/sqrt(numel(allvel))), ';n=', num2str(numel(allvel)),...
           ';median=', num2str(median(allvel))];
str{3} = ['sequence mean=', num2str(mean(seqvel)), ';se=', num2str(std(seqvel)/sqrt(numel(seqvel))), ';n=', num2str(numel(seqvel)),...
           ';median=', num2str(median(seqvel))];  
str{4} = ['sequence zscores: mean=', num2str((mean(seqvel)-mean(allvel))/std(allvel)), '; median=', ...
    num2str((median(seqvel)-mean(allvel))/std(allvel))];
for (i = 1:numel(str))
    text('Parent', hax, 'Interpreter', 'none','String', str{i}, 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
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

function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end

function spikeX = findspikex(postime, posx, spiketime)
spikeX = zeros(size(spiketime));
[postime, iii] = sort(postime); posx = posx(iii);
lastpoint = 1;  %this only for saving time
for (j = 1:numel(spiketime))
    for (k = lastpoint:numel(postime)-1) %find corresponding time in position data
         if (postime(k) <= spiketime(j)) && (postime(k+1) > spiketime(j)) 
             spikeX(j) = posx(k); lastpoint = k; break; 
         end
    end
end
