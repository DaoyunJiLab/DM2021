function DM_Seq_PlotSequencePositionStability
%%%based on sequencing results of a .seqdb (for the Tau mouse project)
%%%      plot how sequence center positions shift lap-by-lap on a linear track
%%%          shift defined as the mean of    given any sequence in a lsp, find the nearest distance from the center of any sequences in another lap
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

%%%%event sequencing results
if ok
evtraj = {'self'; 'samTraj'; 'sameSess'; 'allDiff'}; ntraj = numel(evtraj);
plotval = cell(ngroup, ntraj); mm = NaN*ones(ngroup, ntraj); se = NaN*ones(ngroup, ntraj); 
md = NaN*ones(ngroup, ntraj);  nn = NaN*ones(ngroup, ntraj); 
ppgroup = NaN*ones(ngroup-1, ntraj); pptraj = NaN*ones(ngroup, ntraj-1);
for (i = 1:ngroup)
    for (j = 1:ntraj) 
        plotval{i,j} = NaN*ones(1,numel(cellind{i}));   %%%initial assignments
    end
    tmpname = pinfonow.tmp.sessevname(cellind{i});
    evname = pinfonow.general.eventname(cellind{i}); evtype = pinfonow.parm.eventType(cellind{i});
    for (k = 1:numel(cellind{i}))
        val = NaN*ones(1, numel(evname{k}));
        for (ii = 1:numel(evname{k}))
            if (strcmp(evtype{k}{ii}, 'run'))
                val(ii) = findmeansequencevelocity(datanow, cellind{i}(k), ii, plotparm.significancelevel, ...
                    pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, evname{k}{ii});  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
            end
        end
        [S, T] = strtok(tmpname{k}, '_'); S = lower(S); T = lower(T);
        evnow = evname{k}; evtypenow = evtype{k};
        sessyes = zeros(1, numel(evnow)); trajyes = zeros(1, numel(evnow)); typeyes = zeros(1, numel(evnow));
        for (j = 1:numel(evnow))
                if (~isempty(strfind(lower(evnow{j}), S))) sessyes(j) = 1; end
                if (~isempty(strfind(lower(evnow{j}), T))) trajyes(j) = 1; end
                if strcmp(evtypenow{j}, 'run') typeyes(j) = 1; end
        end
        ii = find(sessyes & trajyes & typeyes); if (numel(ii) == 1) plotval{i,1}(k) = val(ii); end %%%self
        ii = find( (~sessyes) & (trajyes) & typeyes); if (numel(ii) == 1) plotval{i,2}(k) = val(ii); end %%%same traj
        ii = find((sessyes) & (~trajyes) & typeyes); if (numel(ii) == 1) plotval{i,3}(k) = val(ii); end %%%self sess
        ii = find((~sessyes) & (~trajyes) & typeyes); if (numel(ii) == 1) plotval{i,4}(k) = val(ii); end %%%all diff
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
hf = figure('Name', ['SignificantSequences - lap location shift']); wid = 1/ngroup; gap = 0.05;
for (i = 1:ngroup)
        hax(i) = axes('Parent', hf, 'NextPlot', 'add', 'Units', 'normalized', 'Position', [(i-1)*wid+gap 0.05 wid-gap 0.95]);
        for (j = 1:ntraj)
            YY = plotval{i,j}; XX = j*ones(1, nn(i,j)); 
            line(XX, YY, 'Parent', hax(i), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30);
            line([j-0.2 j+0.2], [mm(i,j) mm(i,j)], 'Parent', hax(i), 'LineWidth', 2); %%%mean value line 
        end
        str = [grpnames{i}, ': self, sameTraj, sameSess, allDiff'];
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

function val = findmeansequencevelocity(datanow, cellind, evind, siglevel, pinfo, data, behav, bhdata, finaldirnow, evnamenow)
val = NaN;
seqnow = datanow.data.evseq{cellind}{evind};
seqstart = cell2mat(seqnow.seqstart); seqend = cell2mat(seqnow.seqend); 
posmatchprob = cell2mat(seqnow.posmatchprob); %matchscore = seqnow.matchscore;  negmatchprob = seqnow.negmatchprob;
iii = find( posmatchprob < siglevel ); 
posmatchprob = posmatchprob(iii); seqstart = seqstart(iii); seqend = seqend(iii); %%%only the significant sequences
[newstart, newend, ~, ~] = resolveoverlap(posmatchprob, seqstart, seqend); %%resolve the overlaps
seqcenter = (newstart + newend)/2; %%%this contains all the center times of the significant sequences
%%%%locate event position data
evTime = data.events.eventimes{cellind}{evind};
evSess = identifysession(evTime, pinfo.general.sessionname{cellind}, pinfo.general.sessionstartT{cellind}, pinfo.general.sessionendT{cellind});
posid = []; evid = []; lappostime =[]; lapx = [];
if (~isempty(evSess))
    posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess) );
end
if numel(posid == 1)
     if (isfield(behav.general, 'eventname'))
         evid = find(strcmp(behav.general.eventname{posid}, evnamenow));
     else
         evid = find(strcmp(behav.behavior.eventname{posid}, evnamenow));
     end
end
if (numel(posid)~=1)|(numel(evid)~=1)
      disp(['-------------> sequence lap positions not computed: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evnamenow]);
else
      %%%%get lap-by-lap time/positions
      nlap = numel(evTime.start); lappostime = cell(1, nlap); lapx = cell(1, nlap);
      for (tj = 1:nlap)
            lappostime{tj} = bhdata.event.LapAllPostimestamp{posid}{evid}{tj}; 
            lapx{tj} = bhdata.event.LapAllX{posid}{evid}{tj};
      end
end
if ~isempty(lappostime)
    %%%%%identify sequence laps and locations
    nlap = numel(lapx); lapseqcenter = cell(1, nlap);
    for (i = 1:nlap)
        minT = min(lappostime{i}); maxT = max(lappostime{i});
        iii = find( (seqcenter>=minT) & (seqcenter<=maxT) );
        lapseqcenter{i} = NaN*ones(1,numel(iii));
        for (tj = 1:numel(iii))
            tnow = seqcenter(iii(tj));
            lapseqcenter{i}(tj) = findspikex(lappostime{i}, lapx{i}, tnow);
        end
    end
    %%%%compute mean location shift
    locshift = [];
    for (i = 1:nlap)
        for (j = i+1:nlap)
            shiftnow = findlocsjift(lapseqcenter{i}, lapseqcenter{j});
            if ~isempty(shiftnow)
                locshift = [locshift shiftnow];
            end
        end
    end
    if (numel(locshift)>5) %%%get rid of outliers if too small number of laps
        val = median(locshift); %disp(locshift);
    end
end

function shiftnow = findlocsjift(lap1, lap2)
shiftnow = [];
if (~isempty(lap1)) && (~isempty(lap2))
   if (numel(lap1) > numel(lap2))
       dd = lap2; lap2 = lap1; lap1 = dd; %%%switch of lap1 has more seqs
   end
   val = zeros(1, numel(lap1));
   for (i = 1:numel(lap1))
       val(i) = min(abs(lap1(i)-lap2));
   end 
   shiftnow = mean(val);
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
