function DM_Seq_PlotSequenceInterval
%%%based on sequencing results of a .seqdb (for the Tau mouse project)
%%%      plot the temporal/spatial interval distribution of the detected significant sequences

ok = 1;
evsigval = 1.645; %2; %-10; %1.645; %1; %1.645;
tracklength = 1150; %1150; %2300; %total length (2 trajectories) of of the track
%%%%%%% examples to plot
exptmp{1} = 'Tau21_111021TO_Track1_cw.evt';
%exptmp{2} = 'Tau13_110404TO_Track1_cw.evt';
exptmp{2} = 'Tau25_120118TO_Track1_ncw.evt'; %%%plot event rasters to show periodicity
%exptmp{1} = 'Tau21_111027ND1_Track1_cw.evt';
%exptmp{2} = 'Tau25_120124DOT_Track1_cw.evt';

hf = gcbf; pinfonow = getappdata(hf, 'pinfo'); datanow = getappdata(hf, 'data');
plotparm = getappdata(hf, 'plotparm');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); grpind = find(groupselection == 1); 
ngroup = numel(grpind); cellind = cell(1, ngroup); grpnames = cell(1, ngroup);
for (kk = 1:numel(grpind)) 
    cellind{kk} = datanow.grouplist.groupindex{grpind(kk)}; grpnames{kk} = datanow.grouplist.groupname{grpind(kk)}; 
end
if (plotparm.linkbehav == 0);
    behav = []; bhdata = [];
    disp(['--------> warning: no behavioral data linked; track spatial intervals not computed']); 
else
    behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
end
evtraj = {'selfTracktime'; 'opentime'; 'sleeptime'; 'selftrackspace'; 'selftrackdistance'; 'opendistance'; 'sleepdistance'}; ntraj = numel(evtraj);
           %%%selftrackdistance doesnot work - 
plotval = cell(ngroup, ntraj); %%%here is the interval
allval = cell(ngroup, ntraj); nval = zeros(ngroup, ntraj);%%%%contain real time/space points for auto-correlation plots
mm = NaN*ones(ngroup, ntraj); se = NaN*ones(ngroup, ntraj); 
md = NaN*ones(ngroup, ntraj);  nn = NaN*ones(ngroup, ntraj); 
ppgroup = NaN*ones(ngroup-1, ntraj); pptraj = NaN*ones(ngroup, ntraj-1);

expval = cell(numel(exptmp), ntraj);
for (i = 1:ngroup) %%%%self track running events: This should be done not only on self traj but also on the other traj
    tmpname = pinfonow.tmp.sessevname(cellind{i}); score = pinfonow.seq.evPosMatchNshufsig(cellind{i});
    evname = pinfonow.general.eventname(cellind{i}); evtype = pinfonow.parm.eventType(cellind{i}); tmpID = pinfonow.general.tmpID(cellind{i});
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
            aa = []; bb = []; cc = []; dd = [];
            jjj = find(sessyes & typeyes);
            for (ttt = 1:numel(jjj))
                [aaa, bbb, ccc, ddd] = findeventmeansequencevelocity(datanow, cellind{i}(k), jjj(ttt), plotparm.significancelevel, ...
                    pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, evname{k}{jjj(ttt)}, ttt, tracklength);  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
                aa = [aa aaa]; bb = [bb bbb]; cc = [cc ccc]; dd = [dd ddd];
            end
            if (~isempty(aa)) plotval{i,1} = [plotval{i,1} aa]; end
            if (~isempty(bb)) plotval{i,4} = [plotval{i,4} cc]; end%%%%event spatial intervals for the entire group
            if (~isempty(bb)) 
                 nval(i,1) = nval(i,1) + 1;
                 allval{i,1}{nval(i,1)} = bb; 
            end
            if (~isempty(dd))
                nval(i,4) = nval(i,4) + 1;
                allval{i,4}{nval(i,4)} =  dd; %%%%event locations for each template
            end
        end
        III = find(strcmp(exptmp, tmpID{k})); 
        for (ttt = 1:numel(III))
            expval{III(ttt),1} = bb; expval{III(ttt),4} = dd;
        end
    end
end
for (i = 1:ngroup) %%%%open/sleep sessions
    tmpname = pinfonow.tmp.sessevname(cellind{i}); score = pinfonow.seq.sessPosMatchNshufsig(cellind{i}); fdir = pinfonow.general.finaldir(cellind{i});
    sessname = pinfonow.general.sessionname(cellind{i}); sesstype = pinfonow.parm.sessType(cellind{i}); tmpID = pinfonow.general.tmpID(cellind{i});
    for (k = 1:numel(cellind{i}))
        [S, T] = strtok(tmpname{k}, '_'); S = lower(S);
        T = lower(T);
        sessnow = sessname{k}; sesstypenow = sesstype{k}; scorenow = score{k}; fdirnow = fdir{k};
        for (j = 1:numel(sessnow))
             [allD, allT] = finddistancetime(fdirnow, sessnow{j}, behav, bhdata);
             if (strcmp(sesstypenow{j}, 'open')) && (scorenow(j) > evsigval) 
                [aa, cc, ee, ff] = findsessmeansequencevelocity(datanow, cellind{i}(k), j, plotparm.significancelevel, ...
                    pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, sessname{k}{j}, allD, allT);  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
                if (~isempty(aa)) plotval{i,2} = [plotval{i,2} aa]; end
                if (~isempty(ee)) plotval{i,6} = [plotval{i,6} ee]; end
                if (~isempty(cc)) 
                    nval(i,2) = nval(i,2) + 1;
                    allval{i,2}{nval(i,2)} =  cc; 
                end
                if (~isempty(ff)) 
                    nval(i,6) = nval(i,6) + 1;
                    allval{i,6}{nval(i,6)} =  ff; 
                end
                III = find(strcmp(exptmp, tmpID{k})); %%%assign example-specific locations
                for (ttt = 1:numel(III)) 
                    expval{III(ttt),2} = cc; expval{III(ttt),6} = ff; %disp('******assign open now');
                end
             elseif (strcmp(sesstypenow{j}, 'linear')) && (scorenow(j) > evsigval) && (~isempty(strfind(lower(sessnow{j}), S)))
                [aa, cc, ee, ff] = findsessmeansequencevelocity(datanow, cellind{i}(k), j, plotparm.significancelevel, ...
                    pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, sessname{k}{j}, allD, allT);  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
                %if (~isempty(aa)) plotval{i,1} = [plotval{i,1} aa]; end
                if (~isempty(ee)) plotval{i,5} = [plotval{i,5} ee]; end
                %if (~isempty(cc)) 
                %    nval(i,1) = nval(i,1) + 1;
                %    allval{i,1}{nval(i,1)} =  cc; 
                %end
                if (~isempty(ff)) 
                    nval(i,5) = nval(i,5) + 1;
                    allval{i,5}{nval(i,5)} =  ff; 
                end
                III = find(strcmp(exptmp, tmpID{k})); %%%assign example-specific locations
                for (ttt = 1:numel(III)) 
                    expval{III(ttt),5} = ff; %disp('******assign linear now');
                end
             elseif (strcmp(sesstypenow{j}, 'sleep')) && (scorenow(j) > evsigval) 
                [aa, cc, ee, ff] = findsessmeansequencevelocity(datanow, cellind{i}(k), j, plotparm.significancelevel, ...
                    pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, sessname{k}{j}, allD, allT);  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
                if (~isempty(aa)) plotval{i,3} = [plotval{i,3} aa]; end
                if (~isempty(aa)) plotval{i,7} = [plotval{i,7} ee]; end
                if (~isempty(cc)) 
                    nval(i,3) = nval(i,3) + 1;
                    allval{i,3}{nval(i,3)} =  cc; 
                end
                if (~isempty(cc)) 
                    nval(i,7) = nval(i,7) + 1;
                    allval{i,7}{nval(i,7)} =  ff; 
                end
                III = find(strcmp(exptmp, tmpID{k})); %%%assign example-specific locations
                for (ttt = 1:numel(III)) 
                    expval{III(ttt),3} = cc; expval{III(ttt),7} = ff; %disp('******assign sleep now');
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

%%%%%First do event rasters and autocrr for example events
for (i = 1:numel(exptmp))
    for (j = 1:ntraj)
        spike = expval{i,j}; spike = sort(spike);
        binsize = 1;
        if (j >= 4) binsize = 25; end
        maxlag = 200*binsize;
        np = ceil(maxlag/binsize); lagbin = (-np:1:np); %%%%lagbin is a row vector
        X = lagbin*binsize; nnbin = numel(X);
        [crrok, norm, pp] = computecrrhere(spike, spike, min(spike), max(spike), lagbin, 'rate', binsize); %%%spikes here are column vectors
        ploteventraster(spike, crrok, pp, X, exptmp{i}, evtraj{j});
        %line(X, crrok, 'Parent', hax(i), 'LineWidth', 2, 'Color', rand(1,3));
    end
end

for (j = 1:ntraj)
%     %%%%do histgramps for sequence intervals
%     hf = figure('Name', ['SignificantSequencesIntervals-', evtraj{j}]); wid = 1/ngroup; gap = 0.05;
%     mmx = 0; nbin = 30;
%     for (i = 1:ngroup) mmx = max([mmx max(plotval{i,j})]); end
%     binsize = 2; 
%     %binsize = mmx/nbin; 
%     if (j>=4) binsize =20; end
%     X = 0:binsize:mmx; hax = [];
%     for (i = 1:ngroup)
%         hax(i) = axes('Parent', hf, 'NextPlot', 'add', 'Units', 'normalized', 'Position', [0.05 (ngroup-i)*wid+gap 0.95 wid-gap]);
%         if (~isempty(plotval{i,j}))
%            Y = histc(plotval{i,j}, X); 
%            bar(hax(i), X, Y, 'EdgeColor', [0 0 1], 'FaceColor', [0 0 1]);
%         end
%         str = grpnames{i};
%         text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.96]);
%         str = ['nn=[', num2str(nn(i,:)), ']; mean=[', num2str(mm(i,:)), ']']; 
%         text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.92]);
%         str = ['se=[', num2str(se(i,:)), ']; median=[', num2str(md(i,:)), ']'];
%         text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.88]);
%         str = ['sess/ev comparisons (t-test against the 1st) p=', num2str(pptraj(i,:))]; %, '%10.5e')];
%         text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.84]);
%         thnow = 0.84;
%         for (k = 1:ntraj) %%%%group comparisons
%             thnow = thnow - 0.04;
%             str = [evtraj{k}, ' group comparisons (t-test against the 1st) p=', num2str(ppgroup(:,k)')]; %, '%10.5e')];
%             text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 thnow]);
%         end
%     end

    %%%%do auto-correlations for all the time/space points
    hf = figure('Name', ['SignificantSequencesAutoCrr-', evtraj{j}]); wid = 1/ngroup; gap = 0.05;
    binsize = 1;
    if (j >= 4) binsize = 25; end
    maxlag = 200*binsize;
    np = ceil(maxlag/binsize); lagbin = (-np:1:np); %%%%lagbin is a row vector
    X = lagbin*binsize; nnbin = numel(X);
    for (i = 1:ngroup)
        hax(i) = axes('Parent', hf, 'NextPlot', 'add', 'Units', 'normalized', 'Position', [0.05 (ngroup-i)*wid+gap 0.95 wid-gap], 'YLim', [-0.05 0.3]);
        if (~isempty(allval{i,j}))
            allcrr = NaN*ones(nval(i,j), nnbin);
            for (kkt = 1:nval(i,j))
                 spike = allval{i,j}{kkt}; spike = sort(spike);
                 [crrok, norm] = computecrrhere(spike, spike, min(spike), max(spike), lagbin, 'rate', binsize); %%%spikes here are column vectors
                 allcrr(kkt,:) = crrok;
                 line(X, crrok, 'Parent', hax(i), 'LineWidth', 2, 'Color', rand(1,3));
            end
            crr = NaN*ones(1, nnbin); err = NaN*ones(1, nnbin); nnow = NaN*ones(1, nnbin); pp = NaN*ones(1, nnbin);
            for (kk = 1:nnbin)
                iii = find(~isnan(allcrr(:,kk))); crr(kk) = mean(allcrr(iii,kk)); [hh, pp(kk)] = ttest(allcrr(iii,kk), 0);
                nnow(kk) = numel(iii); 
                if (nnow(kk) > 0) err(kk) = std(allcrr(iii,kk))/sqrt(nnow(kk)); end
            end
            xlabel ('time lag (s)'); ylabel('Mean correlation'); 
            line(X, crr, 'Parent', hax(i), 'LineWidth', 2, 'Color', [0 0 0]);
            %itk = find((pp<0.05/(np+1))&(crr>0)); itt = setdiff(1:nnbin, itk);
            %line(X(itt), crr(itt), 'Parent', hax(i), 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', 'none',...
            %    'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
            %%Drawerrorupdown(X(itt), crr(itt), crr(itt)+err(itt), crr(itt)-err(itt), hax(i), [0 0 0]);
            %line(X(itk), crr(itk), 'Parent', hax(i), 'LineWidth', 2, 'Color', [1 0 0], 'LineStyle', 'none',...
            %    'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
            %%Drawerrorupdown(X(itk), crr(itk), crr(itk)+err(itk), crr(itk)-err(itk), hax(i), [0 0 0]);
        end
        str = grpnames{i};
        text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.96]);
        str = ['nn=[', num2str(nval(i,:)), ']']; %; mean=[', num2str(mm(i,:)), ']']; 
        text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.92]);
%         str = ['se=[', num2str(se(i,:)), ']; median=[', num2str(md(i,:)), ']'];
%         text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.88]);
%         str = ['sess/ev comparisons (t-test against the 1st) p=', num2str(pptraj(i,:))]; %, '%10.5e')];
%         text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.84]);
%         thnow = 0.84;
%         for (k = 1:ntraj) %%%%group comparisons
%             thnow = thnow - 0.04;
%             str = [evtraj{k}, ' group comparisons (t-test against the 1st) p=', num2str(ppgroup(:,k)')]; %, '%10.5e')];
%             text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 thnow]);
%         end
    end
end
disp('*******************************');

function [timeint, alltime, distint, alldist] = findsessmeansequencevelocity(datanow, cellind, sessind, siglevel,...
    pinfo, data, behav, bhdata, finaldirnow, sessnamenow, allD, allT)
seqnow = datanow.data.sessseq{cellind}{sessind};
seqstart = cell2mat(seqnow.seqstart); seqend = cell2mat(seqnow.seqend); 
posmatchprob = cell2mat(seqnow.posmatchprob); %matchscore = seqnow.matchscore;  negmatchprob = seqnow.negmatchprob;
iii = find( posmatchprob < siglevel ); 
posmatchprob = posmatchprob(iii); seqstart = seqstart(iii); seqend = seqend(iii); %%%only the significant sequences
[newstart, newend, ~, ~] = resolveoverlap(posmatchprob, seqstart, seqend); %%resolve the overlaps
 %%%this contains all the center times of the significant sequences
centertime = sort((newstart + newend)/2); timeint = diff(centertime); alltime = centertime;
alldist = finddistance(allD, allT, centertime); distint = diff(alldist); 

function [timeint, alltime, spaceint, allspace] = findeventmeansequencevelocity(datanow, cellind, evind, siglevel,...
    pinfo, data, behav, bhdata, finaldirnow, evnamenow, ntraj, tracklength)
timeint = []; spaceint = []; alltime = []; allspace = [];
seqnow = datanow.data.evseq{cellind}{evind};
seqstart = cell2mat(seqnow.seqstart); seqend = cell2mat(seqnow.seqend); 
posmatchprob = cell2mat(seqnow.posmatchprob); %matchscore = seqnow.matchscore;  negmatchprob = seqnow.negmatchprob;
iii = find( posmatchprob < siglevel ); 
posmatchprob = posmatchprob(iii); seqstart = seqstart(iii); seqend = seqend(iii); %%%only the significant sequences
[newstart, newend, ~, ~] = resolveoverlap(posmatchprob, seqstart, seqend); %%resolve the overlaps
centertime = sort((newstart+newend)/2); alltime = centertime;
 %%%this contains all the center times of the significant sequences
timeint =[timeint diff(sort(centertime))];
 
%%%%locate event position data
evTime = data.events.eventimes{cellind}{evind};
%for (i = 1:numel(evTime.start))
%    iii = find( (centertime>=evTime.start(i)) & (centertime<=evTime.ent(i)) );
%    timeint =[timeint diff(sort(centertime(iii)))];
%end
if (~isempty(timeint)) && (~isempty(behav))
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
      accudis = (ntraj-1) * tracklength/2; %%%starting position of the current trajectory
      for (tj = 1:nlap)
        lappostime = bhdata.event.LapAllPostimestamp{posid}{evid}{tj}; 
        lapx = bhdata.event.LapAllX{posid}{evid}{tj}; 
        minT = min(lappostime); maxT = max(lappostime);
        iii = find( (centertime>=minT) & (centertime<=maxT) );
        lapseqcenter = NaN*ones(1,numel(iii));
        for (tk = 1:numel(iii))
            tnow = centertime(iii(tk));
            lapseqcenter(tk) = accudis + findspikex(lappostime, lapx, tnow);
        end
        spaceint = [spaceint diff(sort(lapseqcenter))];
        allspace = [allspace lapseqcenter];
        accudis = accudis + tracklength; %2*max(lapx); %add whole track distance to gaps between laps
      end
    end
end

function [alldistance, alltime] = finddistancetime(fdirnow, sessName, behav, bhdata)
posid = find( strcmp(behav.general.finaldir, fdirnow) & strcmp(behav.general.sessname, sessName) );
if numel(posid)~=1
    disp(['-------------> warning: no or more than 1 positon files match the session: ', fdirnow, '___', drddName]);
else
    alltime = bhdata.pos.postimestamp{posid}*behav.parm.timeunit(posid); 
    alldistance = bhdata.sess.AllTrajDistance{posid};
end
function spikex = finddistance(allxnow, alltimenow, spikenow)
spikex = NaN*ones(size(spikenow)); spikenow = sort(spikenow);
lastpoint = 1;  %this only for saving time
for (j = 1:numel(spikenow))
         for (k = lastpoint:numel(alltimenow)-1) %find corresponding time in position data
             if (alltimenow(k) <= spikenow(j)) && (alltimenow(k+1) > spikenow(j)) 
                 spikex(j) = allxnow(k); lastpoint = k; break; 
             end
         end
end

function plotveldistribution(allvel, seqvel, fdir, evname)
nbin = 40; maxV = max(allvel); binsize = maxV/nbin;
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

function newcrr = smoothcrr(crr, smoothbin)
newcrr = zeros(size(crr)); %NaN*ones(size(crr)); 
halfbin = floor(smoothbin/2);
for (i = 1:numel(crr))
    kkk = [i-halfbin:i+halfbin]; kkk = kkk( (kkk>=1) & (kkk<=numel(crr)) );
    crrnow = crr(kkk); newcrr(i) = mean(crrnow(~isnan(crrnow)));
end

function [crr, norm, pp] = computecrrhere(spike1, spike2, sT, eT, lagbin, crrmode, binsize)
crr = zeros(size(lagbin)); pp = ones(size(lagbin)); 
norm = 0;
if (~isempty(spike1)) && (~isempty(spike2))
    %iii = find(eT-sT>=binsize); sT = sT(iii); eT = eT(iii); %%%filter out short events
    if (strncmpi(crrmode, 'rate', 2))
        [crr, pp] = computeratecrr(spike1, spike2, sT, eT, lagbin, binsize);
    elseif (strncmpi(crrmode, 'count', 2))
        crr = computecountcrr(spike1, spike2, sT, eT, lagbin, binsize);
    elseif (strncmpi(crrmode, 'normcount', 2))
        [crr, norm] = computenormcountcrr(spike1, spike2, sT, eT, lagbin, binsize);
    end
end

function [ccc,pp] = computeratecrr(t1, t2, sT, eT, lagbin, binsize) %%%t1, t2 are row vectors
nev = numel(sT); timebin = cell(1,nev); nc = numel(lagbin); ccc = zeros(size(lagbin)); pp = ones(size(lagbin)); 
% %%%%%%%%%%need to bin spikes first --- It is specific here, assume t1 = t2!!!! Only for autocrr
% nbin = ceil( (max(t1)-min(t1))/binsize );
% timebin = min(t1) + ( (1:nbin) - 1)*binsize; 
% cnt = histc(t1, timebin);
% for (i = 1:nc)
%     startbin1 = max([1 1-lagbin(i)]);
%     endbin1 = min([nbin nbin-lagbin(i)]);
%     startbin2 = max([1 1+lagbin(i)]);
%     endbin2 = min([nbin nbin+lagbin(i)]);
%     cn1 = cnt(startbin1:endbin1); cn2 = cnt(startbin2:endbin2);
%     [R,P] = corrcoef(cn1, cn2); if (numel(R)>1) ccc(i) = R(1,2); pp(i) = P(1,2); end
% end

%%%%%%%%%%%%%decide not to use the event-based crr below since events already assigned non-overlapping time points
cnt1 = cell(1,nev); cnt2 = cell(1,nev);
for (i = 1:nev)
     nbin = ceil( (eT(i)-sT(i))/binsize );
     timebin{i} = sT(i) + ( (1:nbin) - 1 )*binsize;
end
%tttt = cputime;
for (i = 1:nev)
    t1now = t1( (t1>=sT(i)) & (t1<=eT(i)) ); t2now = t2( (t2>=sT(i)) & (t2<=eT(i)) ); %%%filter here for saving hisc time?
    nbinnow = numel(timebin{i}) - 1; 
    cnt1{i} = zeros(1, nbinnow); cnt2{i} = zeros(nbinnow, 1);
    if (~isempty(t1now))
       cnt = histc(t1now, timebin{i}); cnt1{i} = cnt(1:nbinnow); %output row vector (histc output same row or column vector as t1now)
    end
    if (~isempty(t2now))
       cnt = histc(t2now, timebin{i}); cnt2{i} = cnt(1:nbinnow)'; %output column vector
    end
end
%tttt1 = cputime;
%disp(strcat('*****multiple ep binning time:*****', num2str(tttt1-tttt)));
%move 2nd t2 firing rate time series accroding to time lag: if lagbin <0, move t2 to the right against t1
for (i= 1:nc)
    %tttt = cputime;
    nbin = 0; sum1 = 0; sum2 = 0; sumsq1 = 0; sumsq2 = 0; sumsq12 = 0; %catecated count vector cross all episodes for corrcoeff
    for (j = 1:nev)
        cn1 = cnt1{j}; cn2 = cnt2{j}; %%%%using matrix instead of cell saves a lot of time
        nbinnow = numel(timebin{j}) - 1; 
        startbin1 = max([1 1-lagbin(i)]);
        endbin1 = min([nbinnow nbinnow-lagbin(i)]);
        startbin2 = max([1 1+lagbin(i)]);
        endbin2 = min([nbinnow nbinnow+lagbin(i)]);
        if ( (startbin1<endbin1) & (startbin2<endbin2) ) %if both series have space for the lag in this episode
            nbin = nbin + endbin1 - startbin1 + 1;
            sum1 = sum1 + sum(cn1(startbin1:endbin1)); sum2 = sum2 + sum(cn2(startbin2:endbin2));
            sumsq1 = sumsq1 + sum(cn1(startbin1:endbin1).^2); sumsq2 = sumsq2 + sum(cn2(startbin2:endbin2).^2);
            sumsq12 = sumsq12 + cn1(startbin1:endbin1) * cn2(startbin2:endbin2); %%%row vector times column vector
        end
    end
    %tttt1 = cputime;
    %disp(strcat('*****crr single point computing time:*****', num2str(tttt1-tttt)));
    if (nbin ~= 0)
       co1 = sumsq1-sum1^2/nbin; co2 = sumsq2-sum2^2/nbin;
       if ( (co1 == 0) | (co2 == 0) )
           ccc(i) = 0;
       else
           ccc(i) = (sumsq12 - sum1*sum2/nbin) / sqrt(co1*co2);
       end
    end
end

function count = computecountcrr(spike1, spike2, sT, eT, lagbin, bin)
count = zeros(size(lagbin));   %initial assignment for spike coun
np = lagbin(numel(lagbin)); %%%number of left and right lag points
t1 = []; t2 = [];
for (i = 1:numel(sT))
    t1 = [t1 spike1( (spike1>=sT(i)) & (spike1<=eT(i)) )];
    t2 = [t2 spike2( (spike2>=sT(i)) & (spike2<=eT(i)) )]; %%%do not filter the second train to avoid the edge effect?
                                                            %%%better do it to reflect the true event windows;
                                                            %%%   be aware of the edge effect for short event files
end
for (i = 1:numel(t1))
    % for each spike in first train, count spike number in second train that within a time window defined by bin
%     for (k = -np:1:np)
%         minn = t1(i) + k*bin - bin/2;  %defining the window
%         maxx = t1(i) + k*bin + bin/2;
%         count(k+np+1) = count(k+np+1) + numel(find( (t2>=minn) & (t2<maxx)));
%     end
    k = -np:1:np; wind = t1(i) + k*bin - bin/2;
    count(k+np+1) = count(k+np+1) + histc(t2, wind);
end

function [count, mm] = computenormcountcrr(spike1, spike2, sT, eT, lagbin, bin)
%%%%normalized count: assuming spike1 and spike2 are independent poisson, at each time lag bin, the mean expected spike count is
%%%%       mm = numel(spike1)*numel(spike2)*bin/totallengthofspiketrains; 
%%%% and the std (for Poisson) is ss = sqrt(mm)
%%%% normalized as (count-mm)/ss
%%%% lagbin = (-np:1:np); 
count = zeros(size(lagbin)); mm = 0;   %initial assignment for spike coun = row vector
np = lagbin(numel(lagbin)); %%%number of left and right lag points
t1 = []; t2 = []; tleng = 0;
for (i = 1:numel(sT))
    tleng = tleng + (eT(i)-sT(i));
    t1 = [t1 spike1( (spike1>=sT(i)) & (spike1<=eT(i)) )];
    t2 = [t2 spike2( (spike2>=sT(i)) & (spike2<=eT(i)) )]; %%%%need to filter to compute the norm factor
end
%tttt1 = cputime;
if (~isempty(t2))
for (i = 1:numel(t1)) % for each spike in first train, count spike number in second train that within a time window defined by bin
%     for (k = -np:1:np)
%         minn = t1(i) + k*bin - bin/2;  %defining the window
%         maxx = t1(i) + k*bin + bin/2;
%         count(k+np+1) = count(k+np+1) + numel(find( (t2>=minn) & (t2<maxx)));
%     end
      k = -np:1:np; wind = t1(i) + k*bin - bin/2;
      count(k+np+1) = count(k+np+1) + histc(t2, wind);
end
%tttt2 = cputime;
%disp(strcat('*****crr single point counting time:*****', num2str(tttt2-tttt1)));
%%%%%normalization
mm = numel(t1)*numel(t2)*bin/tleng;
if (mm>0) count = (count-mm)/sqrt(mm); end
end

function ploteventraster(spike, crrok, pp, X, tmpname, trajname)
hfig = figure('Name', strcat(tmpname, '__', trajname), 'NextPlot', 'add', 'Color', [1 1 1]);
posvec = [0.05 0.05 0.90 0.90]; nticktrace = numel(spike); np = (numel(X)-1)/2;
for (i=1:nticktrace) %first tick trace is on the top.
         %position for each episode
         mtickposvec{i} = [posvec(1) posvec(2)+posvec(4)/3+(nticktrace-i)*posvec(4)*2/3/nticktrace posvec(3) posvec(4)*2/3/nticktrace];    % multiple tick traces
         %axes handle for each episode
         hmtickaxes(i) = axes('Parent', hfig, 'Units', 'normalized', 'Position', mtickposvec{i}, 'Visible', 'off', 'NextPlot', 'add', 'XLim', [min(X) max(X)], 'YLim', [0 1.2]);
         %xlim ([min(X) max(X)]); %plot x-axis as linearized track
         %ylim ([0 1.2]);
end
%%% plot rastors for each episode and output spike timing for each episode
for (i=1:nticktrace) %for each episode
        spikenow = spike - spike(i);
        for (k=1:numel(spikenow))
              line([spikenow(k); spikenow(k)], [0; 1], 'Parent', hmtickaxes(i), 'LineWidth', 4, 'Color', [0 0 0]) % plot spikes as ticks
        end
end
%%plot rate (spike count in space grid) trace   %%%only one rate trace!!
rateposvec = [posvec(1) posvec(2) posvec(3) posvec(4)/4]; 
hrateaxes = axes('Parent', hfig, 'Units', 'normalized', 'Position', rateposvec, 'FontSize', 8, 'XLim', [min(X) max(X)], 'YLim', [-0.05 0.3]);
%xlim ([min(X) max(X)]); %plot x-axis as linearized track 
line(X, crrok, 'Parent', hrateaxes, 'LineWidth', 2, 'Color', [0 0 0]);
%itk = find((pp<0.05/(np+1))&(crrok>0)); itt = setdiff(1:numel(X), itk);
%line(X(itt), crrok(itt), 'Parent', hrateaxes, 'LineWidth', 2, 'Color', [0 0 0]); %, 'LineStyle', 'none',...
%               % 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
%line(X(itk), crrok(itk), 'Parent', hrateaxes, 'LineWidth', 2, 'Color', [1 0 0]); %, 'LineStyle', 'none',...
%               % 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);


