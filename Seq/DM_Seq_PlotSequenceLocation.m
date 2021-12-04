function DM_Seq_PlotSequenceLocation
disp(['----------> sequence location display for group(s) not available; try to combine plots from individual events.'])
%%% Based on sequencing/decoding results of a .seqdb, characterize various quantities of sequence locations
%%% !!!!!!This is project specific - not done at this time!!!!!!!!
%%%      test whether sequences on open fields are preferentially located 
%%%           test whther the sequence start -> end vectors of multiple sequences are aligned.
%%% Do this for SelfTrack, Open, Sleep categories 

% ok = 1;
% hf = gcbf; pinfonow = getappdata(hf, 'pinfo'); datanow = getappdata(hf, 'data');
% plotparm = getappdata(hf, 'plotparm');
% hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); grpind = find(groupselection == 1); 
% ngroup = numel(grpind); cellind = cell(1, ngroup); grpnames = cell(1, ngroup);
% for (kk = 1:numel(grpind)) 
%     cellind{kk} = datanow.grouplist.groupindex{grpind(kk)}; grpnames{kk} = datanow.grouplist.groupname{grpind(kk)}; 
% end
% if (plotparm.linkbehav == 0)
%     behav = []; bhdata = [];
%     disp(['--------> warning: no behavioral data linked; aborted']); ok = 0;
% else
%     behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
% end
% 
% if ok
% evsigval = 1.645;
% evtraj = {'selfTrack'; 'open'; 'sleep'}; ntraj = numel(evtraj);
% plotval = cell(ngroup, ntraj); mm = NaN*ones(ngroup, ntraj); se = NaN*ones(ngroup, ntraj); 
% md = NaN*ones(ngroup, ntraj);  nn = NaN*ones(ngroup, ntraj); 
% ppgroup = NaN*ones(ngroup-1, ntraj); pptraj = NaN*ones(ngroup, ntraj-1);
% % for (i = 1:ngroup) %%%%self track running events
% %     tmplength = pinfonow.tmp.length(cellind{i});
% %     tmpname = pinfonow.tmp.sessevname(cellind{i}); score = pinfonow.seq.evPosMatchNshufsig(cellind{i});
% %     evname = pinfonow.general.eventname(cellind{i}); evtype = pinfonow.parm.eventType(cellind{i});
% %     for (k = 1:numel(cellind{i}))
% %         [S, T] = strtok(tmpname{k}, '_'); S = lower(S); T = lower(T);
% %         evnow = evname{k}; evtypenow = evtype{k}; scorenow = score{k};
% %         sessyes = zeros(1, numel(evnow)); trajyes = zeros(1, numel(evnow)); 
% %         typeyes = zeros(1, numel(evnow)); scoreyes = zeros(1, numel(evnow));
% %         for (j = 1:numel(evnow))
% %                 if (~isempty(strfind(lower(evnow{j}), S))) sessyes(j) = 1; end
% %                 if (~isempty(strfind(lower(evnow{j}), T))) trajyes(j) = 1; end
% %                 if strcmp(evtypenow{j}, 'run') typeyes(j) = 1; end
% %                 if (scorenow(j) > evsigval) scoreyes(j) = 1; end
% %         end
% %         ii = find(sessyes & trajyes & typeyes & scoreyes); 
% %         if (numel(ii) == 1) 
% %             aa = findeventmeansequencevelocity(datanow, cellind{i}(k), ii, plotparm.significancelevel, ...
% %                     pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, evname{k}{ii});  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
% %             if (~isempty(aa)) plotval{i,1} = [plotval{i,1} aa]; end
% %         end
% %     end
% % end
% for (i = 1:ngroup) %%%%track/open/sleep sessions
%     tmplength = pinfonow.tmp.length(cellind{i});
%     tmpname = pinfonow.tmp.sessevname(cellind{i}); score = pinfonow.seq.sessPosMatchNshufsig(cellind{i});
%     sessname = pinfonow.general.sessionname(cellind{i}); sesstype = pinfonow.parm.sessType(cellind{i});
%     for (k = 1:numel(cellind{i}))
%         sessnow = sessname{k}; sesstypenow = sesstype{k}; scorenow = score{k};
%         for (j = 1:numel(sessnow))
%              if (strcmp(sesstypenow{j}, 'open')) && (scorenow(j) > evsigval) 
%                 aa = findsessmeansequencevelocity(datanow, cellind{i}(k), j, plotparm.significancelevel, ...
%                     pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, sessname{k}{j});  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
%                 if (~isempty(aa))
%                     plotval{i,2} = [plotval{i,2} aa];
%                 end
%              elseif (strcmp(sesstypenow{j}, 'sleep')) && (scorenow(j) > evsigval) 
%                 aa = findsessmeansequencevelocity(datanow, cellind{i}(k), j, plotparm.significancelevel, ...
%                     pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, sessname{k}{j});  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
%                 if (~isempty(aa))
%                     plotval{i,3} = [plotval{i,3} aa];
%                 end
%              elseif (strcmp(sesstypenow{j}, 'linear')) && (scorenow(j) > evsigval) 
%                 aa = findsessmeansequencevelocity(datanow, cellind{i}(k), j, plotparm.significancelevel, ...
%                     pinfonow, datanow, behav, bhdata, pinfonow.general.finaldir{cellind{i}(k)}, sessname{k}{j});  %%%get the position variation = val(1, nev) but only do it for 'run' type events 
%                 if (~isempty(aa))
%                     plotval{i,1} = [plotval{i,1} aa];
%                 end
%              end
%         end
%     end
% end
% for (i = 1:ngroup) %%%%basic statistics
%         for (j = 1:numel(evtraj))
%             aa = plotval{i,j}; aa = aa(~isnan(aa)); plotval{i,j} = aa;
%             nn(i,j) = numel(aa);
%             if (nn(i,j) > 1)
%                 mm(i,j) = mean(aa); se(i,j) = std(aa)/sqrt(nn(i,j)); md(i,j) = median(aa);
%             end
%         end
% end
% for (i = 1:ngroup-1) %%%%group comparasons: all compare with the first one
%         for (j = 1:ntraj)
%             if (~isempty(plotval{1,j})) && (~isempty(plotval{i+1,j}))
%                [~, ppgroup(i,j)] = ttest2(plotval{1,j}, plotval{i+1,j});
%             end
%         end
% end
% for (j = 1:ntraj-1) %%%%traj comparasons: all compare with the first one
%         for (i = 1:ngroup)
%             if (~isempty(plotval{i,1})) && (~isempty(plotval{i,j+1}))
%                [~, pptraj(i,j)] = ttest2(plotval{i,1}, plotval{i,j+1});
%             end
%         end
% end
% 
% 
%     
%     hf = figure('Name', 'SignificantSequencesLocationVectorReyleighP-'); wid = 1/ngroup; gap = 0.05;
%     for (i = 1:ngroup)
%         hax(i) = axes('Parent', hf, 'NextPlot', 'add', 'Units', 'normalized', 'Position', [(i-1)*wid+gap 0.05 wid-gap 0.95]);
%         for (j = 1:ntraj)
%             YY = plotval{i,j}; XX = j*ones(1, nn(i,j)); 
%             bar(hax(i), j, mm(i,j), 'k');
%             line(XX, YY, 'Parent', hax(i), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30);
%             %line([j-0.2 j+0.2], [mm(i,j) mm(i,j)], 'Parent', hax(i), 'LineWidth', 2); %%%mean value line 
%         end
%         line([0 ntraj+1], [0.05 0.05], 'Parent', hax(i)); %%%significance line
%         str = [grpnames{i}, ': self, sameTraj, sameSess, allDiff'];
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
%     
% %     mmx = 0; nbin = 30;
% %     for (i = 1:ngroup) mmx = max([mmx max(plotval{i,j})]); end
% %     binsize = 0.02; %binsize = mmx/nbin; 
% %     if (j==4) binsize =10; end
% %     X = 0:binsize:mmx; hax = [];
% %     for (i = 1:ngroup)
% %         hax(i) = axes('Parent', hf, 'NextPlot', 'add', 'Units', 'normalized', 'Position', [0.05 (ngroup-i)*wid+gap 0.95 wid-gap]);
% %         if (~isempty(plotval{i,j}))
% %            Y = histc(plotval{i,j}, X); 
% %            bar(hax(i), X, Y, 'EdgeColor', [0 0 1], 'FaceColor', [0 0 1]);
% %         end
% %         str = grpnames{i};
% %         text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.96]);
% %         str = ['nn=[', num2str(nn(i,:)), ']; mean=[', num2str(mm(i,:)), ']']; 
% %         text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.92]);
% %         str = ['se=[', num2str(se(i,:)), ']; median=[', num2str(md(i,:)), ']'];
% %         text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.88]);
% %         str = ['sess/ev comparisons (t-test against the 1st) p=', num2str(pptraj(i,:))]; %, '%10.5e')];
% %         text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 0.84]);
% %         thnow = 0.84;
% %         for (k = 1:ntraj) %%%%group comparisons
% %             thnow = thnow - 0.04;
% %             str = [evtraj{k}, ' group comparisons (t-test against the 1st) p=', num2str(ppgroup(:,k)')]; %, '%10.5e')];
% %             text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 thnow]);
% %         end
% %     end
% 
% end
disp('*******************************');

function val = findsessmeansequencevelocity(datanow, cellind, sessind, siglevel,...
    pinfo, data, behav, bhdata, finaldirnow, sessnamenow)
val = NaN;
seqnow = datanow.data.sessseq{cellind}{sessind};
seqstart = cell2mat(seqnow.seqstart); seqend = cell2mat(seqnow.seqend); 
posmatchprob = cell2mat(seqnow.posmatchprob); %matchscore = seqnow.matchscore;  negmatchprob = seqnow.negmatchprob;
iii = find( posmatchprob < siglevel ); 
posmatchprob = posmatchprob(iii); seqstart = seqstart(iii); seqend = seqend(iii); %%%only the significant sequences
[newstart, newend, ~, ~] = resolveoverlap(posmatchprob, seqstart, seqend); %%resolve the overlaps
 %%%this contains all the start/end times of the significant sequences
%%%%locate session position data
posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, sessnamenow) );
if numel(posid == 1)
    i = posid;
    postimestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i);
    Pmarker = behav.parm.sessPmarker{i}; allposmarker = behav.general.posMarker{i}; 
    posXX = []; posYY = []; ik = find(strcmp(allposmarker, Pmarker)); 
    if (numel(ik) == 1) posXX = bhdata.pos.XX{i}{ik}; posYY = bhdata.pos.YY{i}{ik}; end
   
    %%%%%%%%%sequence location
    nseq = numel(newstart); seqvecX = NaN*ones(1, nseq); seqvecY = NaN*ones(1, nseq);
    for (i = 1:nseq)
        iii = find( (postimestamp>=newstart(i)) & (postimestamp<=newend(i)) );
        seqvecX(i) = posXX(iii(numel(iii))) - posXX(iii(1)); seqvecY(i) = posYY(iii(numel(iii))) - posYY(iii(1));
    end
    ang = cart2pol(seqvecX, seqvecY); val = circ_rtest(ang); 
end


function [timeint, spaceint] = findeventmeansequencevelocity(datanow, cellind, evind, siglevel,...
    pinfo, data, behav, bhdata, finaldirnow, evnamenow)
timeint = []; spaceint = [];
seqnow = datanow.data.evseq{cellind}{evind};
seqstart = cell2mat(seqnow.seqstart); seqend = cell2mat(seqnow.seqend); 
posmatchprob = cell2mat(seqnow.posmatchprob); %matchscore = seqnow.matchscore;  negmatchprob = seqnow.negmatchprob;
iii = find( posmatchprob < siglevel ); 
posmatchprob = posmatchprob(iii); seqstart = seqstart(iii); seqend = seqend(iii); %%%only the significant sequences
[newstart, newend, ~, ~] = resolveoverlap(posmatchprob, seqstart, seqend); %%resolve the overlaps
centertime = (newstart+newend)/2;
 %%%this contains all the center times of the significant sequences
%%%%locate event position data
evTime = data.events.eventimes{cellind}{evind};
for (i = 1:numel(evTime.start))
    iii = find( (centertime>=evTime.start(i)) & (centertime<=evTime.ent(i)) );
    timeint =[timeint newend(iii)-newstart(iii)];
end
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
      nlap = numel(evTime.start); 
      for (tj = 1:nlap)
        lappostime = bhdata.event.LapAllPostimestamp{posid}{evid}{tj}; 
        lapx = bhdata.event.LapAllX{posid}{evid}{tj};
        minT = min(lappostime); maxT = max(lappostime);
        iii = find( (centertime>=minT) & (centertime<=maxT) );
        laplenth = NaN*ones(1,numel(iii));
        for (tk = 1:numel(iii))
            tsnow = newstart(iii(tk)); tenow = newend(iii(tk));
            ssnow = findspikex(lappostime, lapx, tsnow);
            senow = findspikex(lappostime, lapx, tenow);
            laplenth(tk) = senow-ssnow;
        end
        spaceint = [spaceint laplenth];
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
