function DM_Seq_PlotSequenceResults
%%%plot sequencing results of a .seqdb (for the Tau mouse project)
%1. Plot zscores/p-values/laprate/timerate of positive/negative matches for self, same traj, same sess, diff sess diff trajs
%2. Plot zscores/p-values/timerate of positive/negative matches for open, sleep, self/other track sessions

hf = gcbf; pinfonow = getappdata(hf, 'pinfo'); datanow = getappdata(hf, 'data');
plotparm = getappdata(hf, 'plotparm');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); grpind = find(groupselection == 1); 
ngroup = numel(grpind); cellind = cell(1, ngroup); grpnames = cell(1, ngroup);
for (kk = 1:numel(grpind)) 
    cellind{kk} = datanow.grouplist.groupindex{grpind(kk)}; grpnames{kk} = datanow.grouplist.groupname{grpind(kk)}; 
end
evtraj = {'self'; 'samTraj'; 'sameSess'; 'allDiff'}; sesstask = {'open'; 'sleep'; 'self'; 'otherLinear'};
evvar = {'evPosMatchNshufsig'; 'evPosMatchShufP'; 'evPosTimeRate'; 'evPosLapRate'; 'evNegMatchNshufsig'; 'evNegMatchShufP'; 'evNegTimeRate'; 'evNegLapRate'};
evsigval = [1.645 0.05 NaN NaN 1.645 0.05 NaN NaN];
sessvar = {'sessPosMatchNshufsig'; 'sessPosMatchShufP'; 'sessPosRate'; 'sessNegMatchNshufsig'; 'sessNegMatchShufP'; 'sessNegRate'};
sesssigval = [1.645 0.05 NaN 1.645 0.05 NaN];
%%%%event sequencing results
ntraj = numel(evtraj);
for (tt = 1:numel(evvar))
    plotval = cell(ngroup, ntraj); mm = NaN*ones(ngroup, ntraj); se = NaN*ones(ngroup, ntraj); 
    md = NaN*ones(ngroup, ntraj);  nn = NaN*ones(ngroup, ntraj); 
    ppgroup = NaN*ones(ngroup-1, ntraj); pptraj = NaN*ones(ngroup, ntraj-1);
    for (i = 1:ngroup)
        val = pinfonow.seq.(evvar{tt})(cellind{i}); tmpname = pinfonow.tmp.sessevname(cellind{i});
        evname = pinfonow.general.eventname(cellind{i}); evtype = pinfonow.parm.eventType(cellind{i});
        for (j = 1:ntraj)
            plotval{i,j} = NaN*ones(1,numel(cellind{i}));
        end
        for (k = 1:numel(cellind{i}))
            [S, T] = strtok(tmpname{k}, '_'); S = lower(S); T = lower(T);
            evnow = evname{k}; evtypenow = evtype{k};
            sessyes = zeros(1, numel(evnow)); trajyes = zeros(1, numel(evnow)); typeyes = zeros(1, numel(evnow));
            for (j = 1:numel(evnow))
                if (~isempty(strfind(lower(evnow{j}), S))) sessyes(j) = 1; end
                if (~isempty(strfind(lower(evnow{j}), T))) trajyes(j) = 1; end
                if strcmp(evtypenow{j}, 'run') typeyes(j) = 1; end
            end
            ii = find(sessyes & trajyes & typeyes); if (numel(ii) == 1) plotval{i,1}(k) = val{k}(ii); end %%%self
            ii = find( (~sessyes) & (trajyes) & typeyes); if (numel(ii) == 1) plotval{i,2}(k) = val{k}(ii); end %%%same traj
            ii = find((sessyes) & (~trajyes) & typeyes); if (numel(ii) == 1) plotval{i,3}(k) = val{k}(ii); end %%%self sess
            ii = find((~sessyes) & (~trajyes) & typeyes); if (numel(ii) == 1) plotval{i,4}(k) = val{k}(ii); end %%%all diff
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
    hf = figure('Name', ['SignificantSequences - ', evvar{tt}]); wid = 1/ngroup; gap = 0.05;
    for (i = 1:ngroup)
        hax(i) = axes('Parent', hf, 'NextPlot', 'add', 'Units', 'normalized', 'Position', [(i-1)*wid+gap 0.05 wid-gap 0.95]);
        for (j = 1:ntraj)
            YY = plotval{i,j}; XX = j*ones(1, nn(i,j)); 
            bar(hax(i), j, mm(i,j), 'k');
            line(XX, YY, 'Parent', hax(i), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30);
            %line([j-0.2 j+0.2], [mm(i,j) mm(i,j)], 'Parent', hax(i), 'LineWidth', 2); %%%mean value line 
        end
        line([0 ntraj+1], [evsigval(tt) evsigval(tt)], 'Parent', hax(i)); %%%significance line
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

%%%%session sequencing results
ntraj = numel(sesstask);
for (tt = 1:numel(sessvar))
    plotval = cell(ngroup, ntraj); mm = NaN*ones(ngroup, ntraj); se = NaN*ones(ngroup, ntraj); 
    md = NaN*ones(ngroup, ntraj);  nn = NaN*ones(ngroup, ntraj); 
    ppgroup = NaN*ones(ngroup-1, ntraj); pptraj = NaN*ones(ngroup, ntraj-1);
    for (i = 1:ngroup)
        val = pinfonow.seq.(sessvar{tt})(cellind{i}); tmpname = pinfonow.tmp.sessevname(cellind{i});
        sessname = pinfonow.general.sessionname(cellind{i}); sesstype = pinfonow.parm.sessType(cellind{i});
        %for (j = 1:ntraj) plotval{i,j} = NaN*ones(1,numel(cellind{i})); end
        for (k = 1:numel(cellind{i}))
            [S, T] = strtok(tmpname{k}, '_'); S = lower(S);
            sessnow = sessname{k}; sesstypenow = sesstype{k};
            sessyes = zeros(1, numel(sessnow)); 
            for (j = 1:numel(sessnow))
                if (~isempty(strfind(S, lower(sessnow{j})))) sessyes(j) = 1; end
            end
            ii = find(strcmp(sesstypenow, 'open')); 
            for (j = 1:numel(ii)) plotval{i,1} = [plotval{i,1} val{k}(ii(j))]; end %%%open
            ii = find(strcmp(sesstypenow, 'sleep')); 
            for (j = 1:numel(ii)) plotval{i,2} = [plotval{i,2} val{k}(ii(j))]; end %%%sleep
            ii = find(strcmp(sesstypenow, 'linear') & sessyes); 
            for (j = 1:numel(ii)) plotval{i,3} = [plotval{i,3} val{k}(ii(j))]; end %%%track self
            ii = find(strcmp(sesstypenow, 'linear') & (~sessyes)); 
            for (j = 1:numel(ii)) plotval{i,4} = [plotval{i,4} val{k}(ii(j))]; end %%%track other
        end
    end
    for (i = 1:ngroup) %%%%basic statistics
        for (j = 1:numel(sesstask))
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
    hf = figure('Name', ['SignificantSequences - ', sessvar{tt}]); wid = 1/ngroup; gap = 0.05;
    for (i = 1:ngroup)
        hax(i) = axes('Parent', hf, 'NextPlot', 'add', 'Units', 'normalized', 'Position', [(i-1)*wid+gap 0.05 wid-gap 0.95]);
        for (j = 1:ntraj)
            YY = plotval{i,j}; XX = j*ones(1, nn(i,j)); 
            bar(hax(i), j, mm(i,j), 'k'); %, 'Parent', hax, 'LineWidth', 2); %%%mean value line 
            line(XX, YY, 'Parent', hax(i), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30);
            %line([j-0.2 j+0.2], [mm(i,j) mm(i,j)], 'Parent', hax(i), 'LineWidth', 2); %%%mean value line 
        end
        line([0 ntraj+1], [sesssigval(tt) sesssigval(tt)], 'Parent', hax(i)); %%%significance line
        str = [grpnames{i}, ': open, sleep, self, otherlinear'];
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
            str = [sesstask{k}, ' group comparisons (t-test against the 1st) p=', num2str(ppgroup(:,k)')]; %, '%10.5e')];
            text('Parent', hax(i), 'Interpreter', 'none', 'String', str, 'Units', 'normalized', 'Position', [0.05 thnow]);
        end
    end
end
disp('*******************************');

