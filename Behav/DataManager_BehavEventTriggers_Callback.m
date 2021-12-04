function DataManager_BehavEventTriggers_Callback
%compute event triggered behavioral features: location distribution,
%triggered speed, 2D and 1D trajectory overlay
hf = gcbf; behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
grpind = find(groupselection == 1); cellind = bhdata.grouplist.groupindex(grpind); groupname = [];
sessind = []; ok = 1;
for (i = 1:numel(grpind))
    sessind = union(sessind, cellind{i}); groupname = [groupname '+' bhdata.grouplist.groupname{grpind(i)}]; 
end
if isempty(sessind) ok = 0; end
if ok
    input = inputdlg({'Event keyword'; 'Event type'; 'Event keyNOword'; 'Event NOtype'}, 'Event selection', 4, {'head'; ''; ''; ''}); 
    if (~isempty(input))
        parm.evkeyword = input{1}; parm.evkeytype = input{2}; parm.evkeynoword = input{3}; parm.evkeynotype = input{4};
    else
        ok = 0;
    end
end
if ok
   tagmark = get(gcbo, 'Tag'); 
   if (strcmp(tagmark, 'eventlocations'))
      input = inputdlg({'Which time (start/ent/ref/whole)'; 'Time range'}, 'Time parameters', 2, {'whole'; '-0.8 0.8'}); 
      if (~isempty(input))
         parm.timetag = input{1}; parm.range = str2num(input{2}); 
         if numel(parm.range)~=2 ok = 0; end
      else
         ok = 0;
      end
      if ok showeventlocations(behav, bhdata, sessind, parm); end
   elseif (strcmp(tagmark, 'eventtrigbehavparm'))
      input = inputdlg({'Which time (start/ent/ref)'; 'Time bin (s)'; 'Time range'; 'Group averages?'}, 'Time parameters', 4, {'start'; '0.1'; '-3 3'; 'Yes'}); 
      if (~isempty(input))
         parm.timetag = input{1}; parm.binsize = str2num(input{2}); parm.range = str2num(input{3}); parm.avgplotflag = 0;
         if strncmpi(input{4}, 'Yes', 1) parm.avgplotflag = 1; end
         if numel(parm.range)~=2 ok = 0; end
      else
         ok = 0;
      end
      if ok eventtriggeredbehavior(behav, bhdata, sessind, groupname, parm); end
   end
end
disp('********************');

function eventtriggeredbehavior(behav, bhdata, sessind, groupname, parm)
%%%%prepare time bins
nn = floor(parm.range/parm.binsize); tbin = [nn(1):-1 0 1:nn(2)]*parm.binsize; 
Y2dVall = ones(0, numel(tbin)); Y2dHDall = ones(0, numel(tbin));
Y1dXall = ones(0, numel(tbin)); Y1dSall = ones(0, numel(tbin));
for (ttt = 1:numel(sessind))
    i = sessind(ttt); sessname = behav.general.sessID{i};
    evName = behav.general.eventname{i}; evType = behav.parm.eventType{i}; evT = bhdata.event.eventtimes{i}; 
    evsel = checkifcorrectevents(evName, evType, parm.evkeyword, parm.evkeytype, parm.evkeynoword, parm.evkeynotype);
    evind= find(evsel);
    T = [];
    for (j = 1:numel(evind))
        T = union(T, evT{evind(j)}.(parm.timetag)); 
    end
    T = sort(T);
    timestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i); 
%     %%%%%2D trajectories - not sure what to do with this - skip this    
%     Pmarker = behav.parm.sessPmarker{i}; allposmarker = behav.general.posMarker{i}; 
%     posXX = []; posYY = []; ik = find(strcmp(allposmarker, Pmarker)); 
%     if (numel(ik) == 1) posXX = bhdata.pos.XX{i}{ik}*behav.parm.pixelXSize(i); posYY = bhdata.pos.YY{i}{ik}*behav.parm.pixelYSize(i); end
    %%%%%2D speed, headdir
    allV = bhdata.sess.AllVel{i}; allHdir = bhdata.sess.AllHeadDir{i}; 
    [tbin, Y, Yavg] = getandplot2Davgs(tbin, timestamp, allV, T, '2D velocity (cm; mean+-err)', sessname, ~parm.avgplotflag); %%%%2D velocity: always postive
    if parm.avgplotflag Y2dVall = [Y2dVall; Y]; end
    [tbin, Y, Yavg] = getandplot2Davgs(tbin, timestamp, allHdir, T, '2D head dir (o; mean+-var)', sessname, ~parm.avgplotflag); %%%%2D head dir: [0 360] circular
    if parm.avgplotflag Y2dHDall = [Y2dHDall; Y]; end
    %%%%%1D positon, speed
    laptime = bhdata.event.LapAllPostimestamp{i}; %%{ntraj}{nlap}
    iiev = find(strcmp(evType, 'run')); laptime = laptime(iiev); 
    lapX = bhdata.event.LapAllX{i}(iiev);
    speed1d = bhdata.event.LapAll1DSpeed{i}(iiev);%%{ntraj}{nlap}
    [tbin, Y, Yavg] = getandplot1Davg(tbin, laptime, lapX, T, '1D Positon (cm)', sessname, ~parm.avgplotflag); %%%%1D position
    if parm.avgplotflag Y1dXall = [Y1dXall; Y]; end    
    [tbin, Y, Yavg] = getandplot1Davg(tbin, laptime, speed1d, T, '1D speed (cm/s)', sessname, ~parm.avgplotflag); %%%%1D speed 
    if parm.avgplotflag Y1dSall = [Y1dSall; Y]; end        
end
if parm.avgplotflag
    plot2Davgsonly(tbin, Y2dVall, '2D velocity (cm; mean+-err)', groupname);
    plot2Davgsonly(tbin, Y2dHDall, '2D head dir (o; mean+-var)', groupname);
    plot1Davgsonly(tbin, Y1dXall, '1D Positon (cm)', groupname);
    plot1Davgsonly(tbin, Y1dSall, '1D speed (cm/s)', groupname);
end

function showeventlocations(behav, bhdata, sessind, parm)
for (ttt = 1:numel(sessind))
    i = sessind(ttt); sess = behav.general.sessID{i};
    %%%%2D locations to show track
    postimestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i); 
    Pmarker = behav.parm.sessPmarker{i}; allposmarker = behav.general.posMarker{i}; 
    posXX = []; posYY = []; ik = find(strcmp(allposmarker, Pmarker)); 
    if (numel(ik) == 1) posXX = bhdata.pos.XX{i}{ik}*behav.parm.pixelXSize(i); posYY = bhdata.pos.YY{i}{ik}*behav.parm.pixelYSize(i); end
    if (mean(posXX)>0)
        hf = figure('Name', strcat(sess, '---event 2D location')); hax = axes('Parent', hf, 'NextPlot', 'add'); str = sess;
        line(posXX, posYY, 'Parent', hax, 'LineStyle', 'none',...
                            'Marker', '.', 'MarkerSize', 2, 'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', [0.8 0.8 0.8]); % 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
        text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
        xlabel('X (cm)'); ylabel('Y (cm)'); 
    end
    %%%%identify events
    evName = behav.general.eventname{i}; evType = behav.parm.eventType{i}; evT = bhdata.event.eventtimes{i}; 
    evsel = checkifcorrectevents(evName, evType, parm.evkeyword, parm.evkeytype, parm.evkeynoword, parm.evkeynotype);
    evind= find(evsel); T = []; Bf = []; Af = []; nt= 0;
    if ~strncmpi(parm.timetag, 'who', 3)
       for (j = 1:numel(evind)) T = union(T, evT{evind(j)}.(parm.timetag)); end 
       T = sort(T);
       [x, y] = findtime2Dlocations(T, postimestamp, posXX, posYY);
       line(x, y, 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
    else
       for j = 1:numel(evind)
           for k = 1:numel(evT{evind(j)}.start)
               nt = nt+1;
               T{nt} = postimestamp( (postimestamp>=evT{evind(j)}.start(k)) & (postimestamp<=evT{evind(j)}.ent(k)) );
               Bf{nt} = postimestamp( (postimestamp>=evT{evind(j)}.start(k)+parm.range(1)) & (postimestamp<evT{evind(j)}.start(k)) );
               Af{nt} = postimestamp( (postimestamp>evT{evind(j)}.ent(k)) & (postimestamp<=evT{evind(j)}.ent(k)+parm.range(2)) );
               [xt, yt] = findtime2Dlocations(T{nt}, postimestamp, posXX, posYY);
               [xb, yb] = findtime2Dlocations(Bf{nt}, postimestamp, posXX, posYY); [xa, ya] = findtime2Dlocations(Af{nt}, postimestamp, posXX, posYY);
               line(xt, yt, 'Parent', hax, 'LineWidth', 1, 'Color', [0 0 0]);
               line(xb, yb, 'Parent', hax, 'LineWidth', 1, 'Color', [0 1 0]); line(xa, ya, 'Parent', hax, 'LineWidth', 1, 'Color', [1 0 0]);
               line(xt(1), yt(1), 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]);
               line(xt(numel(xt)), yt(numel(yt)), 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
           end
       end
    end
    %%%%1D positions
    laptime = bhdata.event.LapAllPostimestamp{i}; %%{ntraj}{nlap}
    xj = bhdata.event.Xjoint{i}; iiev = find(strcmp(evType, 'run')); evname = evName(iiev); 
    lapX = bhdata.event.LapAllX{i}(iiev); laptime = laptime(iiev); xj = xj(iiev); %%plot only those 'run' events
    plotlinearnow(laptime, lapX, laptime, xj, sess, evname, T, Bf, Af, parm.timetag);
end

function [x, y] = findtime2Dlocations(T, timestamp, xpos, ypos)
lastpoint = 1;  %this only for saving time
x = NaN*size(T); y = NaN*size(T);
for (j = 1:numel(T))
    for (k = lastpoint:numel(xpos)-1) %find corresponding time in position data
         if ((timestamp(k) <= T(j)) && (timestamp(k+1) > T(j))) 
             x(j) = xpos(k); y(j) = ypos(k); 
             lastpoint = k;
             break; 
         end
    end
end

function x = findtime1Dlocations(T, timestamp, xpos)
lastpoint = 1;  %this only for saving time
x = NaN*ones(size(T));
for (j = 1:numel(T))
    for (k = lastpoint:numel(xpos)-1) %find corresponding time in position data
         if ((timestamp(k) <= T(j)) && (timestamp(k+1) > T(j))) 
             x(j) = xpos(k);  
             lastpoint = k;
             break; 
         end
    end
end

function hax = plotlinearnow(lapX, lapY, lapT, xj, sess, evname, T, Bf, Af, timetag)
nev = numel(lapX);
if nev==2
   if (lapT{1}{1}(1) > lapT{2}{1}(1)) %%%reverse the values of the second trajectory
       mV = xj{1}(numel(xj{1})); xj{1} = mV - xj{1}; for (ki = 1:numel(lapY{1})) lapY{1}{ki} = mV - lapY{1}{ki}; end 
   else
       mV = xj{2}(numel(xj{2})); xj{2} = mV - xj{2}; for (ki = 1:numel(lapY{2})) lapY{2}{ki} = mV - lapY{2}{ki}; end 
   end
end
minY = 0; maxY = 1;
for (i = 1:numel(lapY))
    for (j = 1:numel(lapY{i}))  
         if (minY > min(lapY{i}{j})) minY = min(lapY{i}{j}); end
         if (maxY < max(lapY{i}{j})) maxY = max(lapY{i}{j}); end
    end
end
if (nev>0)
    hf = figure('Name', strcat(sess, '--1D event location')); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
    for (i = 1:nev)
         for (j = 1:numel(lapY{i}))
             line(lapX{i}{j}, lapY{i}{j}, 'Parent', hax, 'LineStyle', 'none',...
                            'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', [0.8 0.8 0.8]); % 'Color', col{i}, 'LineWidth', 2);  
             set(hax, 'YLim', [minY maxY]); %set(hax, 'YLim', [-10 500]);
             for (k = 1:numel(xj{i}))
                 if ~isempty(lapX{i}{j})
                      line([min(lapX{i}{j}) max(lapX{i}{j})], [xj{i}(k) xj{i}(k)], 'Parent', hax, 'Color', [0 0 0],'LineWidth', 1); 
                 end
             end
         end
         
    end
    if ~strncmpi(timetag, 'who', 3)
        for (i=1:nev)
            for (j = 1:numel(lapY{i}))
                iii = find( (T>=lapT{i}{j}(1)) & (T<=lapT{i}{j}(numel(lapT{i}{j}))) );
                xnow = findtime1Dlocations(T(iii), lapT{i}{j}, lapY{i}{j}); 
                line(T(iii), xnow, 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
            end
        end
    else
        for (k = 1:numel(T))
            for (i=1:nev)
            for (j = 1:numel(lapY{i}))
                iii = find( (T{k}>=lapT{i}{j}(1)) & (T{k}<=lapT{i}{j}(numel(lapT{i}{j}))) );
                if ~isempty(iii)
                   tnow = T{k}(iii); xnow = findtime1Dlocations(tnow, lapT{i}{j}, lapY{i}{j}); 
                   line(tnow, xnow, 'Parent', hax, 'LineWidth', 1, 'Color', [0 0 0]);
                   line(tnow(1), xnow(1), 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]);
                   line(tnow(numel(tnow)), xnow(numel(xnow)), 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
                   jjj = find( (Bf{k}>=lapT{i}{j}(1)) & (Bf{k}<=lapT{i}{j}(numel(lapT{i}{j}))) );
                   if ~isempty(jjj)
                      tnow = Bf{k}(jjj); xnow = findtime1Dlocations(tnow, lapT{i}{j}, lapY{i}{j}); 
                      line(tnow, xnow, 'Parent', hax, 'LineWidth', 1, 'Color', [0 1 0]);
                   end
                   jjj = find( (Af{k}>=lapT{i}{j}(1)) & (Af{k}<=lapT{i}{j}(numel(lapT{i}{j}))) );
                   if ~isempty(jjj)
                      tnow = Af{k}(jjj); xnow = findtime1Dlocations(tnow, lapT{i}{j}, lapY{i}{j}); 
                      line(tnow, xnow, 'Parent', hax, 'LineWidth', 1, 'Color', [1 0 0]);
                   end
                end
            end
            end
        end
    end
    for (i = 1:nev)
       text('Interpreter', 'none', 'Parent', hax, 'String', evname{i}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
    end 
    xlabel('Time (s)'); ylabel('Position (cm)');
end

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

function [tbin, Y, Yavg] = getandplot2Davgs(tbin, timestamp, allV, T, tag, sessname, plotflag)
ntrace = numel(T); nt = numel(tbin); binsize = mean(diff(tbin));
Y = NaN*ones(ntrace, nt); Yavg = NaN*ones(1, nt); E = NaN*ones(1, nt); refV = NaN*ones(1, nt);
for (i= 1:ntrace) %%%reference - time 0 - values
    iii = find( (timestamp>=T(i)-binsize/2) & (timestamp<T(i)+binsize/2) ); 
    ynow = allV(iii); ynow =ynow(~isnan(ynow));
    if ~isempty(ynow) refV(i) = mean(ynow); end
end
if ~contains(tag, 'dir') %%%%if not 2D head dir, set ref values to zeros - do not do re-setting
    iii = find(~isnan(refV)); refV(iii) = zeros(1, numel(iii)); 
end
if ~contains(tag, 'dir')
   for (j = 1:nt)
    for (i = 1:ntrace)
        tp = T(i)+tbin(j);
        iii = find( (timestamp>=tp-binsize/2) & (timestamp<tp+binsize/2) ); 
        ynow = allV(iii); ynow =ynow(~isnan(ynow));
        if ~isempty(ynow) Y(i,j) = mean(ynow) - refV(i); end
    end
    aa = Y(:,j); aa = aa(~isnan(aa));
    if ~isempty(aa) Yavg(j) = mean(aa); E(j) = std(aa)/sqrt(numel(aa)); end
   end
else %%%%head dir uses circular mean -- also need to reset time 0 dir = 0o
   for (i= 1:ntrace) %%%reference - time 0 - values
        iii = find( (timestamp>=T(i)-binsize/2) & (timestamp<T(i)+binsize/2) ); 
        ynow = allV(iii); ynow =ynow(~isnan(ynow));
        if ~isempty(ynow) refV(i) = meandirnow(ynow); end
   end
   for (j = 1:nt)
    for (i = 1:ntrace)
        tp = T(i)+tbin(j);
        iii = find( (timestamp>=tp-binsize/2) & (timestamp<tp+binsize/2) ); 
        ynow = allV(iii); ynow =ynow(~isnan(ynow));
        if ~isempty(ynow) Y(i,j) = mod((meandirnow(ynow) - refV(i) + 180), 360); end %%%%reference = 180o
    end
    aa = Y(:,j); aa = aa(~isnan(aa));
    if ~isempty(aa) [Yavg(j), E(j), ~, ~] = findphaseprop(aa); end
   end
end
if plotflag
   hf = figure('Name', strcat(sessname, '-2D_', tag)); 
   hax = axes('Parent', hf, 'NextPlot', 'add'); xlabel('Time (s)'); ylabel(tag);
   for (i = 1:ntrace)
       line(tbin, Y(i,:), 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
   end
   line(tbin, Yavg, 'Parent', hax, 'LineWidth', 2, 'Color', [1 0 0]);
   Drawerrorupdown(tbin, Yavg, Yavg+E, Yavg-E, hax, [1 0 0]);
   text('Interpreter', 'none', 'Parent', hax, 'String', sessname, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
   text('Interpreter', 'none', 'Parent', hax, 'String', ['N = '  num2str(numel(find(~isnan(refV))))], 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92]);
end

function plot2Davgsonly(tbin, Y, tag, sessname)
[ntrace, nt] = size(Y); Yavg = NaN*ones(1, nt); E = NaN*ones(1, nt); 
if ~contains(tag, 'dir')
   for (j = 1:nt)
    aa = Y(:,j); aa = aa(~isnan(aa));
    if ~isempty(aa) Yavg(j) = mean(aa); E(j) = std(aa)/sqrt(numel(aa)); end
   end
else %%%%head dir uses circular mean -- also need to reset time 0 dir = 0o
   for (j = 1:nt)
    aa = Y(:,j); aa = aa(~isnan(aa));
    if ~isempty(aa) [Yavg(j), E(j), ~, ~] = findphaseprop(aa); end
   end
end
hf = figure('Name', strcat(sessname, '-2D_', tag)); 
hax = axes('Parent', hf, 'NextPlot', 'add'); xlabel('Time (s)'); ylabel(tag);
for (i = 1:ntrace)
     line(tbin, Y(i,:), 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
end
line(tbin, Yavg, 'Parent', hax, 'LineWidth', 2, 'Color', [1 0 0]);
Drawerrorupdown(tbin, Yavg, Yavg+E, Yavg-E, hax, [1 0 0]);
text('Interpreter', 'none', 'Parent', hax, 'String', sessname, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
text('Interpreter', 'none', 'Parent', hax, 'String', ['N = '  num2str(ntrace)], 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92]);
[str, X] = findsignificance(tbin, Yavg, E, Y);
for i = 1:numel(str)
    text('Interpreter', 'none', 'Parent', hax, 'String', str{i}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92-i*0.04]);
end
%%%save data to figure
setappdata(hf, 'X', X);

function [str, X] = findsignificance(xbin, crr, err, allcrr)
%%%%assess signficance: allcrr(ntrace, nbin);
%%%%%% test peak crr significance (and from zero)
[nsess, nbin] = size(allcrr); parm.basetime = [-1 -0.8]; parm.biastime = [-0.3 0.3];
X.alltroughR = NaN*ones(1, nsess); X.alltroughT = NaN*ones(1, nsess); basetime = parm.basetime; biastime = parm.biastime;
X.allbias = NaN*ones(1, nsess); X.allbase = NaN*ones(1, nsess); X.parm = parm;
valnow = NaN*ones(1, nbin*nsess); grp = NaN*ones(1, nbin*nsess);
for (i = 1:nsess)
     for (j = 1:nbin)
          valnow(j+(i-1)*nsess) = allcrr(i,j); grp(j+(i-1)*nsess) = j;
     end
end
iii = find(~isnan(valnow)); valnow = valnow(iii); grp = grp(iii);
[ppp, table, stats] = anova1(valnow, grp, 'off'); 
str{1} = ['anova1 p=', num2str(ppp, '%10.5e')];
[peakr,iii] = max(crr);[troughr,jjj] = min(crr); 
[~,peakp] = ttest(allcrr(:,iii)); [~,troughp] = ttest(allcrr(:,jjj));
str{2} = ['peak: r=', num2str(peakr), '; se=', num2str(err(iii)), '; t=', num2str(xbin(iii)), '; p(from 0)=', num2str(peakp, '%10.5e')];
str{3} = ['trough: r=', num2str(troughr), '; se=', num2str(err(jjj)),  '; t=', num2str(xbin(jjj)), '; p(from 0)= ', num2str(troughp, '%10.5e')];
%%%%%peaks/troughs
allpeakbin = find( (xbin>=min(biastime)) & (xbin<max(biastime)) ); 
allcrr = allcrr'; %%now (xbin, nsess/ncell/nevent)
crrnow = allcrr(allpeakbin,:);
[X.allpeakR, iii] = max(crrnow); X.allpeakT = xbin(allpeakbin(iii));  
[X.alltroughR, jjj] = min(crrnow); X.alltroughT = xbin(allpeakbin(jjj));
leftbin = find( (xbin>=min(biastime)) & (xbin<mean(biastime)) ); rightbin = find( (xbin<=max(biastime)) & (xbin>mean(biastime)) ); 
if isempty(leftbin) || isempty(rightbin)
   disp('----> warning: at least one bias bin is empty');
else
   for ij = 1:nsess
       aa =  allcrr(leftbin,ij); leftB = mean(aa(~isnan(aa))); 
       aa =  allcrr(rightbin,ij); rightB = mean(aa(~isnan(aa))); 
       X.allbias(ij) = (rightB-leftB)/(rightB+leftB);
   end
   aa = X.allbias; aa = aa(~isnan(aa));
   mm = mean(aa); nn = numel(aa); se = std(aa)./sqrt(nn); [~,biasp] = ttest(aa);
   str{4} = ['bias: mean=', num2str(mm), '; se=', num2str(se), '; n=', num2str(nn), '; p (from 0) = ', num2str(biasp, '%10.5e')];
end
%%%baseline
leftbin = find( (xbin>=min(basetime)) & (xbin<=max(basetime)) ); rightbin = find( (xbin>=min(-1*basetime)) & (xbin<=max(-1*basetime)) ); 
allbasebin = union(leftbin, rightbin);
if ~isempty(allbasebin)
   if numel(allbasebin)>1
       for ij = 1:nsess
           aa = allcrr(allbasebin,ij); aa = mean(aa(~isnan(aa))); 
           X.allbase(ij) = aa;
       end
   else
      X.allbase = allcrr(allbasebin,:);
   end
   aa = X.allbase; aa = aa(~isnan(aa)); 
   mm = mean(aa); nn = numel(aa); se = std(aa)./sqrt(nn); [~,basep] = ttest2(aa, X.allpeakR);
   str{5} = ['baseline: mean=', num2str(mm), '; se=', num2str(se), '; n=', num2str(nn), '; p (from peak) = ', num2str(basep, '%10.5e')];
end

function [tbin, Y, Yavg] = getandplot1Davg(tbin, laptime, lapX, T, tag, sessname, plotflag)
ntrace = numel(T); nt = numel(tbin); binsize = mean(diff(tbin));
Y = NaN*ones(ntrace, nt); Yavg = NaN*ones(1, nt); E = NaN*ones(1, nt); refV = NaN*ones(1, ntrace);
%%%reference - time 0 - values
for (i= 1:ntrace) %%%reference - time 0 - values
     ynow = findvalues(laptime, lapX, T(i), binsize);
     if ~isempty(ynow) refV(i) = mean(ynow); end
end
if ~contains(tag, 'Pos') %%%%if not 1D postion, set ref values to zeros - do not do re-setting
    iii = find(~isnan(refV)); refV(iii) = zeros(1, numel(iii)); 
end
%%%values for all time points
for (j = 1:nt)
    for (i = 1:ntrace)
        ynow = findvalues(laptime, lapX, T(i)+tbin(j), binsize);
        if ~isempty(ynow) Y(i,j) = mean(ynow) - refV(i); end
    end
    aa = Y(:,j); aa = aa(~isnan(aa));
    if ~isempty(aa) Yavg(j) = mean(aa); E(j) = std(aa)/sqrt(numel(aa)); end
end
if contains(tag, 'Pos') && contains(tag, '1D') %%%for 1D position: plot as always get in to ref from bottom left 
   for (i=1:ntrace)
       aa = Y(i,1:5); aa = aa(~isnan(aa)); 
       if ~isempty(aa) && mean(aa)>0 Y(i,:) = -Y(i,:); end
   end
end
if plotflag
   hf = figure('Name', strcat(sessname, '-1D_', tag)); 
   hax = axes('Parent', hf, 'NextPlot', 'add'); xlabel('Time (s)'); ylabel(tag);
   for (i = 1:ntrace)
       line(tbin, Y(i,:), 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
   end
   line(tbin, Yavg, 'Parent', hax, 'LineWidth', 2, 'Color', [1 0 0]);
   Drawerrorupdown(tbin, Yavg, Yavg+E, Yavg-E, hax, [1 0 0]);
   text('Interpreter', 'none', 'Parent', hax, 'String', sessname, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
   text('Interpreter', 'none', 'Parent', hax, 'String', ['N = ' num2str(numel(find(~isnan(refV))))], 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92]);
end
function plot1Davgsonly(tbin, Y, tag, sessname)
[ntrace, nt] = size(Y); Yavg = NaN*ones(1, nt); E = NaN*ones(1, nt); 
%%%values for all time points
for (j = 1:nt)
    aa = Y(:,j); aa = aa(~isnan(aa));
    if ~isempty(aa) Yavg(j) = mean(aa); E(j) = std(aa)/sqrt(numel(aa)); end
end
if contains(tag, 'Pos') && contains(tag, '1D') %%%for 1D position: plot as always get in to ref from bottom left 
   for (i=1:ntrace)
       aa = Y(i,1:5); aa = aa(~isnan(aa)); 
       if ~isempty(aa) && mean(aa)>0 Y(i,:) = -Y(i,:); end
   end
end
hf = figure('Name', strcat(sessname, '-1D_', tag)); 
hax = axes('Parent', hf, 'NextPlot', 'add'); xlabel('Time (s)'); ylabel(tag);
for (i = 1:ntrace)
     line(tbin, Y(i,:), 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
end
line(tbin, Yavg, 'Parent', hax, 'LineWidth', 2, 'Color', [1 0 0]);
Drawerrorupdown(tbin, Yavg, Yavg+E, Yavg-E, hax, [1 0 0]);
text('Interpreter', 'none', 'Parent', hax, 'String', sessname, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
text('Interpreter', 'none', 'Parent', hax, 'String', ['N = ' num2str(ntrace)], 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92]);
[str, X] = findsignificance(tbin, Yavg, E, Y);
for i = 1:numel(str)
    text('Interpreter', 'none', 'Parent', hax, 'String', str{i}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92-i*0.04]);
end
%%%save data to figure
setappdata(hf, 'X', X);

function ynow = findvalues(laptime, lapX, T, binsize) %%%laptime, lapX contain multiple trajectories {ntraj}{nlap}(npoints)
iii = []; ynow = [];
for (i = 1:numel(lapX))
    for (j = 1:numel(lapX{i}))
        jjj = find( (laptime{i}{j}>=T-binsize/2) & (laptime{i}{j}<T+binsize/2) );
        if numel(jjj)>numel(iii)
           ynow = lapX{i}{j}(jjj); ynow =ynow(~isnan(ynow));
        end
    end
end

function [meanphase, phasevar, Rayleng, RayP] = findphaseprop(phasenow)
meanphase = NaN; phasevar = NaN; Rayleng = NaN; RayP = NaN;
nspike = numel(phasenow); %phase in [0 360)
if (nspike > 1)
   CC = exp(1i*phasenow*2*pi/360); %%now a complex number
   %first do the corcular statistics 
   finalV = sum(CC);
   Rayleng = abs(finalV); meanphase = angle(finalV); %%here phase is in [-pi pi), need to change to [0 360)
   if (meanphase >= 0)
       meanphase = 180 * meanphase / pi;
   else
       meanphase = 360 + (180 * meanphase / pi);
   end
   phasevar = 1 - Rayleng/nspike;
   %next do the Rayleigh test
   dist = zeros(1, nspike); %distance from origin
   for (k = 1:nspike)
       dist(k) = abs( sum(CC(1:k)) ); %samples of distance
   end
   b = raylfit(dist); %estimate of Rayleigh distribution parameter
   RayP = 1 - raylcdf(Rayleng, b); %final vector length's probability from Rayleigh distribution
end

function meanphase = meandirnow(phasenow) %phasenow in [0 360]o
CC = exp(1i*phasenow*2*pi/360); %%now a complex number
meanphase = angle(sum(CC)); %%here phase is in [-pi pi), need to change to [0 360)
if (meanphase >= 0)
   meanphase = 180 * meanphase / pi;
else
   meanphase = 360 + (180 * meanphase / pi);
end

