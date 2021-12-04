function DataManager_EEGShowEEGdata_Callback
%%Plot the computed eegdata
hf = gcbf; eeg = getappdata(hf, 'eeg'); eegdata = getappdata(hf, 'eegdata');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
hfield = getappdata(hf, 'hfield'); plotparm = getappdata(hf, 'plotparm');
hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); 
dateind = find(spikeselection == 1); %%%%%dateind here is actually eegind (eeg IDs)
plotparm = getappdata(hf, 'plotparm');
grpind = find(groupselection == 1); grpname = eegdata.grouplist.groupname(grpind);
hfield = getappdata(hf, 'hfield');

tagmark = get(gcbo, 'Tag'); 
if (strcmp(tagmark, 'spectrogram'))
    plotPowerSpectrogram(eeg, eegdata, dateind, 'PSD', plotparm);
elseif (strcmp(tagmark, 'PSD'))
    plotPSD(eeg, eegdata, dateind, 'Power', plotparm);
elseif (strcmp(tagmark, 'signaldist'))
    plotsignaldist(eeg, eegdata, dateind, 'Raw signal', plotparm);
elseif (strcmp(tagmark, 'crosscorrelation'))
    plotcrosscorrelation(eeg, eegdata, dateind, 'crosscorrelation', plotparm);
elseif (strcmp(tagmark, 'triggeraverage'))
    plottriggeredaverages(hfield, eeg, eegdata, dateind, 'triggeredaverage', plotparm, grpname);
elseif (strcmp(tagmark, 'triggerincidence'))
    plottriggeredincidence(hfield, eeg, eegdata, dateind, 'triggeredincidence', plotparm, grpname);    
elseif (strcmp(tagmark, 'showeegdata')) || (strcmp(tagmark, 'showcluster0')) %%%for all other options to plot data
     showalldata(eeg, eegdata, dateind, 'eegdata');
end
disp('************************');

function plotcrosscorrelation(eeg, eegdata, dateind, tmark, plotparm)
p = {'Auto/cross?'; 'Bin size (s)'; 'Max lag (s)'}; ok = 1;
d = {'cross'; '1'; '50'};
II = inputdlg(p, 'Correlation parameters', 3, d, 'on'); %%%resizable window
if ~isempty(II)
      crrmode = II{1}; binsize = str2num(II{2}); maxlag = str2num(II{3}); 
else
      ok = 0;
end
if (plotparm.evselect == 1) && (~isempty(dateind))%%if do it for events
   i = dateind(1);
   evnames = eeg.general.eventname{i}; 
   [sss, ok] = listdlg('ListString', evnames, 'PromptString', 'Choose one or two events to show', 'SelectionMode', 'multiple');
   if ok
      eventnames = evnames(sss);
   end
end
if ok
    if strncmpi(crrmode, 'cross', 2) && (numel(dateind)<2)
        disp('----> no more than two traces are selected'); ok = 0;
    else
        ind = 1:min([2 numel(dateind)]);
        dateind = dateind(ind); %%%only do it for the first 2 selections
        if strncmpi(crrmode, 'auto', 2)
            for (j = 1:numel(dateind))
                i = dateind(j);
                if ~strcmp(eeg.parm.band{i}, 'cluster0')
                    fname = eeg.general.eegfile{i};
                    [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
                else
                    dat = eegdata.cluster0.cnt{i}; timestamp = eegdata.cluster0.timepoint{i};
                    fs = 1/eeg.parm.cl0Binsize(i); fname = eeg.general.eegfile{i};
                end
                if (plotparm.evselect == 0)
                    sT = eeg.general.sessstartT{i}; eT = eeg.general.sessendT{i}; 
                    plotcrrnow(dat, dat, timestamp, fs, binsize, maxlag, tmark, sT, eT, 'Whole session', fname, fname)
                else
                    for (k = 1:numel(eventnames))
                        evind = find(strcmp(eeg.general.eventname{i}, eventnames{k})); 
                        if numel(evind)==1
                           evtTimes = eegdata.event.eventtimes{i}{evind};
                           plotcrrnow(dat, dat, timestamp, fs, binsize, maxlag, tmark, evtTimes.start, evtimes.ent, eventnames{k}, fname, fname);
                        end
                    end
                end
            end
        else
            fdir1 = eeg.general.finaldir{dateind(1)}; fdir2 = eeg.general.finaldir{dateind(2)};
            sess1 = eeg.general.sessname{dateind(1)}; sess2 = eeg.general.sessname{dateind(2)};
            if strcmp(fdir1, fdir2) && strcmp(sess1, sess2)
                for (j = 1:numel(dateind))
                    i = dateind(j); fname{j} = eeg.general.eegfile{i};
                    if ~strcmp(eeg.parm.band{i}, 'cluster0')
                       [timestamp{j}, dat{j}, gain, fs(j)] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
                    else
                       dat{j} = eegdata.cluster0.cnt{i}; fs(j) = 1/eeg.parm.cl0Binsize(i);  timestamp{j} = eegdata.cluster0.timepoint{i};
                    end
                end
                if (numel(dat{1})==numel(dat{2})) && (fs(1)==fs(2))
                    if (plotparm.evselect == 0)
                        sT = eeg.general.sessstartT{i}; eT = eeg.general.sessendT{i}; 
                        plotcrrnow(dat{1}, dat{2}, timestamp{1}, fs(1), binsize, maxlag, tmark, sT, eT, 'Whole session', fname{1}, fname{2});
                    else
                        for (k = 1:numel(eventnames))
                            evind = find(strcmp(eeg.general.eventname{i}, eventnames{k})); 
                            if numel(evind)==1                           
                               evtTimes = eegdata.event.eventtimes{i}{evind};                              
                               plotcrrnow(dat{1}, dat{2}, timestamp{1}, fs(1), binsize, maxlag, tmark, evtTimes.start, evtTimes.ent, eventnames{k}, fname{1}, fname{2});
                            end
                        end
                    end
                else
                    disp('----> number of data points or sampling frequency does not match'); ok = 0; 
                end
            else
                 disp('----> first two traces not in the same session'); ok = 0;
            end
        end
    end
end

function plottriggeredaverages(hfield, eeg, eegdata, dateind, tmark, plotparm, grpname)
catname = fieldnames(eeg); ntimevar = 0; triggertime = []; triggername = []; triggerstr = []; ok = 1; nonevent = 0;
%%%determine time varibales: selected event times or selected time variables
for (i = 1:numel(catname))
    subfield = fieldnames(eeg.(catname{i}));
    fieldselection = getappdata(hfield(i), 'selection');
    for (j = 1:numel(subfield))
        if (fieldselection(j) == 1) %%if a subfield is selected
            varnow = eeg.(catname{i}).(subfield{j});
            if (~strcmp(subfield{j}, 'eventname')) || (~strcmp(catname{i}, 'general'))
                for (k = 1:numel(varnow))
                    if isnumeric(varnow{k}) && (~isempty(varnow{k}))
                        ntimevar = ntimevar + 1; triggername{ntimevar} = subfield{j}; nonevent = 1;
                        triggertime{ntimevar} = varnow; triggerstr{ntimevar} = []; break;
                    end
                end
            else %%%if select an event file
                unievent = [];
                for (k = 1:numel(varnow))
                    if ~isempty(varnow{k})
                       unievent = unique(union(unievent, lower(varnow{k})));
                    end
                end
                [sel,ok] = listdlg('ListString', unievent, 'SelectionMode', 'single', 'PromptString', 'Select an event to average');
                if (ok)
                    ntimevar = ntimevar + 1; triggername{ntimevar} = subfield{j}; 
                    triggertime{ntimevar} = unievent(sel); %this is just an event name, no time data here yet
                    lststr = {'Start time'; 'End time'; 'Reference time'};
                    [sel,ok] = listdlg('ListString', lststr, 'SelectionMode', 'single', 'PromptString', 'Select which time to average');
                    if ok
                       triggerstr{ntimevar} = lststr{sel};
                    end
                end
            end
        end
    end
end
%%%determine average parameters
if ok
  if ntimevar == 0
    disp('----> no time variables selected'); ok = 0;
  else
    p = {'Binsize/Maxlag (s)'; 'Smooth bins/times (bins/0=no smoothing)'; 'Bias interval: min/max (s)]'; 'Baseline min/max (s)'; 'session average?'};
    d = {'0.05 5'; '10 10'; '-5 -3'; '-2 2'; 'no'};
    II = inputdlg(p, 'Average parameters', 5, d, 'on'); %%%resizable window
    if ~isempty(II)
       aa = str2num(II{1}); parm.binsize =aa(1); parm.maxlag = aa(2); 
       aa = str2num(II{2}); parm.smoothbin =aa(1); parm.smoothtimes = aa(2);
       parm.biastime = str2num(II{3}); parm.basetime = str2num(II{4}); parm.isavg = 0; if strncmpi(II{5}, 'yes', 1) parm.isavg = 1; end
    else
       ok = 0;
    end 
  end
end
if ok
    evkeyword = []; evkeytype = [];
  %if nonevent %%%if at least one of the trigger variable is not selected from an event file
    input = inputdlg({'Event keyword';'Event type'}, 'Filtering through another event?', 2, {''; ''}); 
    %%%%%3rd question: yes = if average over all trigger times first (therefore generate one trace/sample per session/file)
    if (~isempty(input))
        evkeyword = input{1}; evkeytype = input{2}; 
    else
        ok = 0;
    end
    if ok && ((~isempty(evkeyword)) && (~isempty(evkeytype)))
        disp(['-------> filtering through events ', evkeyword, ' + ', evkeytype]);
    elseif ok
        disp(['-------> no further event filtering ']);
    end
end
%%%%determine what to average
if ok
   spikedbfile = [];
   button = questdlg('Average a spike database?');
   if (strcmp(button, 'Yes')) %%%ask if want to average individual or groups of spikes
      cpathname = fullfile(cd, '*.spikedb'); 
      [fname, pname] = uigetfile(cpathname, 'Select a spike database file to open:');
      if (fname ~= 0)
          spikedbfile = fullfile(pname, fname); 
          S = load(spikedbfile, '-mat'); pinfo = S.pinfo; data = S.data; S = []; %load the plot structure
          allgroupname = data.grouplist.groupname; 
          [sel,ok] = listdlg('ListString', allgroupname, 'SelectionMode', 'single', 'PromptString', 'Select a group to average');
          if (ok)
             dataindthere = data.grouplist.groupindex{sel};
          end
      else
          ok = 0;
      end
   elseif strcmp(button, 'Cancel') || (isempty(button))
        ok = 0;
   end
end
if ok
for (tt = 1:ntimevar)
    avg = []; xbin = [];
    for (j = 1:numel(dateind))
         i = dateind(j);
         fdir = eeg.general.finaldir{i}; sess = eeg.general.sessname{i}; 
         evName = eeg.general.eventname{i}; evType = eeg.parm.eventtype{i}; evT = eegdata.event.eventtimes{i}; 
         if strcmp(triggername{tt}, 'eventname') && (~isempty(triggerstr{tt})) %if select an event file to trigger
            triggertimenow = matchEventtimeout(evT, evType, evName, triggertime{tt}, triggerstr{tt});
            seleventname = strcat('Trigger: ', triggertime{tt}, '+', triggerstr{tt}, '; in events:', evkeyword, '_', evkeytype); 
            %disp(triggertimenow(1:5)');
         else %%if select a time variable to trigger
            triggertimenow = matchtimeout(eeg, eegdata, fdir, sess, triggertime{tt});
            seleventname = strcat('Trigger: ', triggername{tt}, '+', triggerstr{tt}, '; in events:', evkeyword, '_', evkeytype); 
            %disp(triggertimenow(1:5));
         end
         triggertimenow = filtertriggertime(triggertimenow,evkeyword, evkeytype, evName, evType, evT); %%% if ketwords/types specified, filter through these events
         if isempty(triggertimenow)
             disp(['-------> warning: trigger times not found in ', fdir, ' -> ', sess, ': ', triggername{tt}]);
         else
              if ~isempty(spikedbfile)
                 [avgnow, xbinnow] = findotherdbavgs(pinfo, data, dataindthere, fdir, sess, triggertimenow, parm.binsize, parm.maxlag, parm.isavg);
              else
                 fname = eeg.general.eegfile{i};
                 if ~strcmp(eeg.parm.band{i}, 'cluster0')
                    [timestamp, dat, gain, fs] = ReadEEGFile(fname, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
                 else
                    dat = eegdata.cluster0.cnt{i}; fs = 1/eeg.parm.cl0Binsize(i);  timestamp = eegdata.cluster0.timepoint{i};
                 end
                 [avgnow, xbinnow] = findtriggeredaverage(dat, timestamp, fs, parm.binsize, parm.maxlag, triggertimenow, parm.isavg);
              end
              if parm.smoothbin>=3 avgnow = smoothavgs(avgnow, parm.smoothbin, parm.smoothtimes); end 
              if ~isempty(avgnow)
                 if isempty(avg)
                     avg = [avg; avgnow]; xbin = xbinnow;
                 else
                     if (numel(xbinnow) == numel(xbin))
                         avg = [avg; avgnow]; 
                     else
                         disp(['-------> warning: trace does not match with the previous one: ', fdir, ' -> ', sess]);
                     end
                 end
              end
         end
    end
    [mm,nn] = size(avg);
    if (mm>1)
       %se = std(avg)/sqrt(mm); avg = mean(avg);
       %strnow{1} = strcat(grpname, '_', evkeyword, '_', evkeytype);
       %strnow{2} = ['Number of triggers: ' num2str(mm)];
       %plottrigavg(avg, xbin, se, triggername{tt}, strnow);
       recarea = []; 
       allgrpname = []; for i=1:numel(grpname) allgrpname = [allgrpname '+' grpname{i}]; end
       plotaverageandmore(xbin, avg, recarea, allgrpname, seleventname, parm);
    else
       disp('-------> nothing to average');
    end
end
end
% function avgnow = smoothavgs(avgnow, smoothbin)%%%gaussian smoothing
% [mm,nn] = size(avgnow);
% for i = 1:mm
%     A = avgnow(i,:)'; A = OneDGaussianSimpleSmooth(A, smoothbin); avgnow(i,:) = A'; 
% end
function crr = smoothavgs(crr, smoothbin, nav) %%%moving averages nav times
for k = 1:nav
    crr = smoothnow(crr, smoothbin);
end
function crr = smoothnow(crr, smoothbin)
[mm,nn] = size(crr);
halfbin = floor(smoothbin/2);
for j=1:mm
    crrnow = crr(j,:); crrnow = movingaveragessmooth(crrnow, halfbin); crr(j,:) = crrnow;
end
function crrthen = movingaveragessmooth(crrnow, halfbin)
crrthen= zeros(size(crrnow));
for i = 1:numel(crrnow)
    kkk = [i-halfbin:i+halfbin]; kkk = kkk( (kkk>=1) & (kkk<=numel(crrnow)) );
    crrok = crrnow(kkk); crrthen(i) = mean(crrok(~isnan(crrok)));
end

function plottriggeredincidence(hfield, eeg, eegdata, dateind, tmark, plotparm, grpname)
catname = fieldnames(eeg); ntimevar = 0; triggername = []; triggerstr = []; ok = 1; nonevent = 0;
%%% selected incidence (time variables) triggered by selected events
for (i = 1:numel(catname)) %%%get all time variables to compute incidence
    subfield = fieldnames(eeg.(catname{i}));
    fieldselection = getappdata(hfield(i), 'selection');
    for (j = 1:numel(subfield))
        if (fieldselection(j) == 1) %%if a subfield is selected
            varnow = eeg.(catname{i}).(subfield{j});
            if (~strcmp(subfield{j}, 'eventname')) || (~strcmp(catname{i}, 'general'))
                for (k = 1:numel(varnow))
                    if isnumeric(varnow{k}) && (~isempty(varnow{k}))
                        ntimevar = ntimevar + 1; incidencename{ntimevar} = subfield{j}; nonevent = 1;
                        incidencetime{ntimevar} = varnow; incidencestr{ntimevar} = []; break;
                    end
                end
            end
        end
    end
end
varnow = eeg.general.eventname; %%% select an event file to trigger
unievent = [];
for (k = 1:numel(varnow))
     if ~isempty(varnow{k})
        unievent = unique(union(unievent, lower(varnow{k})));
     end
end
[sel,ok] = listdlg('ListString', unievent, 'SelectionMode', 'single', 'PromptString', 'Select an event to average');
if (ok)
   triggername = unievent{sel}; %this is just an event name, no time data here yet
   lststr = {'Start time'; 'End time'; 'Reference time'};
   [sel,ok] = listdlg('ListString', lststr, 'SelectionMode', 'single', 'PromptString', 'Select which time to average');
   if ok
       triggerstr = lststr{sel};
   end
end

%%%determine average parameters
if ok
  if ntimevar == 0
    disp('----> no time variables selected'); ok = 0;
  else
    p = {'Bin size (s)'; 'Max lag (s)'; 'Smooth bins (bins, 0=no smoothing)'; 'Smooth times'};
    d = {'0.05'; '10'; '4'; '5'};
    II = inputdlg(p, 'Average parameters', 4, d, 'on'); %%%resizable window
    if ~isempty(II)
        binsize = str2num(II{1}); maxlag = str2num(II{2}); smoothbin = str2num(II{3}); smoothtimes = str2num(II{4});
    else
        ok = 0;
    end
  end
end
%%%%determine what to average
if ok
for (tt = 1:ntimevar)
    incidence = cell(0,1); maxlagbin = round(maxlag/binsize); xbin = (-maxlagbin:maxlagbin)*binsize; 
    for (j = 1:numel(dateind))
         i = dateind(j);
         fdir = eeg.general.finaldir{i}; sess = eeg.general.sessname{i}; 
         evName = eeg.general.eventname{i}; evType = eeg.parm.eventtype{i}; evT = eegdata.event.eventtimes{i}; 
         triggertimenow = matchEventtimeout(evT, evType, evName, triggername, triggerstr);
            %disp(triggertimenow(1:5)');
         incidencetimenow = matchtimeout(eeg, eegdata, fdir, sess, incidencetime{tt});
            %disp(incidencetimenow(1:5));
         if isempty(triggertimenow)
             disp(['-------> warning: trigger times not found in ' fdir ' -> ' sess ': ' triggername]);
         else
              incidencenow = findtriggeredincidence(incidencetimenow, xbin, binsize, triggertimenow);   
              if ~isempty(incidencenow)
                 incidence = [incidence; incidencenow]; 
              end
         end
    end
    [avg, se] = findincidenceavgs(incidence, xbin, binsize, smoothbin, smoothtimes);
    if ~isempty(avg)  
       strnow{1} = strcat(grpname, '_', triggername);
       strnow{2} = ['Number of triggers: ' num2str(numel(incidence))];
       %plottrigavg(avg, xbin, se, incidencename{tt}, strnow);
       plottrigincidence(avg, se, incidence, xbin, incidencename{tt}, strnow);
    else
       disp(['--------> nothing to average: ' incidencename{tt}])
    end
end
end
function incidence = findtriggeredincidence(ripT, xbin, binsize, triggertime)
incidence = cell(numel(triggertime),1);
for i=1:numel(triggertime)
    tnow = ripT-triggertime(i);
    ttt = find( (tnow>=min(xbin)-binsize/2) & (tnow<=max(xbin)+binsize/2) ); 
    incidence{i} = sort(tnow(ttt));
end
function [avg, se] = findincidenceavgs(incidence, xbin, binsize, smoothbin, smoothtimes)
incicnt = zeros(numel(incidence), numel(xbin));
for i = 1:numel(incidence)
    incicnt(i,:) = hist(incidence{i}, xbin); %%%center bin times
end
incicnt = incicnt ./ binsize; %%%% now count is changed to rate
if smoothbin>=3 incicnt = smoothavgs(incicnt, smoothbin, smoothtimes); end  
if numel(incidence)>1
   avg = mean(incicnt); se = std(incicnt)/sqrt(numel(incidence));
else
   avg = incicnt; se = NaN*ones(1, numel(xbin));
end
function plottrigincidence(avg, se, ripT, X, incidencename, strnow)
f1name = ['Triggerd averages - ' incidencename];
hf1 = figure('Name', f1name); ninc = numel(ripT);
ha1 = axes('Parent', hf1, 'NextPlot', 'add', 'Position', [0.1 0.4 0.85 0.68], 'XLim', [min(X) max(X)], 'YLim', [0 ninc+1]);
for (i = 1:ninc)
    cl = [0 0 0];
    line(ripT{i}, i*ones(1, numel(ripT{i})), 'Parent', ha1, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'MarkerFaceColor', cl, 'MarkerEdgeColor', cl);
end
ylabel('Trigger event number');
ha2 = axes('Parent', hf1, 'NextPlot', 'add', 'Position', [0.1 0.05 0.85 0.27], 'XLim', [min(X) max(X)]);
line(X, avg, 'Parent', ha2, 'Color', [1 0 0], 'LineWidth', 2);
line(X, avg+se, 'Parent', ha2, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5); 
line(X, avg-se, 'Parent', ha2, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
xlabel('Trigger time(s)'); ylabel('Trigger average (/s)');
text('Interpreter', 'none', 'Parent', ha2, 'String', strnow{1}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
text('Interpreter', 'none', 'Parent', ha2, 'String', strnow{2}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.86]);

function [avgout, xbinnow] = findotherdbavgs(pinfo, data, dataind, fdir, sessname, triggertimenow, binsize, maxlag, isavg)
avgout = []; xbinnow = [];
ncell = numel(dataind); sel = strcmp(pinfo.general.finaldir(dataind), fdir); 
for (i = 1:ncell)
    if isempty(find(strcmp(pinfo.general.sessionname{dataind(i)}, sessname))) sel(i) = 0; end
end
nsel = numel(find(sel==1));
if nsel > 0
   spiketime = data.spike.spiketime(dataind); spiketime = spiketime(sel==1);
   for (i = 1:numel(spiketime)) spiketime{i} = spiketime{i}'; end
   spiketime = cell2mat(spiketime);
   maxlagbin = round(maxlag/binsize); xbinnow = (-maxlagbin:maxlagbin)*binsize;
   ntime = numel(triggertimenow); nbin = numel(xbinnow); avgnow = zeros(ntime, nbin);
   for (i = 1:ntime)
       bintimes = xbinnow + triggertimenow(i);
       spiketimenow = spiketime( (spiketime>=min(bintimes)) & (spiketime<=max(bintimes)) );
       cnt = hist(spiketimenow, bintimes); %%%center bin times
       avgnow(i,:) = cnt/nsel/binsize;
   end
   avgout = avgnow;
   if isavg
      if ntime > 1
         avgout = mean(avgnow);
      end
   end
end

function plotPSD(eeg, eegdata, dateind, tmark, plotparm)
neeg = numel(dateind);
for (k = 1:neeg)
    i = dateind(k); %disp(i);
    x = eeg.spec.freq{i}; sess = eeg.general.eegfile{i};
    if (plotparm.evselect == 0)
        y{1} = eeg.spec.sessPSD{i}; tinfix{1} = 'Whole session';
        y{2} = eeg.spec.sessNormPSD{i}; tinfix{2} = 'Whole session normalized';
    else
        y{1} = eeg.spec.evtPSD{i}; tinfix{1} = 'Selected events';
        y{2} = eeg.spec.evtNormPSD{i}; tinfix{2} = 'Selected events normalized';
    end
    for (iik = 1:numel(y))
        if (plotparm.setlog >0) 
            y{iik} = log10(y{iik}); stry = strcat('Log10(', tmark, ')');
        else
            stry = tmark;
        end
        if (~isempty(y{iik}))
           hf = figure('Name', strcat(sess, '---', tmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = tinfix{iik};
           line(x, y{iik}, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]); xlabel('Frequency (Hz)'); ylabel(stry);
           text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
        else
           str = tinfix{iik}; 
           disp(['----------------> PSD not computed yet: ' strcat(sess, '---', str)]);
        end
    end
end

function plotsignaldist(eeg, eegdata, dateind, tmark, plotparm)
neeg = numel(dateind);
binvector = [plotparm.range(1):plotparm.bin:plotparm.range(2)]; 
for (k = 1:neeg)
    i = dateind(k); 
    if ~strcmp(eeg.parm.band{i}, 'cluster0')
       %%%read the eegfile data
       [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
       dat = dat-mean(dat); sess = eeg.general.eegfile{i}; 
    else
       dat = eegdata.cluster0.cnt{i}; sess = eeg.general.eegfile{i}; 
       timestamp = eegdata.cluster0.timepoint{i};
    end
    if (plotparm.evselect == 0)
        y{1} = dat; tinfix{1} = 'Whole session';
    else
        for (iik = 1:numel(eeg.general.eventname{i}))
            tinfix{iik} = eeg.general.eventname{i}{iik}; iii = [];
            startT = eegdata.event.eventtimes{i}{iik}.start; endT = eegdata.event.eventtimes{i}{iik}.ent;
            for (tt = 1:numel(startT))
                ii = find( (timestamp>=startT(tt)) & (timestamp<=endT(tt)) ); iii = union(iii, ii);
            end
            y{iik} = dat(iii);
        end
    end
    for (iik = 1:numel(y))
        if (plotparm.setlog == 1) | (plotparm.setlog == 4) 
            y{iik} = log10(abs(y{iik})); stry = strcat('Log10(', tmark, ')');
        else
            stry = tmark;
        end
        hf = figure('Name', strcat(sess, '---', tmark)); hax = axes('Parent', hf, 'NextPlot', 'add');
        Y = histc(y{iik}, binvector); bar(hax, binvector, Y, 'EdgeColor', [0 0 1], 'FaceColor', [0 0 1]); 
        xlabel(stry); ylabel('count');
        text('Interpreter', 'none', 'Parent', hax, 'String', tinfix{iik}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
        strnow = strcat('nn=', num2str(numel(y{iik})), ';mean=', num2str(mean(y{iik})), ';std=', num2str(std(y{iik})), 'mV');
        text('Interpreter', 'none', 'Parent', hax, 'String', strnow, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92]);
    end
end

function plotPowerSpectrogram(eeg, eegdata, dateind, tmark, plotparm)
smooth = 0; Tsig = 0; NTsig = 0; Fsig = 0; NFsig = 0; ok = 1; 
showtimepoint = 100;  %default number of winodws to plot (1 window ~= 1 second as defined in parm.specWinShift)
SS = questdlg(['Smoothing the spectrogram?']);
if strcmp(SS, 'Yes')
   smooth = 1; 
   p = {'Time sigma (bin)'; 'Freq sigma (bin)'}; 
   d = {'1'; '1'};
   II = inputdlg(p, 'Smoothing parameters', 4, d, 'on'); %%%resizable window
   if ~isempty(II)
      Tsig = str2num(II{1}); NTsig = 5; Fsig = str2num(II{2}); NFsig = 5;
   else
      ok = 0;
   end
elseif strcmp(SS, 'Cancel')
   ok = 0;
end


if ok 
   p = {'Min freq (Hz)'; 'Max freq (Hz):'}; 
   d = {'1'; '400'};
   II = inputdlg(p, 'Frequency range', 2, d, 'on'); %%%resizable window
   if ~isempty(II)
      Fmin = str2num(II{1}); Fmax = str2num(II{2}); 
   else
      ok = 0;
   end
elseif strcmp(SS, 'Cancel')
   ok = 0;
end

if ok
if (~isempty(dateind))
    [MCroot, MCname, DAname, DEname] = CurrentVersion;
    hmain = figure('Name', DAname, 'NumberTitle', 'off', 'NextPlot', 'add', 'MenuBar', 'figure', 'Unit', 'normalized',...
       'OuterPosition', [0.204 0.425 0.779 0.576], 'Position', [0.207 0.429 0.773 0.496]);
    hf = uimenu(hmain, 'Label','Add Plot');
    uimenu(hf,'Label','Add New Spike','Callback','DataAnimator_AddSpike_Callback');
    uimenu(hf,'Label','Add New Position','Callback','DataAnimator_AddPosition_Callback');
    uimenu(hf,'Label','Add New EEG','Callback','DataAnimator_AddEEG_Callback');
    uimenu(hf,'Label','Add New Powergram/ratio','Callback','DataAnimator_AddPower_Callback', 'Separator', 'on');
    uimenu(hf,'Label','Add New PSD','Callback','DataAnimator_AddPSD_Callback');
    uimenu(hf,'Label','Add New SleepClass','Callback','DataAnimator_AddSleepClass_Callback');
    ndatafile = 0; displaysetting = []; haxes = {}; data = {};
    setappdata(hmain, 'ndatafile', ndatafile); setappdata(hmain, 'displaysetting', displaysetting);
    setappdata(hmain, 'haxes', haxes); setappdata(hmain, 'data', data);
    ndatafile = 0;
    ndisplay = 0; 
end    

for (j = 1:numel(dateind))
    k = dateind(j); filename = eeg.general.eegfile{k}; freqY = eeg.spec.freq{k}; timewinX = eeg.spec.timewin{k};
    iii = find( (freqY>=Fmin) & (freqY<=Fmax) ); freqY = freqY(iii);
    winpower = eegdata.spec.winpower{k}(iii,:); %psd[fy timewin]
    if (plotparm.setlog>0) winpower = log10(abs(winpower)); end%%%this change the plot into log-scale
    [mm,nn] = size(winpower);
    for (im = 1:mm)
        for (in = 1:nn)
            if (winpower(im,in) == -Inf) winpower(im,in) = NaN; end
        end
    end
 
    ntimepoint = numel(timewinX); nfreqpoint = numel(freqY);
    showtimepoint = min([showtimepoint ntimepoint]);
if (ntimepoint ~= 0)
    ndatafile = ndatafile + 1;
    displaysetting{ndatafile}{1} = ndatafile;
    displaysetting{ndatafile}{2} = filename;
    displaysetting{ndatafile}{3} = 'd'; % data type
    displaysetting{ndatafile}{4} = [0 0]; %default plot, continuous mode
    displaysetting{ndatafile}{5} = [ntimepoint nfreqpoint 0 ndisplay showtimepoint 2 nfreqpoint+1]; % 
    displaysetting{ndatafile}{6} = [0.5 0.5 0.5]; %default backgrnd
    displaysetting{ndatafile}{7} = [0 0 1]; %default line color
    displaysetting{ndatafile}{8} = 0.5; %default dot line width
    displaysetting{ndatafile}{9} = 'none'; %default marker
    set(0, 'ShowHiddenHandles', 'on');
    H = get(hmain, 'Children');
    TypeH = get(H, 'Type');
    for (i = 1:size(TypeH))
        if (strcmp(TypeH{i}, 'uicontrol')) delete(H(i)); end
    end
    delete(findobj(hmain, 'Tag', 'colorbar'));
    set(0, 'ShowHiddenHandles', 'off');
    
    setappdata(hmain, 'ndatafile', ndatafile);
    setappdata(hmain, 'displaysetting', displaysetting);  

    %%re-allocacte spaces for all the plots
    [axposvector, uiposvector] = DataAnimator_PlotSpaceAuto(ndatafile,0.8);
    %%re-plot all axes and uicontrols
    DataAnimator_Replot(haxes, ndatafile-1, axposvector, uiposvector);
    
    %%add new pSD plot
    displaysetting = getappdata(hmain, 'displaysetting');
    data = getappdata(hmain, 'data');

%%read settings
idnum = ndatafile;
fname = displaysetting{idnum}{2};
ntimepoint = displaysetting{idnum}{5}(1); %total number of time points
nfreqpoint = displaysetting{idnum}{5}(2); %total number of frequency points
finitpos = displaysetting{idnum}{5}(3); %binary data start postion after the header
ndisplay = displaysetting{idnum}{5}(4); %start display number
showtimepoint = displaysetting{idnum}{5}(5); %display time period in time points
startfreq = displaysetting{idnum}{5}(6); %starting frequency index
endfreq = displaysetting{idnum}{5}(7); %end frequency index
maxdisplay = ceil(ntimepoint/showtimepoint);

data{idnum}(1,1) = 0; %first element of the data matrix
data{idnum}(1,2:ntimepoint+1) = timewinX; %read time stamps to first row
data{idnum}(2:nfreqpoint+1,1) = freqY; %read frequency points to first column
data{idnum}(2:nfreqpoint+1, 2:ntimepoint+1) = winpower; 

minfreq = min(data{idnum}(2:nfreqpoint+1,1)); 
maxfreq = max(data{idnum}(2:nfreqpoint+1,1)); 

freqindex = find((data{idnum}(:,1)>=minfreq)&(data{idnum}(:,1)<=maxfreq));
startfreq = max([min(freqindex) 2]); 
endfreq = max(freqindex);
displaysetting{idnum}{5}(6) = startfreq;
displaysetting{idnum}{5}(7) = endfreq;

A = data{idnum}(startfreq:endfreq, 2:ntimepoint+1); 
% %%%%smoothing the power spectrum
if smooth == 1
   disp('----------> smoothing the spectrogram ......');  
   A = TwoDSmooth_separate(A, Tsig, NTsig, Fsig, NFsig);
   %B = TwoDSmooth_fast(A, Tsig, Fsig, NTsig, NFsig);
   %occu = ones(size(A)); B = TwoDSmooth_new(A, occu, Tsig, Fsig, NTsig, NFsig); 
   data{idnum}(startfreq:endfreq, 2:ntimepoint+1) = A; 
end
minvalue = 0; maxvalue = 1;
if ~isempty(A)
   [mm,nn] = size(A); BB = reshape(A, mm*nn, 1); BB = BB(~isnan(BB));
   %ma = mean(BB); ss = std(BB);
   %minvalue = ma-5*ss; maxvalue =ma+5*ss;
   minvalue = prctile(BB, 2); maxvalue = prctile(BB, 98);
end
%minvalue = min(min(data{idnum}(startfreq:endfreq, 2:ntimepoint+1)));
%maxvalue = max(max(data{idnum}(startfreq:endfreq, 2:ntimepoint+1)));
if (isempty(minvalue)) || (isempty(maxvalue)) || (minvalue == maxvalue) || isnan(minvalue) || isnan(maxvalue)
    minvalue = 0; maxvalue = 1; disp('----------> spectrogram data seem not meaningful!'); 
end
    
posvec = axposvector{idnum};
haxes{idnum} = axes('Parent', hmain, 'Units', 'normalized', 'Position', posvec, 'FontSize', 8, 'Visible', 'off', 'Color', [1 1 1],...
    'CLim', [minvalue maxvalue], 'CLimMode', 'manual', 'XTickMode', 'manual', 'XTickLabelMode', 'manual');
displaysetting{idnum}{10} = [minvalue maxvalue];
%displaysetting{idnum}{10} = [-3 6];
ndatapoint=size(data{idnum},1);
   disp('----------> display the spectrogram ......'); 
   starttimepoint = ndisplay*showtimepoint+2; % PSD time starts from second column
   endtimepoint = min(starttimepoint+showtimepoint-1, ntimepoint);
   hp = pcolor(haxes{idnum}, data{idnum}(1,starttimepoint:endtimepoint), data{idnum}(startfreq:endfreq,1), data{idnum}(startfreq:endfreq, starttimepoint:endtimepoint));
   xlabel('Time (s)'); ylabel('Frequency (Hz)');
   set(hp, 'EdgeColor', 'none', 'LineStyle', 'none');
   colormap(jet);  
   hc = colorbar('vert', 'peer', haxes{idnum}); %, 'CLim', [minvalue maxvalue]); %, 'CLimMode', 'manual');
   set(hc, 'Tag', 'colorbar');
   
   tickdis = (data{idnum}(1,endtimepoint)-data{idnum}(1,starttimepoint))/5; %10 ticks
   labelnumber = data{idnum}(1,starttimepoint):tickdis:data{idnum}(1,endtimepoint); 
   for (i = 1:numel(labelnumber))
       labelstr{i} = num2str(labelnumber(i), 8);
   end
   set(haxes{idnum}, 'XTickMode', 'manual', 'XTick', labelnumber, 'XTickLabelMode', 'manual', 'XTickLabel', labelstr);
   set(gca, 'FontSize', 8);
   
   titletext = strcat(num2str(idnum), ' = ', fname,' -PSD');
   htitle=text('String', titletext, 'Interpreter', 'none', 'Parent', haxes{idnum}, 'Tag', 'title',...
              'Position', [0.01 0.95], 'Units', 'normalized', 'Color', [1 0 0], 'FontSize', 8); % add a title
   %displaysetting{idnum}{10} = [0 50]; %
   %displaysetting{idnum}{10} = get(haxes, 'CLim');
   set(haxes{idnum}, 'CLim', displaysetting{idnum}{10});
   set(haxes{idnum}, 'CLimMode', 'manual');

%%plot a control bar with the axes (%%%main figure (hmain))
barpos = [posvec(1) posvec(2)-0.20*posvec(4)/0.75 posvec(3) posvec(4)*0.06];
hbar = uicontrol('Style', 'slider', 'Parent', hmain, 'Units', 'normalized', 'Position', barpos,...
                 'Min', 0, 'Max', maxdisplay, 'SliderStep', [1/maxdisplay, 10/maxdisplay], 'Callback', 'DataAnimator_PSDBar_Callback',...
                 'Tag', strcat('slider_',num2str(idnum)));
set(hbar, 'Value', ndisplay); %%initial bar position
setappdata(hmain, 'displaysetting', displaysetting);
setappdata(hmain, 'data', data);
setappdata(haxes{idnum}, 'hcolorbar', hc); setappdata(haxes{idnum}, 'hpcolor', hp);

    %%add new control panel
    DataAnimator_PSDControl(hmain, haxes{ndatafile}, ndatafile, uiposvector{ndatafile});
    setappdata(hmain, 'haxes', haxes);
end
end
end

function showalldata(eeg, eegdata, dateind, tmark)
if (~isempty(dateind))
    [MCroot, MCname, DAname, DEname] = CurrentVersion;
    hmain = figure('Name', DAname, 'NumberTitle', 'off', 'NextPlot', 'add', 'MenuBar', 'figure', 'Unit', 'normalized',...
       'OuterPosition', [0.204 0.425 0.779 0.576], 'Position', [0.207 0.429 0.773 0.496]);
    hf = uimenu(hmain, 'Label','Add Plot');
    uimenu(hf,'Label','Add New Spike','Callback','DataAnimator_AddSpike_Callback');
    uimenu(hf,'Label','Add New Position','Callback','DataAnimator_AddPosition_Callback');
    uimenu(hf,'Label','Add New EEG','Callback','DataAnimator_AddEEG_Callback');
    uimenu(hf,'Label','Add New Powergram/ratio','Callback','DataAnimator_AddPower_Callback', 'Separator', 'on');
    uimenu(hf,'Label','Add New PSD','Callback','DataAnimator_AddPSD_Callback');
    uimenu(hf,'Label','Add New SleepClass','Callback','DataAnimator_AddSleepClass_Callback');
    ndatafile = 0; displaysetting = []; haxes = {}; data = {};
    setappdata(hmain, 'ndatafile', ndatafile); setappdata(hmain, 'displaysetting', displaysetting);
    setappdata(hmain, 'haxes', haxes); setappdata(hmain, 'data', data);
    nread = 0; nreadbuffer = 4;  %default number of buffers being read and plot (4 buffs ~= 1 second data)
    defaultbuffsize = 512; %%%this is the default buffersize -EEG data will use actual buffersize
    defaultfreq = 2034; %%%this is the default sampling frequency - EEG data will use actual vale
end
for (j = 1:numel(dateind))
    k = dateind(j);
    if ~strcmp(eeg.parm.band{k}, 'cluster0') %%%if EEG traces
      filename = eeg.general.eegfile{k}; if exist(filename, 'file')~=2 filename = strrep(filename, 'I:', 'J:'); end
      displaysetting{j}{2} = filename; 
      gain = eeg.general.eeggain{k}; %%gain already in mV
      freq = eeg.general.freq{k}; timeunit = eeg.parm.timeunit(k); buffersize = eeg.parm.buffersize(k);
      fid = fopen(filename);  %first open a file for read
      while 1               %find where header ends
          tline = fgets(fid);
          if (contains(tline, 'ADMaxValue')) 
             [str, ent] = strtok(tline, ' '); ent = ent(2:numel(ent));  maxADvalue = str2num(ent); %maximum AD value in the data
          end
          if (strncmpi(tline, '%%ENDHEADER', 8)), break, end
      end
      finitpos = ftell(fid);       %header end position, initial data reading position
      ok = fclose(fid);
      %get total number of buffers in the file
      maxread = EEG_ReadTotalBuffer(filename, nreadbuffer, finitpos, buffersize);
      %%%%assign parameters to display setting file and variables
      displaysetting{j}{1} = j;  % id number
      displaysetting{j}{3} = 'e'; % data type
      displaysetting{j}{4} = [0 0]; %default plot, continuous mode
      displaysetting{j}{5} = [nread nreadbuffer finitpos buffersize gain maxread freq maxADvalue];
      displaysetting{j}{6} = [0.5 0.5 0.5]; %default backgrnd
      displaysetting{j}{7} = [0 0 1]; %default line color
      displaysetting{j}{8} = 0.5; %default dot line width
      displaysetting{j}{9} = 'none'; %default marker
    
    else %%%if cluster0 power data
      filename = eeg.general.eegfile{k};
      displaysetting{j}{2} = filename;
      freq = NaN; finitpos = NaN;
      wstime = eeg.parm.cl0Binsize(k); %%%window size time
      shifttime = wstime; %%window shift time
      ntimestamp = numel(eegdata.cluster0.timepoint{k});
      %%%%assign parameters to display setting file and variables
      ndisplay = 0; %default start display: first display
      showtimepoint = (nreadbuffer*defaultbuffsize/defaultfreq)/wstime; %computed to align with EEG display time axes
      displaysetting{j}{1} = j;  % id number
      displaysetting{j}{3} = 'w'; % data type
      displaysetting{j}{4} = [0 0]; %default plot, continuous mode
      displaysetting{j}{5} = [wstime shifttime finitpos ndisplay showtimepoint ntimestamp]; % 
      displaysetting{j}{6} = [0.5 0.5 0.5]; %default backgrnd
      displaysetting{j}{7} = [0 0 1]; %default line color
      displaysetting{j}{8} = 0.5; %default dot line width
      displaysetting{j}{9} = 'none'; %default marker
    end     
        
      %%%%%get the band-specific data
      startT{j} = []; endT{j} = []; peakT{j} = []; startA{j} = []; endA{j} = []; peakA{j} = [];
      %switch eeg.parm.band{k}
      %case 'ripple' 
           if isfield(eeg, 'ripple')
               if isfield(eeg.ripple, 'sessStartT')
                  startT{j} = eeg.ripple.sessStartT{k}; endT{j} = eeg.ripple.sessEndT{k}; peakT{j} = eeg.ripple.sessPeakT{k};
                  nrip = numel(eeg.ripple.sessStartT{k}); peakA{j} = eeg.ripple.sessAmp{k};
               elseif isfield(eeg.ripple, 'StartT')
                  startT{j} = eeg.ripple.StartT{k}; endT{j} = eeg.ripple.EndT{k}; peakT{j} = eeg.ripple.PeakT{k};
                  nrip = numel(eeg.ripple.StartT{k}); peakA{j} = eeg.ripple.Amp{k};
               end
               if isfield(eeg.parm, 'rippPeakMode')
                   if (~isempty(eeg.parm.rippPeakMode{k})) && contains(eeg.parm.rippPeakMode{k}, 'neg')
                       startA{j} = -1*eeg.ripple.startThreshold{k}*ones(1, nrip); endA{j} = -1*eeg.ripple.startThreshold{k}*ones(1, nrip); 
                   elseif (~isempty(eeg.parm.rippPeakMode{k}))
                       startA{j} = eeg.ripple.startThreshold{k}*ones(1, nrip); endA{j} = eeg.ripple.startThreshold{k}*ones(1, nrip); 
                   end
               end
           end
      %case 'theta'
           if (contains(filename, 'smooth')) && (isfield(eeg, 'theta')) && (strcmp(eeg.parm.band{k}, 'theta'))
              startT{j} = eeg.theta.maxTime{k}; endT{j} = eeg.theta.minTime{k}; 
              startA{j} = eeg.theta.maxAmp{k}; endA{j} = eeg.theta.minAmp{k}; 
           end
      %case 'spindle'
           if (isfield(eeg, 'spindle')) && (strcmp(eeg.parm.band{k}, 'spindle'))
              startT{j} = eeg.spindle.maxTime{k}; endT{j} = eeg.spindle.minTime{k}; 
              startA{j} = eeg.spindle.maxAmp{k}; endA{j} = eeg.spindle.minAmp{k}; 
           end
      %case 'hvs'
           if (isfield(eeg, 'hvs')) && (strcmp(eeg.parm.band{k}, 'hvs'))
              startT{j} = eeg.hvs.maxTime{k}; endT{j} = eeg.hvs.minTime{k}; 
              startA{j} = eeg.hvs.maxAmp{k}; endA{j} = eeg.hvs.minAmp{k}; 
           elseif (isfield(eeg, 'hvs')) && (contains(eeg.parm.band{k}, 'gamma'))
              startT{j} = eeg.hvs.sessStartT{k}; endT{j} = eeg.hvs.sessEndT{k}; nhvs = numel(eeg.hvs.sessStartT{k});
              peakT{j} = eeg.hvs.sessPeakT{k}; peakA{j} = eeg.hvs.sessAmp{k};  
              if isfield(eeg.parm, 'hvsPeakMode')
                   if (~isempty(eeg.parm.hvsPeakMode{k})) && contains(eeg.parm.hvsPeakMode{k}, 'neg')
                       startA{j} = -1*eeg.hvs.startThreshold{k}*ones(1, nhvs); endA{j} = -1*eeg.hvs.startThreshold{k}*ones(1, nhvs); 
                   elseif (~isempty(eeg.parm.hvsPeakMode{k}))
                       startA{j} = eeg.hvs.startThreshold{k}*ones(1, nhvs); endA{j} = eeg.hvs.startThreshold{k}*ones(1, nhvs); 
                   end
              end
           end
      %case 'cluster0'
           if (isfield(eeg, 'cluster0')) %&& (strcmp(eeg.parm.band{k}, 'cluster0'))
               if isfield(eeg, 'slowtheta')
                  ncl0 = numel(eeg.slowtheta.sessStartT{k}); 
                  startT{j} = eeg.slowtheta.sessStartT{k}; endT{j} = eeg.slowtheta.sessEndT{k}; %peakT{j} = eeg.slowtehta.sessPeakT{k}; 
                  startA{j} = eeg.slowtheta.startThreshold{k}*ones(1, ncl0); endA{j} = eeg.slowtheta.startThreshold{k}*ones(1, ncl0); 
                  %peakA{j} = eeg.slowtehta.sessAmp{k};
               else
                  ncl0 = numel(eeg.cluster0.sessStartT{k}); 
                  startT{j} = eeg.cluster0.sessStartT{k}; endT{j} = eeg.cluster0.sessEndT{k}; peakT{j} = eeg.cluster0.sessPeakT{k}; 
                  startA{j} = eeg.cluster0.startThreshold{k}*ones(1, ncl0); endA{j} = eeg.cluster0.startThreshold{k}*ones(1, ncl0); 
                  peakA{j} = eeg.cluster0.sessAmp{k};
               end
           end
      %end
end

set(0, 'ShowHiddenHandles', 'on');
H = get(hmain, 'Children');
TypeH = get(H, 'Type');
for (i = 1:size(TypeH))
     if (strcmp(TypeH{i}, 'uicontrol')) delete(H(i)); end
end
delete(findobj(hmain, 'Tag', 'colorbar'));
set(0, 'ShowHiddenHandles', 'off');

ndatafile = numel(dateind);
setappdata(hmain, 'ndatafile', ndatafile);
setappdata(hmain, 'displaysetting', displaysetting);  
%%re-allocacte spaces for all the plots
[axposvector, uiposvector] = DataAnimator_PlotSpaceAuto(ndatafile,0.8);

for (j = 1:numel(dateind))
    k = dateind(j);
    if ~strcmp(eeg.parm.band{k}, 'cluster0') %%%if EEG traces add new EEG plot/control
       [haxes{j}, hbar] = DataAnimator_EEGDisplay(j, hmain, axposvector{j});
       DataAnimator_EEGControl(hmain, haxes{j}, hbar, j, uiposvector{j});
    else %%%if EEG traces add new EEG plot/control
%         size(eegdata.cluster0.timepoint{k})
%         size(eegdata.cluster0.cnt{k})
       haxes{j} = DataAnimator_PowergramDisplay(j, hmain, axposvector{j},...
           eegdata.cluster0.timepoint{k}', eegdata.cluster0.cnt{k}');
       DataAnimator_PowergramControl(hmain, haxes{j}, j, uiposvector{j}); 
    end
    setappdata(haxes{j}, 'startT', startT{j}); setappdata(haxes{j}, 'endT', endT{j}); setappdata(haxes{j}, 'peakT', peakT{j}); 
    setappdata(haxes{j}, 'startA', startA{j}); setappdata(haxes{j}, 'endA', endA{j}); setappdata(haxes{j}, 'peakA', peakA{j});
    %%%%%%%%%%plot events
    TT = get(haxes{j}, 'XLim'); minT = TT(1); maxT = TT(2);  VV = get(haxes{j}, 'YLim'); minV = VV(1); maxV = VV(2);
    indstart = find( (startT{j}>=minT) & (startT{j}<=maxT) );
    indend = find( (endT{j}>=minT) & (endT{j}<=maxT) );
    indpeak = find( (peakT{j}>=minT) & (peakT{j}<=maxT) );
    if ~strcmp(eeg.parm.band{k}, 'cluster0') %%%if EEG traces add new EEG plot/control
        for (k = 1:numel(indstart))
            line(startT{j}(indstart(k)), startA{j}(indstart(k)), 'Parent', haxes{j}, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20,...
                  'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]);
        end
        for (k = 1:numel(indend))
            line(endT{j}(indend(k)), endA{j}(indend(k)), 'Parent', haxes{j}, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20,...
                  'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
        end    
    else
        for (k = 1:numel(indstart))
            line([startT{j}(indstart(k)) startT{j}(indstart(k))], [0.8*minV 0.8*maxV], 'Parent', haxes{j}, 'LineWidth', 1, 'Color', [0 1 0]);
        end
        for (k = 1:numel(indend))
            line([endT{j}(indend(k)) endT{j}(indend(k))], [0.8*minV 0.8*maxV], 'Parent', haxes{j}, 'LineWidth', 1, 'Color', [1 0 0]);
        end  
    end
    for (k = 1:numel(indpeak))
         line(peakT{j}(indpeak(k)), peakA{j}(indpeak(k)), 'Parent', haxes{j}, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30,...
                  'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
    end
end
setappdata(hmain, 'haxes', haxes);

function plotcrrnow(dat1, dat2, timestamp, fs, binsize, maxlag, tmark, st, et, eventname, f1name, f2name)
binnow = round(binsize*fs); %%%%need to re-sample the data, otherwise the computation is too slow
ind = 1:binnow:numel(dat1); dat1 = dat1(ind); dat2 = dat2(ind); timestamp = timestamp(ind);
maxlagbin = round(maxlag/binsize); lags = (-maxlagbin:maxlagbin)*binsize;
STind = ones(size(st)); ETind = ones(size(et));
for (i = 1:numel(st))
    iii = find( (timestamp>=st(i)) & (timestamp<=et(i)) );
    if ~isempty(iii)
       STind(i) = min(iii); ETind(i) = max(iii);
    end
end
%%%%%%For the function below, dat1 dat2 are row vectors with same length
[mm,nn] = size(dat1); if mm>1 dat1 = dat1'; end
[mm,nn] = size(dat2); if mm>1 dat2 = dat2'; end
ccc = Utilities_FindXCorrCoef_self(dat1, dat2, maxlagbin, STind, ETind);
hg = figure('Name', tmark); 
hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
xlabel ('time lag (s)'); ylabel('Correlation'); 
line(lags, ccc, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]);
text('Interpreter', 'none', 'Parent', hax, 'String', f1name, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.94]);
text('Interpreter', 'none', 'Parent', hax, 'String', f2name, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.91]);
text('Interpreter', 'none', 'Parent', hax, 'String', eventname, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.88]);

function plottrigavg(avg, xbin, se, tmark, grpname)
hg = figure('Name', tmark); 
hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
xlabel ('time lag (s)'); ylabel('Triggered average'); 
line(xbin, avg, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]);
line(xbin, avg-se, 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
line(xbin, avg+se, 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
for (i = 1:numel(grpname))
     text('Parent', hax, 'Interpreter', 'none', 'String', grpname{i}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
end

function triggertimenow = matchEventtimeout(evT, evType, evName, eventname, timepoint)
triggertimenow = []; %%%the key issue here is that evT (like track1_CW.evts) could be associated with multiple EEG traces that belong to the same session
for (i = 1:numel(evName))
    if strcmpi(eventname, evName{i})
      if ~isempty(evT{i}.start)
        if strcmp(timepoint, 'Start time')
           triggertimenow = [triggertimenow; evT{i}.start];
        elseif strcmp(timepoint, 'End time')
           triggertimenow = [triggertimenow; evT{i}.ent];
        elseif strcmp(timepoint, 'Reference time')
           triggertimenow = [triggertimenow; evT{i}.ref];
        end
      end
    end
end
function triggertimenow = matchtimeout(eeg, eegdata, fdir, sess, triggertime)
triggertimenow = [];  %%%the key issue here is that triggertime (like ripple start times) could be determined from multiple EEG traces that belong to the same session
iii = find(strcmp(eeg.general.finaldir, fdir) & strcmp(eeg.general.sessname, sess) ); %%%this is the session index
nm = numel(iii);
if nm>0
    triggertime = triggertime(iii); sel = zeros(1, nm);
    for (i = 1:nm)
        if (isnumeric(triggertime{i}) && (~isempty(triggertime{i}))) sel(i) = 1; end
    end
    ij = find(sel==1); %disp(iii(ij))
    if numel(ij) == 1
        triggertimenow = triggertime{ij};
    elseif numel(ij) > 1
        disp(['-------> warning: multiple matches in ', fdir, ' -> ', sess, '; use the first one']);
        triggertimenow = triggertime{ij(1)};
    else 
        disp(['-------> warning: no matches in ', fdir, ' -> ', sess]);
    end
end
function triggertime = filtertriggertime(triggertime,evkeyword, evkeytype, evName, evType, evT) 
if ~(isempty(evkeyword) && isempty(evkeytype) ) %%if only keywords/keytypes specified
evsel = ones(1,numel(evName));
if ~isempty(evkeyword)
   for (i = 1:numel(evName))
       if isempty(strfind(lower(evName{i}), lower(evkeyword))) evsel(i) = 0; end 
   end
end
if ~isempty(evkeytype)
   for (i = 1:numel(evName))
       if isempty(strfind(lower(evType{i}), lower(evkeytype))) evsel(i) = 0; end 
   end
end
evpos = find(evsel == 1); evT = evT(evpos);
iii = [];
for (ij = 1:numel(evT))
    startT = evT{ij}.start; entT = evT{ij}.ent;
    if ~isempty(startT)
        for (i = 1:numel(startT))
            iiok = find( (triggertime>=startT(i)) & (triggertime<=entT(i)) );
            iii = union(iii, iiok);
        end
    end
end
triggertime = triggertime(iii);
end

function [avgnow, xbinnow] = findtriggeredaverage(dat, timestamp, fs, binsize, maxlag, triggertimenow, isavg)
binnow = round(binsize*fs); %%%%need to re-sample the data, otherwise the computation is too slow
ind = 1:binnow:numel(dat); dat = dat(ind); timestamp = timestamp(ind);
maxlagbin = round(maxlag/binsize); xbinnow = (-maxlagbin:maxlagbin)*binsize; xbin = -maxlagbin:maxlagbin;
iii = find( (triggertimenow>timestamp(1)+maxlag) & (triggertimenow<timestamp(numel(timestamp))-maxlag) );
triggertimenow = triggertimenow(iii);
ntime = numel(triggertimenow); nbin = numel(xbinnow); avgnow = zeros(ntime, nbin);
for (i = 1:ntime)
    [~,iinow] = min(abs(timestamp-triggertimenow(i)));
    iii = xbin + iinow;
    avgnow(i,:) = dat(iii);
end
if isavg
   if ntime > 1
      avgnow = mean(avgnow);
   end
end
% for (i = 1:ntime)
%     iii = find( (timestamp>=triggertimenow(i)-maxlag) & (timestamp<=triggertimenow(i)+maxlag) );
%     numel(iii)
%     if (numel(iii)==nbin+1)
%         iii = iii(1:nbin);
%     elseif (numel(iii)==nbin-1)
%         iii = [iii; max(iii)];
%     end
%     size(avgnow)
%     size(dat(iii))
%     avgnow(i,:) = dat(iii);
% end

function X = plotaverageandmore(xbin, allcrr, dbfile, grpname, tag, parm)
[nsess, ~] = size(allcrr); X.allpeakR = NaN*ones(1, nsess); X.allpeakT = NaN*ones(1, nsess); 
X.alltroughR = NaN*ones(1, nsess); X.alltroughT = NaN*ones(1, nsess); basetime = parm.basetime; biastime = parm.biastime;
X.allbias = NaN*ones(1, nsess); X.allbase = NaN*ones(1, nsess); X.parm = parm;
nbin = numel(xbin); crr = NaN*ones(1, nbin); err = NaN*ones(1, nbin); nn = NaN*ones(1, nbin);
for (j = 1:nbin)
     iii = find(~isnan(allcrr(:,j)));
     %crr(j) = median(allcrr(iii,j)); 
     crr(j) = mean(allcrr(iii,j)); 
     nn(j) = numel(iii); 
     if (nn(j) > 0) err(j) = std(allcrr(iii,j))/sqrt(nn(j)); end
end
mn = mean(nn); sn = std(nn)/sqrt(nbin);
f2name = 'Triggered averages';
hg = figure('Name', f2name); 
hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
xlabel ('Trigger time (s)'); ylabel('Average'); 
line(xbin, crr+err, 'Parent', hax, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
line(xbin, crr-err, 'Parent', hax, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
% for i =1:nsess
%    line(xbin, allcrr(i,:), 'Parent', hax, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]); 
% end
line(xbin, crr, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]);
%Drawerrorupdown(xbin, crr, crr+err, crr-err, hax, [1 0 0]);
text('Interpreter', 'none', 'Parent', hax, 'String', tag, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.97]);
text('Interpreter', 'none', 'Parent', hax, 'String', dbfile, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.94]);
text('Interpreter', 'none', 'Parent', hax, 'String', grpname, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.91]);
str = ['N (mean, se): ', num2str(mn), ', ', num2str(sn)];
text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.88]);
% if (~isempty(dd)) saveas(hf2, fullfile(dd, strcat(f2name, '.fig'))); end
%%%%%% test peak crr significance (and from zero)
valnow = NaN*ones(1, nbin*nsess); grp = NaN*ones(1, nbin*nsess);
for (i = 1:nsess)
     for (j = 1:nbin)
          valnow(j+(i-1)*nsess) = allcrr(i,j); grp(j+(i-1)*nsess) = j;
     end
end
iii = find(~isnan(valnow)); valnow = valnow(iii); grp = grp(iii);
[ppp, table, stats] = anova1(valnow, grp, 'off'); 
str = ['anova1 p=', num2str(ppp, '%10.5e')];
text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.84]);
[peakr,iii] = max(crr);[troughr,jjj] = min(crr); 
[~,peakp] = ttest(allcrr(:,iii)); [~,troughp] = ttest(allcrr(:,jjj));
str1 = ['peak: r=', num2str(peakr), '; se=', num2str(err(iii)), '; t=', num2str(xbin(iii)), '; p(from 0)=', num2str(peakp, '%10.5e')];
str2 = ['trough: r=', num2str(troughr), '; se=', num2str(err(jjj)),  '; t=', num2str(xbin(jjj)), '; p(from 0)= ', num2str(troughp, '%10.5e')];
text('Interpreter', 'none', 'Parent', hax, 'String', str1, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.80]);
text('Interpreter', 'none', 'Parent', hax, 'String', str2, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.76]);
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
   str3 = ['bias: mean=', num2str(mm), '; se=', num2str(se), '; n=', num2str(nn), '; p (from 0) = ', num2str(biasp, '%10.5e')];
   text('Interpreter', 'none', 'Parent', hax, 'String', str3, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.72]);
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
   str4 = ['baseline: mean=', num2str(mm), '; se=', num2str(se), '; n=', num2str(nn), '; p (from peak) = ', num2str(basep, '%10.5e')];
   text('Interpreter', 'none', 'Parent', hax, 'String', str4, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.68]);
end
%%%save data to figure
setappdata(hg, 'X', X);



