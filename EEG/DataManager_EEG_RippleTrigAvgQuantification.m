function DataManager_EEG_RippleTrigAvgQuantification
%%quantify the complexity ripple-triggered spike averages
hf = gcbf; eeg = getappdata(hf, 'eeg'); eegdata = getappdata(hf, 'eegdata');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
hfield = getappdata(hf, 'hfield'); plotparm = getappdata(hf, 'plotparm');
hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); 
dateind = find(spikeselection == 1); %%%%%dateind here is actually eegind (eeg IDs)
plotparm = getappdata(hf, 'plotparm');
grpind = find(groupselection == 1); grpname = eegdata.grouplist.groupname(grpind);
hfield = getappdata(hf, 'hfield');
plottriggeredaverages(hfield, eeg, eegdata, dateind, 'triggeredaverage', plotparm, grpname);
disp('************************');

function plottriggeredaverages(hfield, eeg, eegdata, dateind, tmark, plotparm, grpname)
catname = fieldnames(eeg); ntimevar = 0; triggertime = []; triggername = []; ok = 1;
for (i = 1:numel(catname))
    subfield = fieldnames(eeg.(catname{i}));
    fieldselection = getappdata(hfield(i), 'selection');
    for (j = 1:numel(subfield))
        if (fieldselection(j) == 1) %%if a subfield is selected, delete it
            varnow = eeg.(catname{i}).(subfield{j});
            for (k = 1:numel(varnow))
                if isnumeric(varnow{k})
                    ntimevar = ntimevar + 1; 
                    triggertime{ntimevar} = varnow; triggername{ntimevar} = subfield{j}; break;
                end
            end
        end
    end
end
if ntimevar == 0
    disp('----> no time variables selected'); ok = 0;
else
    p = {'Bin size (s)'; 'Max lag (s)'};
    d = {'0.01'; '0.1'};
    II = inputdlg(p, 'Average parameters', 3, d, 'on'); %%%resizable window
    if ~isempty(II)
        binsize = str2num(II{1}); maxlag = str2num(II{2});
    else
        ok = 0;
    end
end
if ok
    input = inputdlg({'Event keyword';'Event type'; 'Plot single average over session/event?'}, 'Event slelection parms', 3, {'sws'; 'sws'; 'no'}); 
    if (~isempty(input))
        evkeyword = input{1}; evkeytype = input{2}; 
        if strncmpi(input{3}, 'yes', 1)
            isavg = 1;
        else
            isavg = 0;
        end
    else
        ok = 0;
    end
end
if ok
   spikedbfile = [];
   button = questdlg('Average a spike database?');
   if (strcmp(button, 'Yes'))
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
         triggertimenow = matchtimeout(eeg, eegdata, fdir, sess, triggertime{tt});
         triggertimenow = filtertriggertime(triggertimenow,evkeyword, evkeytype, evName, evType, evT); 
         if isempty(triggertimenow)
             disp(['-------> warning: trigger times not found in ', fdir, ' -> ', sess, ': ', triggername{tt}]);
         else
              if ~isempty(spikedbfile)
                 [avgnow, xbinnow] = findotherdbavgs(pinfo, data, dataindthere, fdir, sess, triggertimenow, binsize, maxlag, isavg);
              else
                 fname = eeg.general.eegfile{i};
                 if ~strcmp(eeg.parm.band{i}, 'cluster0')
                    [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
                 else
                    dat = eegdata.cluster0.cnt{i}; fs = 1/eeg.parm.cl0Binsize(i);  timestamp = eegdata.cluster0.timepoint{i};
                 end
                 [avgnow, xbinnow] = findtriggeredaverage(dat, timestamp, fs, binsize, maxlag, triggertimenow, isavg);
              end
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
       %%%get primary/secondary peak height values 
       [pval, sval, pheight, sheight, pratio, sratio] = findpeakheightvalues(avg, xbin, mm, nn); %, xbin, se, triggername{tt}, strcat(grpname, '_', evkeyword, '_', evkeytype));
       allv = []; allv.pval = pval; allv.sval = sval; allv.pheight = pheight; allv.sheight = sheight; allv.pratio = pratio; allv.sratio = sratio;
       grpnamenow = [];
       for (kt = 1:numel(grpname))
           grpnamenow = strcat(grpnamenow, '_', grpname{kt});
       end
       se = std(avg)/sqrt(mm); avg = mean(avg); ssnow = strcat('RipTrigAvg', grpnamenow, '_', evkeyword, '_', evkeytype);
       plottrigavg(avg, xbin, se, triggername{tt}, strcat(grpname, '_', evkeyword, '_', evkeytype));
       %display values as text messages
       valuelistpos = [0.001 0.001 0.98-0.001 0.95-0.001];
       valuepos = [0.001 valuelistpos(2) valuelistpos(3) valuelistpos(4)-0.001];
       valuetitle = {'Primary value', ; 'Secondary value'; 'Primary height'; 'Secondary height'; 'Primary ratio'; 'Secondary ratio'};
       valuecol{1} = cell(1,mm); valuecol{2} = cell(1,mm); valuecol{3} = cell(1,mm); valuecol{4} = cell(1,mm); valuecol{5} = cell(1,mm); valuecol{6} = cell(1,mm); 
       for (ij = 1:mm)
           valuecol{1}{ij} = num2str(pval(ij)); valuecol{2}{ij} = num2str(sval(ij));
           valuecol{3}{ij} = num2str(pheight(ij)); valuecol{4}{ij} = num2str(sheight(ij));
           valuecol{5}{ij} = num2str(pratio(ij)); valuecol{6}{ij} = num2str(sratio(ij));
       end
       hfnow = figure('Name', ssnow, 'NumberTitle', 'off', 'NextPlot', 'add',...
              'MenuBar', 'figure', 'Units', 'normalized', 'Position', [0.05 0.2 0.9 0.7]);
       TextDisplayer_multiple(hfnow, valuepos, valuecol, valuetitle, 'normalized'); setappdata(hfnow, 'allv', allv);
    else
       disp('-------> nothing to average');
    end
end
end

function [pval, sval, pheight, sheight, pratio, sratio] = findpeakheightvalues(avg, xbin, mm, nn)
pheight = zeros(1,mm); sheight = zeros(1, mm); pratio = zeros(1,mm); sratio = zeros(1, mm);
pval = zeros(1,mm); sval = zeros(1, mm); 
ppbin = find( (xbin>= -0.0001) & (xbin<=0.0011) ); %%primary peak bins - take averages of values at [0 1] ms
ptbin = find( (xbin>=0.0039) & (xbin<=0.0051) ); %%primary trough bins - take averages of values at [4 5] ms
spbin = find( (xbin>= -0.0061) & (xbin<= -0.0049) ); %%secondary peak bins - take averages of values at [-6 -5] ms
stbin = find( (xbin>= -0.0031) & (xbin<= -0.0019) ); %%secondary trough bins - take averages of values at [-3 -2] ms
%centerleftbin = round(nn/4); centerrightbin = round(3*nn/4);
for (i = 1:mm)
    curve = avg(i,:); 
    pval(i) = mean( curve(ppbin) ); ptroughval = mean( curve(ptbin) );
    sval(i) = mean( curve(spbin) ); stroughval = mean( curve(stbin) );
    pheight(i) = pval(i)-ptroughval; pratio(i) = pval(i)/ptroughval;
    sheight(i) = sval(i)-stroughval; sratio(i) = sval(i)/stroughval;
    %%%%%%%%%%%Below is to find the real peak of each curve, but most of them do not have real peaks
%     curvenow = curve; %smoothcrr(curve,3); 
%     %hf = figure; plot(1:numel(curve), curve);
%     [minindex, maxindex] = FindLocal(curvenow);
%     minindex = minindex( (minindex>=centerleftbin) & (minindex<=centerrightbin) );
%     maxindex = maxindex( (maxindex>=centerleftbin) & (maxindex<=centerrightbin) );
%     [minindex, maxindex] = indexcorrection(minindex, maxindex, curvenow, 2); %%%need to correct it by comparing with its neighbor to precisely determine the peak/trough indices
%     if ~isempty(maxindex)
%        [pval(i), iii] = max(curve(maxindex)); pind = maxindex(iii); ptroughind = findnexttrough(minindex, pind);
%        if ~isempty(ptroughind)
%           ptroughval = curve(ptroughind); pheight(i) = pval(i)-ptroughval; pratio(i) = pval(i)/ptroughval;
%        end
%        maxindex = setdiff(maxindex, pind);
%        if ~isempty(maxindex)
%           [sval(i), iii] = max(curve(maxindex)); sind = maxindex(iii); stroughind = findnexttrough(minindex, sind);
%           if ~isempty(stroughind)
%              stroughval = curve(stroughind); sheight(i) = sval(i)-stroughval; sratio(i) = sval(i)/stroughval;
%           end
%        end
%        %if ((i==1)||(i==mm)) disp([pind ptroughind sind stroughind]); end
%    end
end
function troughind = findnexttrough(minindex, pind)
troughind = []; iii = find(minindex-pind>0); 
if ~isempty(iii) 
    troughind = min(minindex(iii));
end
function newcrr = smoothcrr(crr, smoothbin)
newcrr = zeros(size(crr)); %NaN*ones(size(crr)); 
halfbin = floor(smoothbin/2);
for (i = 1:numel(crr))
    kkk = [i-halfbin:i+halfbin]; kkk = kkk( (kkk>=1) & (kkk<=numel(crr)) );
    crrnow = crr(kkk); newcrr(i) = mean(crrnow(~isnan(crrnow)));
end
%hf = figure; plot(1:numel(newcrr), newcrr);
function [minindex, maxindex] = indexcorrection(minindex, maxindex, curve, nnei)
nn = numel(curve);
for (i = 1:numel(minindex))
    if ((minindex(i)-nnei)>=1) && (minindex(i)+nnei<=nn) 
        %iii = [minindex(i)-2 minindex(i)-1 minindex(i) minindex(i)+1 minindex(i)+2]; 
        iii = minindex(i)-nnei:minindex(i)+nnei;
        curvenow = curve(iii); [~,ij] = min(curvenow); minindex(i) = iii(ij);
    end
end
for (i = 1:numel(maxindex))
    if ((maxindex(i)-1)>=1) && (maxindex(i)+1<=nn) 
        iii = [maxindex(i)-1 maxindex(i) maxindex(i)+1]; curvenow = curve(iii); 
        [~,ij] = max(curvenow); maxindex(i) = iii(ij);
    end
end

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

function triggertimenow = matchtimeout(eeg, eegdata, fdir, sess, triggertime)
triggertimenow = [];
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
    else
        disp(['-------> warning: multiple matches in ', fdir, ' -> ', sess]);
    end
end
function triggertime = filtertriggertime(triggertime,evkeyword, evkeytype, evName, evType, evT) 
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
startT = []; entT = [];
for (i = 1:numel(evT))
    startT = [startT evT{i}.start]; entT = [entT evT{i}.ent];
end
if ~isempty(startT)
    iii = [];
    for (i = 1:numel(startT))
        iiok = find( (triggertime>=startT(i)) & (triggertime<=entT(i)) );
        iii = union(iii, iiok);
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

function [minindex, maxindex] = FindLocal(dat)
%%NOW IT WORKS REALLY WELL!!!!
%%the problem with the derivative polarity change search is that EEG traces
%%are not that smooth, it gets more rugged on peaks and troughs
%%solved by smoothing EEG traces
deriv = diff(dat); %difference of EEG points
deriva = deriv(1:numel(deriv)-1);
derivb = deriv(2:numel(deriv)); %%deriva and derivb shift by one point
pola = deriva .* derivb; %chech for polarity change (product <0)
mindex = find(pola <= 0); %mindex are indices in pola, deriva and derivb, mindex+1 is the real singular point index in dat
%disp(strcat('-----> totally:', num2str(numel(mindex)), ' found'));
if (numel(mindex) < 2)
    minindex = []; maxindex = []; %if only o or 1 singular point: error
else
    %for each polarity change, check if minima or maxima
    mvalue = ClassPol(mindex, deriva);
    %now classify singular values to minima and maxima
    maxindex = find(mvalue == 1); %return real maximum indices in dat
    minindex = find(mvalue == -1); %return real minimum indices in dat
end
nm = numel(minindex); nM = numel(maxindex);
if (abs(nm-nM)>1) disp('-----------> Warning: minima and maxima do not match!'); end
nn = min([nm nM]); minindex = minindex(1:nn); maxindex = maxindex(1:nn);

function mvalue = ClassPol(mindex, deriva)
%have to classify in detail to get an accurate peaks and troughs
%mvalue = 1 if maxima, =-1 if minima, otherwise =0
%strategy: deal with point one by one, slope1~=0, then move to next until another
%slope2~=0, mvalue and position decide by slope1*slope2
mvalue = zeros(1, numel(deriva)); %initial assignment =0
pointnow = 1;
pointend = numel(mindex);
while (pointnow < pointend)
    slope1 = deriva(mindex(pointnow));
    zeropoint = 0; %how many zeros on peaks or troughs
    if (slope1 > 0) %if up maximum candidate
        while (mindex(pointnow)+1+zeropoint < numel(deriva))
            if (deriva(mindex(pointnow)+1+zeropoint) == 0) %if next point flat (zero) go on getting next point
                zeropoint = zeropoint + 1;
            else
                break
            end
        end
        if (mindex(pointnow)+1+zeropoint < numel(deriva)) %if still valid index
            if (deriva(mindex(pointnow)+1+zeropoint) < 0) %if next none-zero point is down then a maximum
               mvalue(mindex(pointnow)+floor((zeropoint+1)/2)) = 1; %take middle point in the zero points
            end
        end
    elseif (slope1 < 0) %if down minimum candidate
        while (mindex(pointnow)+1+zeropoint < numel(deriva))
            if (deriva(mindex(pointnow)+1+zeropoint) == 0) %if next point flat (zero) go on getting next point
                zeropoint = zeropoint + 1;
            else
                break
            end
        end
        if (mindex(pointnow)+1+zeropoint < numel(deriva)) %if still valid index
            if (deriva(mindex(pointnow)+1+zeropoint) > 0) %if next none-zero point is down then a maximum
                mvalue(mindex(pointnow)+floor((zeropoint+1)/2)) = -1; %take middle point in the zero points
            end
        end
    end
    pointnow = pointnow + zeropoint + 1;
end



