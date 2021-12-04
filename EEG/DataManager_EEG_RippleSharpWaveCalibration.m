function DataManager_EEG_RippleSharpWaveCalibration
%%Plot the computed eegdata
hf = gcbf; eeg = getappdata(hf, 'eeg'); eegdata = getappdata(hf, 'eegdata');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
grpind = find(groupselection == 1); grpname = eegdata.grouplist.groupname(grpind);
cellind = [];
for (k = 1:numel(grpind))
    cellind = union(cellind, eegdata.grouplist.groupindex{grpind(k)});
end
%%%sort ripple band/broad band traces
ok = 1;
[riptraceind, broadtraceind, ok] = findripplebroadpairs(eeg, cellind);
%%%set up parameters
if ok
   parm = [];
   p = {'Bin size (s)'; 'Max lag (s)'};
   d = {'0.001'; '0.05'};
   II = inputdlg(p, 'Ripple average parameters', 3, d, 'on'); %%%resizable window
   if ~isempty(II)
        ripbinsize = str2num(II{1}); ripmaxlag = str2num(II{2});
   else
        ok = 0;
   end
   if ok
      p = {'Bin size (s)'; 'Max lag (s)'};
      d = {'0.01'; '0.1'};
      II = inputdlg(p, 'Sharp wave (broad) average parameters', 3, d, 'on'); %%%resizable window
      if ~isempty(II)
         brdbinsize = str2num(II{1}); brdmaxlag = str2num(II{2});
      else
         ok = 0;
      end 
   end
end
if ok
    input = inputdlg({'Event keyword';'Event type'; 'Average over session?'}, 'Event slelection parms', 3, {'sws'; 'sws'; 'yes'}); 
    if (~isempty(input))
        parm.evkeyword = input{1}; parm.evkeytype = input{2}; 
        if strncmpi(input{3}, 'yes', 1)
            parm.isavg = 1;
        else
            parm.isavg = 0;
        end
    else
        ok = 0;
    end
end
%%%check if ripple.peakT exist
if ok
    if ~isfield(eeg, 'ripple')
        ok = 0; disp('----> ripples not defined');
    else
        if isfield(eeg.ripple, 'PeakT')
            varnow = eeg.ripple.PeakT; varstr = 'PeakT';
        elseif isfield(eeg.ripple, 'sessPeakT')
            varnow = eeg.ripple.sessPeakT; varstr = 'sessPeakT';
        else
            ok = 0; disp('----> ripple peak times not defined');
        end
    end
end
if ok 
   rippleamps = cell(1, numel(riptraceind)); sharpwaveamps = cell(1, numel(broadtraceind));
   ripavg = cell(1, numel(riptraceind)); brdavg = cell(1, numel(broadtraceind)); ripSE = cell(1, numel(riptraceind)); brdSE = cell(1, numel(broadtraceind));
   for (i = 1:numel(riptraceind))
       disp(['----------> ripple file now (', num2str(i), 'out of ', num2str(numel(riptraceind)), '):', eeg.general.eegfile{riptraceind(i)}]); 
       [rippleamps{i}, ripavg{i}, ripSE{i}, ripxbin] = gettriggeredaverages(varnow, eeg, eegdata, riptraceind(i), parm, ripbinsize, ripmaxlag);
       disp(['--------------> paired broad eeg file now: ', eeg.general.eegfile{broadtraceind(i)}]); 
       [sharpwaveamps{i}, brdavg{i}, brdSE{i}, brdxbin] = gettriggeredaverages(varnow, eeg, eegdata, broadtraceind(i), parm, brdbinsize, brdmaxlag);
   end
   rippleamps = cell2mat(rippleamps); sharpwaveamps = cell2mat(sharpwaveamps);
   if numel(rippleamps) ~= numel(sharpwaveamps)
       disp('----> numbers of ripple and sharpwave data points do not match');
   else
       plotresults(rippleamps, sharpwaveamps, varstr, grpname);
   end
   for (i = 1:numel(riptraceind))
       plotavgtraces(ripavg{i}, ripSE{i}, ripxbin, 'ripple', eeg.general.eegfile{riptraceind(i)});
       plotavgtraces(brdavg{i}, brdSE{i}, brdxbin, 'broad', eeg.general.eegfile{broadtraceind(i)});
   end
   ffname = strcat('RipTriggeredSharpwaves', grpname{1}, '.mat');
   save(ffname, 'rippleamps', 'sharpwaveamps', 'ripavg', 'ripSE', 'ripxbin', 'brdavg', 'brdSE', 'brdxbin', '-mat');
end
disp('************************');

function [ripind, broadind, ok] = findripplebroadpairs(eeg, cellind)
ripind = []; broadind = []; ok = 1; nmatch = 0;
allvind = find(strcmp(eeg.parm.band(cellind), 'ripple'));
allhind = find(strcmp(eeg.parm.band(cellind), 'broad'));
for (i = 1:numel(allvind))
     vsess = eeg.general.sessname{cellind(allvind(i))};
     vfdir = eeg.general.finaldir{cellind(allvind(i))};
     for (j = 1:numel(allhind))
            hsess = eeg.general.sessname{cellind(allhind(j))}; 
            hfdir = eeg.general.finaldir{cellind(allhind(j))};
            if strcmp(vsess, hsess) && strcmp(vfdir, hfdir) 
                nmatch = nmatch + 1; ripind(nmatch) = cellind(allvind(i)); broadind(nmatch) = cellind(allhind(j));
                break
            end
     end
end
if nmatch == 0
    ok = 0; disp(['----------> no pairs of ripple and broad band traces found']);
end

function [amps, avg, SE, xbin] = gettriggeredaverages(varnow, eeg, eegdata, dateind, parm, binsize, maxlag)
amps = []; avg = []; SE = []; xbin = [];
i = dateind;
fdir = eeg.general.finaldir{i}; sess = eeg.general.sessname{i}; 
evName = eeg.general.eventname{i}; evType = eeg.parm.eventtype{i}; evT = eegdata.event.eventtimes{i}; 
triggertimenow = matchtimeout(eeg, eegdata, fdir, sess, varnow);
triggertimenow = filtertriggertime(triggertimenow, parm.evkeyword, parm.evkeytype, evName, evType, evT); 
if isempty(triggertimenow)
    disp(['-------> warning: trigger times not found in ', fdir, ' -> ', sess]);
else
    fname = eeg.general.eegfile{i};
    if ~strcmp(eeg.parm.band{i}, 'cluster0')
        if exist(fname, 'file')~=2
           fname = strrep(fname, 'I:', 'J:'); 
        end
        [timestamp, dat, gain, fs] = ReadEEGFile(fname, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
    else
        dat = eegdata.cluster0.cnt{i}; fs = 1/eeg.parm.cl0Binsize(i);  timestamp = eegdata.cluster0.timepoint{i};
    end
    [avg, SE, xbin] = findtriggeredaverage(dat, timestamp, fs, binsize, maxlag, triggertimenow, parm.isavg);
    [mm,nn] = size(avg); amps = NaN*ones(1, mm);
    for (i = 1:mm)
        [~, iii] = max(abs(avg(i,:))); amps(i) = avg(i, iii);
    end
end

function plotresults(rippleamps, sharpwaveamps, varstr, grpname)
hg = figure('Name', 'Ripple vs sharpwave amplitudes'); 
hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
xlabel ('Sharp wave amplitude (mV)'); ylabel('Ripple amplitude (mV)'); 
line(sharpwaveamps, rippleamps, 'Parent', hax, 'Line', 'none', 'Marker', '.', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1], 'MarkerSize', 20);
%%%now need to do a regression to see whether there is a significant correlation between the two
minsharp = -0.8; maxsharp = 0.8; str1 =[]; str2 =[];
iii = find( (sharpwaveamps >= minsharp) & (sharpwaveamps <= maxsharp) ); 
if ~isempty(iii)
   yy = rippleamps(iii); xx = sharpwaveamps(iii);
   [str1, str2, YY, XX] = DoCrrRegress(xx', yy'); %%%inputs are both column vectors
   line(XX, YY, 'Parent', hax, 'Color', [1 0 0]);
end
text('Parent', hax, 'Interpreter', 'none', 'String', str1, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
text('Parent', hax, 'Interpreter', 'none', 'String', str2, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92]);
text('Parent', hax, 'Interpreter', 'none', 'String', varstr, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.88]);
for (i = 1:numel(grpname))
     text('Parent', hax, 'Interpreter', 'none', 'String', grpname{i}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.84-(i-1)*0.04]);
end
function [str1, str2, YY, xx] = DoCrrRegress(xx, yy) %%%xx, yy column vectors
ind = find( (~isnan(xx)) & (~isnan(yy)) ); xx= xx(ind); yy = yy(ind); %%%get rid of the non-numbers
[RR, PP] = corrcoef(xx, yy); r = RR(1,2); p = PP(1,2); %%correlation and p value
str1 = strcat('R=', num2str(RR(1,2)), '; p= ', num2str(PP(1,2)));
aaa = size(xx); if (aaa(1)==1) xx = xx'; end
aaa = size(yy); if (aaa(1)==1) yy = yy'; end
n = numel(xx); Xparm = [ones(n,1) xx];
[BBB, BBBint, R, Rint, stat] =  regress(yy, Xparm, 0.05);
YY = Xparm * BBB; %regression result: line in (X, Y) = (xx, L4)
str2 = strcat('regression: b=',num2str(BBB(1)), '(', num2str(BBBint(1,1)), ',', num2str(BBBint(1,2)), ')', ...
    '; k=',num2str(BBB(2)), '(', num2str(BBBint(2,1)), ',', num2str(BBBint(2,2)), ')' );

function plotavgtraces(avg, se, xbin, band, filename)
hg = figure('Name', 'EEG average'); 
hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
xlabel ('Trigger time (s)'); ylabel('EEG amplitude (mV)');
line(xbin, avg, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]);
line(xbin, avg+se, 'Parent', hax, 'LineWidth', 0.5, 'Color', [1 0 0]);
line(xbin, avg-se, 'Parent', hax, 'LineWidth', 0.5, 'Color', [1 0 0]);
text('Parent', hax, 'Interpreter', 'none', 'String', filename, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
text('Parent', hax, 'Interpreter', 'none', 'String', band, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92]);

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
if (~isempty(evkeyword)) || (~isempty(evkeytype))
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
   for (i = 1:numel(evT))
    for (j = 1:numel(evT{i}.start))
        iiok = find( (triggertime>=evT{i}.start(j)) & (triggertime<=evT{i}.ent(j)) );
        iii = union(iii, iiok);
    end
   end
   triggertime = triggertime(iii);
end

function [avg, SE, xbinnow] = findtriggeredaverage(dat, timestamp, fs, binsize, maxlag, triggertimenow, isavg)
binnow = round(binsize*fs); %%%%need to re-sample the data, otherwise the computation is too slow
ind = 1:binnow:numel(dat); dat = dat(ind); timestamp = timestamp(ind);
maxlagbin = round(maxlag/binsize); xbinnow = (-maxlagbin:maxlagbin)*binsize; xbin = -maxlagbin:maxlagbin;
iii = find( (triggertimenow>timestamp(1)+maxlag) & (triggertimenow<timestamp(numel(timestamp))-maxlag) );
triggertimenow = triggertimenow(iii);
ntime = numel(triggertimenow); nbin = numel(xbinnow); avgnow = zeros(ntime, nbin);
avg = NaN*ones(1, nbin); SE = NaN*ones(1, nbin);
for (i = 1:ntime)
    [~,iinow] = min(abs(timestamp-triggertimenow(i)));
    iii = xbin + iinow;
    avgnow(i,:) = dat(iii);
end
if isavg
   if ntime > 1
      avg = mean(avgnow); SE = std(avgnow)/sqrt(ntime);
   elseif ntime == 1
      avg = avgnow;
   end
end





