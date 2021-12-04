function [eeg, eegdata] = DataManager_FindRippleProp(eeg, eegdata, eegind, vv)

%%%%first add filtered ripple-band EEG traces
ttt = find( strcmp(eeg.parm.band(eegind), 'ripple') ); 
if numel(ttt) == 0
   SS = questdlg(['Found no ripple files. Add/create ripple files?']);
   if (strcmp(SS, 'Yes'))
    %%%%%read hvs filer
    [MCroot, MCname, DAname, DEname] = CurrentVersion;
    filterfile = fullfile(MCroot, 'Filters', 'EEGripple.mat'); %filter = EEG ripple, keep original sampling rate
    s = load(filterfile); num = s.Num; den = s.Den; s = []; %%denumerator 
    %%%%%filter the broad-band EEG traces
    iii = find( ~strcmp(eeg.parm.band(eegind), 'ripple') ); 
    for (i = 1:numel(iii))
        eidnow = eegind(iii(i));
        eegfilenow = eeg.general.eegfile{eidnow}; buffersize = eeg.parm.buffersize(eidnow); infix = 'ripple';
        [pp, nn, ee] = fileparts(eegfilenow); ID = strcat(nn, '_', infix); newfilename = fullfile(pp, strcat(ID, ee));
        if (exist(newfilename, 'file') ~= 2)
            disp(['--------> filter file: ', eegfilenow]);
            EEG_FilterData(num, den, eegfilenow, buffersize, infix, vv); %filter data and output to same directory with name+infix
        end
        [eeg, eegdata] = appendfile(eeg, eegdata, eegfilenow, infix, eidnow);
    end
    sss = find( strcmp(eeg.parm.band, 'ripple') );
    disp(['--------> Added ', num2str(numel(sss)-numel(ttt)), ' ripple files.']);
   elseif (strcmp(SS, 'Cancel'))
    disp(['--------> action cancelles']);
   elseif (strcmp(SS, 'No'))
     TT = questdlg(['Detect ripple-like events in non-ripple files?']);
     if (strcmp(TT, 'Yes'))  
        [eeg, eegdata] = DoRippleAnalysis(eeg, eegdata, eegind);
     end
   end
else
   [eeg, eegdata] = DoRippleAnalysis(eeg, eegdata, eegind);
end
disp('******************');

function [eeg, eegdata] = DoRippleAnalysis(eeg, eegdata, eegind)
%%%detect ripple events and compute their parameters
neeg = numel(eeg.general.eegfile);
%%%the following variables are assigned
% if (~isfield(eegdata, 'ripple')) eegdata.ripple = []; end
% if (~isfield(eegdata.ripple, 'allripples')) eegdata.ripple.allripples = cell(1, neeg); end %%all computed results are here
%%%%displayed variables
if (~isfield(eeg, 'ripple')) eeg.ripple = []; end
if (~isfield(eeg.ripple, 'sessstd')) eeg.ripple.sessstd = cell(1, neeg); end
if (~isfield(eeg.ripple, 'startThreshold')) eeg.ripple.startThreshold = cell(1, neeg); end
if (~isfield(eeg.ripple, 'peakThreshold')) eeg.ripple.peakThreshold = cell(1, neeg); end

if (~isfield(eeg.ripple, 'sessNum')) eeg.ripple.sessNum = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessRate')) eeg.ripple.sessRate = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessMeanAmp')) eeg.ripple.sessMeanAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessMeanDur')) eeg.ripple.sessMeanDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessMeanEng')) eeg.ripple.sessMeanEng = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessMeanFreq')) eeg.ripple.sessMeanFreq = cell(1, neeg); end

if (~isfield(eeg.ripple, 'evtNum')) eeg.ripple.evtNum = cell(1, neeg); end
if (~isfield(eeg.ripple, 'evtRate')) eeg.ripple.evtRate = cell(1, neeg); end
if (~isfield(eeg.ripple, 'evtMeanAmp')) eeg.ripple.evtMeanAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'evtMeanDur')) eeg.ripple.evtMeanDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'evtMeanEng')) eeg.ripple.evtMeanEng = cell(1, neeg); end
if (~isfield(eeg.ripple, 'evtMeanFreq')) eeg.ripple.evtMeanFreq = cell(1, neeg); end

if (~isfield(eeg.ripple, 'StartT')) eeg.ripple.StartT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'EndT')) eeg.ripple.EndT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'PeakT')) eeg.ripple.PeakT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'Dur')) eeg.ripple.Dur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'Amp')) eeg.ripple.Amp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'Eng')) eeg.ripple.Eng = cell(1, neeg); end
if (~isfield(eeg.ripple, 'Freq')) eeg.ripple.Freq = cell(1, neeg); end

for (iiik = 1:numel(eegind))
i = eegind(iiik);
%if strcmp(eeg.parm.band{i}, 'ripple') %%only do this for filtered EEG traces - remove this now
    disp(['--------> detecting ripple events: ', eeg.general.eegfile{i}]);
    ripnorm = eeg.parm.rippNorm{i}; pth = eeg.parm.rippPeakThres(i); sth = eeg.parm.rippStartThres(i); 
    maxgap = eeg.parm.rippMaxGap(i)*eeg.general.freq{i}; %%%split an event if gap in datapoints longer than this (e.g. 0.03s - missing 3 cycles at 100Hz)
    mingap = eeg.parm.rippMaxGap(i)*eeg.general.freq{i}; %%%%%%merge events if gap between ripple events in datapoints smaller than this
    mindur = eeg.parm.rippMinDur(i)*eeg.general.freq{i};
    maxdur = eeg.parm.rippMaxDur(i)*eeg.general.freq{i}; maxamp = eeg.parm.rippMaxAmp(i);
    peakmode = eeg.parm.rippPeakMode{i};
    %%%get eeg data
    [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
    dat = dat-mean(dat);
    %%%find parameters for ripple detection
    ssst = std(dat); maxamp = maxamp*ssst;
    if (strcmp(ripnorm, 'self'))
        pth = pth*ssst; sth = sth*ssst; disp(['---------------> self std: ', num2str(ssst)]);
    elseif (strcmp(ripnorm, 'run'))
        runstd = findrunsessstd(eeg, eegdata, eeg.general.finaldir{i}, 'ripple', eeg.general.recarea{i}); 
        if (~isempty(runstd))
           pth = pth*runstd; sth = sth*runstd; disp(['---------------> run std: ', num2str(runstd)]);
        else
           disp('---------------> Warning: run eeg ripple data not found; use self std instead as threshold unit');
           pth = pth*std(dat); sth = sth*std(dat); disp(['---------------> self std: ', num2str(std(dat))]);
        end
    end
    eeg.ripple.sessstd{i} = ssst;
    eeg.ripple.startThreshold{i} = sth; eeg.ripple.peakThreshold{i} = pth;
    %%%detect ripples
    ripp = findallripples(timestamp, dat, pth, sth, mingap, maxgap, mindur, maxdur, maxamp, peakmode); %%ripp.startT(i), .peakT, .endT, .dur, .amp
    
    if (~isempty(eeg.general.sesslength{i})) && (eeg.general.sesslength{i}>0)
       %%%session variables 
       eeg.ripple.sessNum{i} = numel(ripp.startT); eeg.ripple.sessRate{i} = numel(ripp.startT)/eeg.general.sesslength{i};
       eeg.ripple.sessMeanAmp{i} = mean(ripp.amp); eeg.ripple.sessMeanDur{i} = mean(ripp.dur); 
       eeg.ripple.sessMeanEng{i} = mean(ripp.eng); eeg.ripple.sessMeanFreq{i} = mean(ripp.freq);
       eeg.ripple.StartT{i} = ripp.startT; eeg.ripple.PeakT{i} = ripp.peakT; eeg.ripple.EndT{i} = ripp.endT; 
       eeg.ripple.Amp{i} = ripp.amp; eeg.ripple.Dur{i} = ripp.dur; eeg.ripple.Eng{i} = ripp.eng; eeg.ripple.Freq{i} = ripp.freq;
    
       %%%event variables
       evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i}; 
       for (j = 1:numel(evTime))
           eeg.ripple.evtNum{i}(j) = 0; eeg.ripple.evtRate{i}(j) = 0; eeg.ripple.evtMeanAmp{i}(j) = NaN; 
           eeg.ripple.evtMeanDur{i}(j) = NaN; eeg.ripple.evtMeanEng{i}(j) = NaN;
           st = evTime{j}.start; et = evTime{j}.ent; nev = numel(st); evtid = []; leng = 0;
           for (ttt = 1:nev)
               ik = find( (ripp.peakT>=st(ttt)) & (ripp.peakT<et(ttt)) ); evtid = union(evtid, ik); ik = []; leng = leng + et(ttt)-st(ttt);
           end
           if (leng > 0)
               eeg.ripple.evtNum{i}(j) = numel(ripp.startT(evtid)); eeg.ripple.evtRate{i}(j) = numel(ripp.startT(evtid))/leng;
               eeg.ripple.evtMeanAmp{i}(j) = mean(ripp.amp(evtid)); eeg.ripple.evtMeanDur{i}(j) = mean(ripp.dur(evtid));
               eeg.ripple.evtMeanEng{i}(j) = mean(ripp.eng(evtid)); eeg.ripple.evtMeanFreq{i}(j) = mean(ripp.freq(evtid));
           end
       end
    end
%end
end

function ripp = findallripples(timestamp, dat, pth, sth, mingap, groupspan, mindur, maxdur, maxamp, peakmode) %%ripp.startT(i), .peakT, .endT, .dur, .amp
ripp.startT = []; ripp.peakT = []; ripp.endT = []; ripp.dur = []; ripp.amp = []; ripp.eng = []; ripp.freq = [];
if strncmpi(peakmode, 'abs', 3) ||  strncmpi(peakmode, 'both', 3)
   index = find(abs(dat) >= sth); %all threshold crossing for both positive crossing and negative crossing
elseif strncmpi(peakmode, 'neg', 3)   
   index = find(dat <= -sth); %all threshold negative crossing ---Use this one because MUA agrees with max negative peaks more than the max positive peaks
elseif strncmpi(peakmode, 'pos', 3)   
   index = find(dat >= sth); %all threshold positive crossing
else
   index = [];
end
if (~isempty(index))
    [startindex, endindex] = findeventindex(index, groupspan); %%start/end indices for start threshold crossing
    %%%%now do a filtering to select those with durations within [mindur maxdur]
    alldur = endindex-startindex+1; iii = find( (alldur>=mindur) & (alldur<=maxdur) );
    startindex = startindex(iii); endindex = endindex(iii);
    if (~isempty(startindex))
       %%%%now join events if event gaps are too small 
       [startindex, ttt] = sort(startindex); endindex = endindex(ttt);
       gapnow = startindex(2:numel(startindex)) - endindex(1:numel(startindex)-1);
       smallgapind = find( gapnow < mingap ); largegapind = setdiff( [1:numel(gapnow)], smallgapind );
       newstart = startindex(largegapind+1); newend = endindex(largegapind);
       addstart = []; addend = [];
       if (~isempty(smallgapind))
          [sind, eind] = findeventindex(smallgapind, 1); %%start/end indices for start threshold crossing
          addstart = startindex(sind); addend = endindex(eind+1); 
       end
       startindex = union(union(startindex(1), newstart), addstart);
       endindex = union(union(endindex(numel(endindex)), newend), addend);
       startindex = sort(startindex); endindex = sort(endindex);
       %%%%%%%compute parameters in the events
       ripp = findamp(dat, timestamp, startindex, endindex, pth, maxamp, peakmode);
       %disp([startindex(1:100)' endindex(1:100)']);
    end
end

function [startindex, endindex] = findeventindex(index, groupspan)
%%index = a vector of integer indices of EEG file, in ascending order 
%%groupspan = any points with gap smaller than group span belong to same event
%%[startindex endindex] = detected event start index and end index, real numbers in index
startindex = []; endindex = [];
if (numel(index) == 1)
   startindex = index(1); endindex = index(1); 
else
   index = sort(index); %% with indexing [1,NN]
   indexdiff = diff(index); %get difference with indexing [1, NN-1] = original [2, NN]
   spanindex = find(indexdiff > groupspan); %%get all gaps bigger than groupspan, spanindex are indices in indexdiff
   if (isempty(spanindex)) %%if all within groupspan = only one event
       startindex = index(1); endindex = index(numel(index));
   else
       startindex(1) = index(1);
       endindex(1) = index(spanindex(1));
       nindex = numel(spanindex);
       for (i = 2:nindex)
           startindex(i) = index(spanindex(i-1)+1);
           endindex(i) = index(spanindex(i));
       end
       startindex(nindex+1) = index(spanindex(nindex)+1);
       endindex(nindex+1) = index(numel(index));
   end
end

function ripp = findamp(dat, timestamp, startindex, endindex, pth, maxamp, peakmode)
%%individual event indicated by [startindex, endindex], indices across startthreshold
%%peaktime, amplitude = maximum values within [startindex endindex] (minimum for slow wave detection)
%%starttime = index in [startindex-groupspan startindex] but cross the startthreshold
%%endtime = index in [endindex endindex+groupspan]but cross the startthreshold
%%area = integrate through [tarttime endtime]
ripp.startT = []; ripp.peakT = []; ripp.endT = []; ripp.dur = []; ripp.amp = []; ripp.eng = []; ripp.freq = [];
nevent = 0;
for (i = 1:numel(startindex))
     datnow = dat(startindex(i):endindex(i)); timenow = timestamp(startindex(i):endindex(i));
     if strncmpi(peakmode, 'pos', 3)
        [peaknow, peakindex] = max(datnow); peak1 = peaknow; peak2 = peaknow; %%if either positive crossing
     elseif strncmpi(peakmode, 'neg', 3)
        [peaknow, peakindex] = min(datnow); peak1 = -peaknow; peak2 = -peaknow; %%if either negative crossing
     elseif strncmpi(peakmode, 'abs', 3) % if either
        [peaknow, peakindex] = max(abs(datnow)); peak1 = peaknow; peak2 = peaknow;
     elseif strncmpi(peakmode, 'both', 3)
        [peaknow, peakindex] = max(abs(datnow)); peak1 = -1*min(datnow); peak2 = max(datnow); %%if both positive and negative crossing
     else
         peaknow = NaN; peak1 = NaN; peak2 = NaN; peakindex = NaN;
     end
     if (peak1 >= pth) && (peak2 >= pth) && (abs(peaknow) < maxamp)
         nevent = nevent + 1;
         ripp.amp(nevent) = datnow(peakindex); 
         ripp.peakT(nevent) = timenow(peakindex);
         ripp.startT(nevent) = timestamp(startindex(i)); ripp.endT(nevent) = timestamp(endindex(i));
         ripp.dur(nevent) = ripp.endT(nevent) - ripp.startT(nevent);
         ripp.eng(nevent) = mean(abs(datnow))* ripp.dur(nevent);
         nm = findlocal(datnow);
         ripp.freq(nevent) = nm/2/ripp.dur(nevent);
     end
end

function runstd = findrunsessstd(eeg, eegdata, finaldir, bandinfix, recarea)
runstd = []; dat = [];
%%%look for all the run events and select the data durign these run events
iii = find( strcmp(eeg.general.finaldir, finaldir) & strcmp(eeg.parm.band, bandinfix) & strcmp(eeg.general.recarea, recarea) ); % &(strcmp(eeg.parm.sessType, 'linear')) );
%%%%%%%%%%first look for run session data filtered thru the non-zero speed
for (kkk = 1:numel(iii))
    i = iii(kkk); datok = [];
    sessType = eeg.parm.sessType{i};
    if strncmpi(sessType, 'linear', 4) || strncmpi(sessType, 'open', 4) %%%%if linear or open 
        [timestamp, datnow, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
        datnow = datnow-mean(datnow); 
        [runstart, runend] = findrunepisodes(eeg.general.finaldir{i}, eeg.general.sessname{i});
        if (~isempty(runstart))
           for (k = 1:numel(runstart))
               datok = [datok; datnow( (timestamp>=runstart(k)) & (timestamp<=runend(k)) )];
           end
        end
        dat = [dat; datok]; datnow = [];
    end
end
if (~isempty(dat))
    runstd = std(dat);
else %%%%if non-zero speed run session data can not be defined, look for just run session data without the filtering
    for (kkk = 1:numel(iii))
        i = iii(kkk); datok = [];
        sessType = eeg.parm.sessType{i};
        if strncmpi(sessType, 'linear', 4) || strncmpi(sessType, 'open', 4) %%%%if linear or open 
            [timestamp, datnow, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
            datok = datnow-mean(datnow);
        end
        dat = [dat; datok]; datnow = [];
    end
    if (~isempty(dat))
        runstd = std(dat);
        disp('---------------> Warning: motion power file not found or no running detected; use whole run session for computing thresholds');
    end
end

function [runstart, runend] = findrunepisodes(finaldir, sessname)
runstart = []; runend = []; powerfile = [];
posdir{1} = strcat(finaldir, filesep, 'pos');
[filepath, filename] = GetAllFile(posdir, 'motion', '.pwr'); 
for  (j = 1:numel(filename)) %%%green diodes first: better signal in most of recordings
     if (~isempty(strfind(filename{j}, 'green'))) && (~isempty(strfind(filename{j}, sessname)))
        powerfile = fullfile(filepath{j}, filename{j}); break
     end
end
if isempty(powerfile)
   for (j = 1:numel(filename))
        if (~isempty(strfind(filename{j}, sessname)))
            powerfile = fullfile(filepath{j}, filename{j}); break
        end
   end
end
if (~isempty(powerfile))
    [data, ~, ~, ~, ~] = ReadPowergramData(powerfile);
    timenow = 0.0001*data(:,1); pwrnow = data(:,2); %%%% timestamps in power data are in 0.1ms
    state = ones(size(pwrnow)); %%%initial assignment: intermediate state
    iii = find(pwrnow >=20); state(iii) = 2*ones(size(iii)); %%%%running state: >50 fixed threshold in pixels per second
    %jjj = find(pwrnow <=30); state(jjj) = zeros(size(jjj)); %%%%stopp state: <30 fixed threshold in pixels per second
    [stateStart, stateEnd] = findstateepisode(state, 1, numel(state), 2);
    runstart = timenow(stateStart); runend = timenow(stateEnd);
    iii = find(runend-runstart>1); runstart = runstart(iii); runend = runend(iii); %%%select events longer than 1s
end
%disp([runstart runend]);
    
function [stateStart, stateEnd] = findstateepisode(state, starttime, endtime, value)
%find continous sleep episode with the same sleep state value (0=wake, 1=tran, 2=sws, 3=rem)
%state = classified point by point sleep state, 
%starttime/endtime = start/end time point of sleep classification
stateStart = []; stateEnd = [];
istate= 0; %number of episode
i = starttime;
while (i < endtime)
    if (i == 1) 
        if (state(i) == value)&&(state(i+1)==value)
           istate = istate + 1;
           stateStart(istate) = i; %%assign the starting point
           j = i + 1;
           endflag = 0;  %%endflag = 0 continue to search for an ending point
           while ((j <= endtime) && (endflag == 0))
             if (j == endtime)   %%if this is the last point of the whole data
                stateEnd(istate) = j;
                i = j + 1;
                endflag = 1;
             else
                if ((state(j-1) == value) && (state(j) == value) && (state(j+1) ~= value))   %%detect state ending point
                    stateEnd(istate) = j; %%assign the ending point
                    endflag = 1;
                    i = j + 1;
                else
                    j = j +1;
                end
             end
           end
        else
           i = i + 1;
        end
    else
        if ((state(i) == value) && (state(i-1) ~=value) && (state(i+1) == value))   %%detect state starting point
           istate = istate + 1;
           stateStart(istate) = i; %%assign the starting point
           j = i + 1;
           endflag = 0;  %%endflag = 0 continue to search for an ending point
           while ((j <= endtime) && (endflag == 0))
             if (j == endtime)   %%if this is the last point of the whole data
                stateEnd(istate) = j;
                i = j + 1;
                endflag = 1;
             else
                if ((state(j-1) == value) && (state(j) == value) && (state(j+1) ~= value))   %%detect state ending point
                    stateEnd(istate) = j; %%assign the ending point
                    endflag = 1;
                    i = j + 1;
                else
                    j = j +1;
                end
             end
           end
        else
           i = i + 1;
        end
    end
end

function nm = findlocal(dat)
%%%find number of all local minima and maxima in dat
deriv = diff(dat); %difference of EEG points
deriva = deriv(1:numel(deriv)-1);
derivb = deriv(2:numel(deriv)); %%deriva and derivb shift by one point
pola = deriva .* derivb; %chech for polarity change (product <0)
mindex = find(pola <= 0); %mindex are indices in pola, deriva and derivb, mindex+1 is the real singular point index in dat
nm = numel(mindex);

function [eeg, eegdata] = appendfile(eeg, eegdata, eegfilenow, infix, index)
nfile = numel(eeg.general.eegfile); 
[pp, nn, ee] = fileparts(eegfilenow); ID = strcat(nn, '_', infix); newname = fullfile(pp, strcat(ID, ee));
eeg.general.eegfile{nfile+1} = newname; eeg.general.eegID{nfile+1} = ID;
eeg.parm.band{nfile+1} = infix;
genvar = fieldnames(eeg.general);
for (i = 1:numel(genvar))
    if (~strcmp(genvar{i}, 'eegfile')) && (~strcmp(genvar{i}, 'eegID'))
        eeg.general.(genvar{i}){nfile+1} = eeg.general.(genvar{i}){index};
    end
end
parmvar = fieldnames(eeg.parm);
for (i = 1:numel(parmvar))
    if (~strcmp(parmvar{i}, 'band'))
        if (iscell(eeg.parm.(parmvar{i})))
           eeg.parm.(parmvar{i}){nfile+1} = eeg.parm.(parmvar{i}){index};
        elseif (isnumeric(eeg.parm.(parmvar{i})))
           eeg.parm.(parmvar{i})(nfile+1) = eeg.parm.(parmvar{i})(index);
        end
    end
end
eegdata.event.eventtimes{nfile+1} = eegdata.event.eventtimes{index}; 
grpname = eegdata.grouplist.groupname; zeroind = strcmp(grpname, 'List0'); 
eegdata.grouplist.groupindex{zeroind} = 1:numel(eeg.general.eegfile);
%%%%%Also need to add empty variables for other types of analysis: ripple/theta etc.
catnames = fieldnames(eeg);
for (i = 1:numel(catnames))
    if (~strcmp(catnames{i}, 'general')) && (~strcmp(catnames{i}, 'parm'))
        if (~isempty(eeg.(catnames{i})))
            subvar = fieldnames(eeg.(catnames{i}));
            for (j = 1:numel(subvar))
                if (~isempty(eeg.(catnames{i}).(subvar{j})))
                    %disp([catnames{i}, '.', subvar{j}]);
                    eeg.(catnames{i}).(subvar{j}){nfile+1} = []; 
                end
            end
        end
    end
end


