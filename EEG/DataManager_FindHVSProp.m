function [eeg, eegdata] = DataManager_FindHVSProp(eeg, eegdata, eegind, vv)

%%%%first add filtered ripple-band EEG traces
ttt = find( strcmp(eeg.parm.band(eegind), 'hvs') ); nohvsflag = false;
if numel(ttt) == 0
   SS = questdlg('Found no hvs files. Add/create HVS files?');
   if (strcmp(SS, 'Yes'))
    %%%%%read hvs filer
    [MCroot, MCname, DAname, DEname] = CurrentVersion;
    filterfile = fullfile(MCroot, 'Filters', 'EEGhvsresample_spikes.mat'); %filter = EEG hvs, sampling freq needs to be ~200Hz
    s = load(filterfile); num = s.Num; den = s.Den; s = []; %%denumerator 
    %%%%%filter the broad-band EEG traces
    iii = find( ~strcmp(eeg.parm.band(eegind), 'hvs') ); 
    for (i = 1:numel(iii))
        eidnow = eegind(iii(i));
        eegfilenow = eeg.general.eegfile{eidnow}; buffersize = eeg.parm.buffersize(eidnow); infix = 'hvs';
        [pp, nn, ee] = fileparts(eegfilenow); ID = strcat(nn, '_', infix); newfilename = fullfile(pp, strcat(ID, ee));
        if (exist(newfilename, 'file') ~= 2)
            disp(['--------> filter file: ', eegfilenow]);
            EEG_ResampleFilterData(1, 10, num, den, eegfilenow, buffersize, infix, vv);
        end
        [eeg, eegdata] = appendfile(eeg, eegdata, eegfilenow, infix, eidnow);
    end
    sss = find( strcmp(eeg.parm.band, 'hvs') );
    disp(['--------> Added ', num2str(numel(sss)-numel(ttt)), ' hvs files.']);
   elseif (strcmp(SS, 'Cancel'))
    disp(['--------> action cancelles']);
   elseif (strcmp(SS, 'No'))
     TT = questdlg('Detect hvs-like events in non-hvs files?');
     if (strcmp(TT, 'Yes'))  
        nohvsflag = true;
        [eeg, eegdata] = DohvsAnalysis(eeg, eegdata, eegind, nohvsflag);
     end
   end
else
   [eeg, eegdata] = DohvsAnalysis(eeg, eegdata, eegind, nohvsflag);
end
disp('******************');

%%%%%% The following is the old working version
% %%%%first add filtered hvs-band EEG traces
% ttt = find( strcmp(eeg.parm.band, 'hvs') );
% SS = questdlg(['Found ', num2str(numel(ttt)), ' hvs files. Add/create more hvs files?']);
% ok = 1;
% if (strcmp(SS, 'Yes'))
%     %%%%%read hvs filer
%     [MCroot, MCname, DAname, DEname] = CurrentVersion;
%     filterfile = fullfile(MCroot, 'Filters', 'EEGhvsresample_spikes.mat'); %filter = EEG hvs, sampling freq needs to be ~200Hz
%     s = load(filterfile); num = s.Num; den = s.Den; s = []; %%denumerator 
%     %%%%%filter the broad-band EEG traces
%     iii = find( strcmp(eeg.parm.band, 'broad') & ~strcmp(eeg.general.recarea, 'MSC') ); 
%     neeg = numel(eeg.general.eegfile);
%     for (i = 1:numel(iii))
%         eegfilenow = eeg.general.eegfile{iii(i)}; buffersize = eeg.parm.buffersize(iii(i)); infix = 'hvs';
%         [pp, nn, ee] = fileparts(eegfilenow); ID = strcat(nn, '_', infix); newfilename = fullfile(pp, strcat(ID, ee));
%         if (exist(newfilename, 'file') ~= 2)
%             disp(['--------> filter file: ', eegfilenow]);
%             EEG_ResampleFilterData(1, 10, num, den, eegfilenow, buffersize, infix, vv);
%         end
%         [eeg, eegdata] = appendfile(eeg, eegdata, eegfilenow, infix, iii(i));
%     end
%     eegind = 1:numel(eeg.general.eegfile);
%     sss = find( strcmp(eeg.parm.band, 'hvs') );
%     disp(['Added ', num2str(numel(sss)-numel(ttt)), ' hvs files.']);
% elseif (strcmp(SS, 'Cancel'))
%     ok = 0; disp(['--------------> action cancelles']);
% end
% if ok
%    [eeg, eegdata] = DohvsAnalysis(eeg, eegdata, eegind);
% end

function [eeg, eegdata] = DohvsAnalysis(eeg, eegdata, eegind, nohvsflag)
%%%detect hvs events and compute their parameters
neeg = numel(eeg.general.eegfile);
%%%the following variables are assigned
% if (~isfield(eegdata, 'hvs')) eegdata.hvs = []; end
% if (~isfield(eegdata.hvs, 'allhvss')) eegdata.hvs.allhvss = cell(1, neeg); end %%all computed results are here
%%%%displayed variables
if (~isfield(eeg, 'hvs')) eeg.hvs = []; end
if (~isfield(eeg.hvs, 'startThreshold')) eeg.hvs.startThreshold = cell(1, neeg); end
if (~isfield(eeg.hvs, 'peakThreshold')) eeg.hvs.peakThreshold = cell(1, neeg); end

if (~isfield(eeg.hvs, 'sessNum')) eeg.hvs.sessNum = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessRate')) eeg.hvs.sessRate = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessMeanAmp')) eeg.hvs.sessMeanAmp = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessMeanFreq')) eeg.hvs.sessMeanFreq = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessMeanDur')) eeg.hvs.sessMeanDur = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessMeanEng')) eeg.hvs.sessMeanEng = cell(1, neeg); end

if (~isfield(eeg.hvs, 'evtNum')) eeg.hvs.evtNum = cell(1, neeg); end
if (~isfield(eeg.hvs, 'evtRate')) eeg.hvs.evtRate = cell(1, neeg); end
if (~isfield(eeg.hvs, 'evtMeanAmp')) eeg.hvs.evtMeanAmp = cell(1, neeg); end
if (~isfield(eeg.hvs, 'evtMeanFreq')) eeg.hvs.evtMeanFreq = cell(1, neeg); end
if (~isfield(eeg.hvs, 'evtMeanDur')) eeg.hvs.evtMeanDur = cell(1, neeg); end
if (~isfield(eeg.hvs, 'evtMeanEng')) eeg.hvs.evtMeanEng = cell(1, neeg); end

if (~isfield(eeg.hvs, 'sessStartT')) eeg.hvs.sessStartT = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessEndT')) eeg.hvs.sessEndT = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessPeakT')) eeg.hvs.sessPeakT = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessDur')) eeg.hvs.sessDur = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessAmp')) eeg.hvs.sessAmp = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessEng')) eeg.hvs.sessEng = cell(1, neeg); end
if (~isfield(eeg.hvs, 'sessFreq')) eeg.hvs.sessFreq = cell(1, neeg); end

if (~isfield(eeg.hvs, 'minTime')) eeg.hvs.minTime = cell(1, neeg); end
if (~isfield(eeg.hvs, 'maxTime')) eeg.hvs.maxTime = cell(1, neeg); end
if (~isfield(eeg.hvs, 'minAmp')) eeg.hvs.minAmp = cell(1, neeg); end
if (~isfield(eeg.hvs, 'maxAmp')) eeg.hvs.maxAmp = cell(1, neeg); end

for (iiik = 1:numel(eegind))
i = eegind(iiik);
if strcmp(eeg.parm.band{i}, 'hvs') %%only do this for filtered EEG traces remove this now
    disp(['--------> detecting HVS events: ', eeg.general.eegfile{i}]);
    ripnorm = eeg.parm.hvsNorm{i}; pth = eeg.parm.hvsPeakThres(i); sth = eeg.parm.hvsStartThres(i); 
    mingap = eeg.parm.hvsMaxGap(i); %*eeg.general.freq{i}; %%%now mingap between HVS events in datapoints
    maxgap = mingap;
    %if nohvsflag maxgap = resetmaxgap(eeg.parm.band{i}, eeg.general.eegfile{i}); end
    disp(['---------------> groupspan set as <: ', num2str(maxgap), 's']);
    maxgap = maxgap*eeg.general.freq{i}; %%%now maxgap in datapoints 
    mingap = mingap*eeg.general.freq{i}; %%%now mingap in datapoints 
    minPeakNum = eeg.parm.hvsMinPeakNum(i); %%%*eeg.general.freq{i}; ---!!!This changed from ripple detection = min number of peaks passing the peak threshold
    maxdur = eeg.parm.hvsMaxDur(i)*eeg.general.freq{i}; peakmode = eeg.parm.hvsPeakMode{i}; 
    %%%get eeg data
    [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
    dat = dat-mean(dat);
    %%%find parameters for hvs detection
    if (strcmp(ripnorm, 'self'))
        pth = pth*std(dat); sth = sth*std(dat); disp(['---------------> self std: ', num2str(std(dat))]);
    elseif (strcmp(ripnorm, 'run'))
        runstd = findrunsessstd(eeg, eegdata, eeg.general.finaldir{i}, eeg.parm.band{i}, eeg.general.recarea{i}); 
        if (~isempty(runstd))
           pth = pth*runstd; sth = sth*runstd; disp(['---------------> run std: ', num2str(runstd)]);
        else
           disp('---------------> Warning: run eeg HVS data not found; use self std instead as threshold unit');
           pth = pth*std(dat); sth = sth*std(dat); disp(['---------------> self std: ', num2str(std(dat))]);
        end
    end
    eeg.hvs.startThreshold{i} = sth; eeg.hvs.peakThreshold{i} = pth;
    %%%session variables
    if (~isempty(eeg.general.sesslength{i})) && (eeg.general.sesslength{i}>0)
        %%%detect hvss
        hvs = findallhvss(timestamp, dat, pth, sth, mingap, maxgap, minPeakNum, maxdur, peakmode); %%hvs.startT(i), .peakT, .endT, .dur, .amp
        %%%%%%%%%%% here need to filter thru the peak/trough times: choose the ones truly spike (negative/positive) peaks
        %%%%%%%%%%%       then match the positive/negative peaks to the negative/positive peaks
        hvs = adjustminmaxtimes(hvs, dat, timestamp, sth);
        %disp('-----ignore adjusting min max times');
        %%%%%%%%%%%%now assign session variables        
        eeg.hvs.sessNum{i} = numel(hvs.startT); eeg.hvs.sessRate{i} = numel(hvs.startT)/eeg.general.sesslength{i};
        eeg.hvs.sessMeanAmp{i} = mean(hvs.amp); eeg.hvs.sessMeanDur{i} = mean(hvs.dur); eeg.hvs.sessMeanEng{i} = mean(hvs.eng);
        eeg.hvs.sessStartT{i} = hvs.startT; eeg.hvs.sessPeakT{i} = hvs.peakT; eeg.hvs.sessEndT{i} = hvs.endT;
        eeg.hvs.sessAmp{i} = hvs.amp; eeg.hvs.sessDur{i} = hvs.dur; eeg.hvs.sessEng{i} = hvs.eng;
        eeg.hvs.sessFreq{i} = hvs.freq;
        eeg.hvs.minTime{i} = hvs.minTime; eeg.hvs.minAmp{i} = hvs.minAmp;
        eeg.hvs.maxTime{i} = hvs.maxTime; eeg.hvs.maxAmp{i} = hvs.maxAmp;
        eeg.hvs.sessMeanFreq{i} = mean(hvs.freq);
        
        %%%event variables
        evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i};
        for (j = 1:numel(evTime))
            eeg.hvs.evtNum{i}(j) = 0; eeg.hvs.evtRate{i}(j) = 0; eeg.hvs.evtMeanAmp{i}(j) = NaN;
            eeg.hvs.evtMeanFreq{i}(j) = NaN; eeg.hvs.evtMeanDur{i}(j) = NaN; eeg.hvs.evtMeanEng{i}(j) = NaN;
            st = evTime{j}.start; et = evTime{j}.ent; nev = numel(st); evtid = []; leng = 0;
            for (ttt = 1:nev)
                ik = find( (hvs.peakT>=st(ttt)) & (hvs.peakT<et(ttt)) ); evtid = union(evtid, ik); ik = []; leng = leng + et(ttt)-st(ttt);
            end
            if (leng > 0)
                eeg.hvs.evtNum{i}(j) = numel(hvs.startT(evtid)); eeg.hvs.evtRate{i}(j) = numel(hvs.startT(evtid))/leng;
                eeg.hvs.evtMeanAmp{i}(j) = mean(hvs.amp(evtid)); eeg.hvs.evtMeanDur{i}(j) = mean(hvs.dur(evtid));
                eeg.hvs.evtMeanEng{i}(j) = mean(hvs.eng(evtid)); eeg.hvs.evtMeanFreq{i}(j) = mean(hvs.freq(evtid));
            end
        end
    end
else %%%if not a hvs file, get rid of those parm assignments only for hvs files
     eeg.parm.hvsNorm{i} = []; eeg.parm.hvsPeakThres(i) = NaN; eeg.parm.hvsStartThres(i) = NaN; 
     eeg.parm.hvsMaxGap(i) = NaN; eeg.parm.hvsMinDur(i) = NaN; eeg.parm.hvsMaxDur(i) = NaN;
end
end

function hvs = findallhvss(timestamp, dat, pth, sth, mingap, groupspan, minPeakNum, maxdur, peakmode) %%hvs.startT(i), .peakT, .endT, .dur, .amp
hvs.startT = []; hvs.peakT = []; hvs.endT = []; hvs.dur = []; hvs.amp = []; hvs.eng = []; 
hvs.minTime = []; hvs.maxTime = []; hvs.minAmp = []; hvs.maxAmp = []; hvs.freq = [];
if strncmpi(peakmode, 'abs', 3)
   index = find(abs(dat) >= sth); %all threshold crossing for both positive crossing and negative crossing
elseif strncmpi(peakmode, 'neg', 3)   
   index = find(dat <= -sth); %all threshold negative crossing ---Use this one because MUA agrees with max negative peaks more than the max positive peaks
elseif strncmpi(peakmode, 'pos', 3)   
   index = find(dat >= sth); %all threshold positive crossing
else
   index = [];
end
if (~isempty(index))
    %[startindex, endindex] = findeventindex(index, maxgap); %%start/end indices for start threshold crossing
    [startindex, endindex] = findeventindex(index, groupspan); %%start/end indices for start threshold crossing
    %%%%now do a filtering to select those with durations within [mindur maxdur]
    alldur = endindex-startindex+1; iii = find( (alldur>=groupspan) & (alldur<=maxdur) );
    startindex = startindex(iii); endindex = endindex(iii); %disp('****duraton filtering: '); disp(numel(startindex));
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
       startindex = sort(startindex); endindex = sort(endindex); %disp('****gap filtering: '); disp(numel(startindex));
       %%%%%%%compute parameters in the events
       hvs = findamp(dat, timestamp, startindex, endindex, pth, minPeakNum, sth, peakmode);
       %disp([startindex(1:100)' endindex(1:100)']);
       %disp('****peak amp/number filtering: '); disp(numel(hvs.startT)); disp('&&&&&end');
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


function hvs = findamp(dat, timestamp, startindex, endindex, pth, minPeakNum, sth, peakmode)
%%individual event indicated by [startindex, endindex], indices across startthreshold
%%peaktime, amplitude = maximum values within [startindex endindex] (minimum for slow wave detection)
%%starttime = index in [startindex-groupspan startindex] but cross the startthreshold
%%endtime = index in [endindex endindex+groupspan]but cross the startthreshold
%%area = integrate through [tarttime endtime]
hvs.startT = []; hvs.peakT = []; hvs.endT = []; hvs.dur = []; hvs.amp = []; hvs.eng = []; 
hvs.minTime = []; hvs.maxTime = []; hvs.minAmp = []; hvs.maxAmp = []; hvs.freq = [];
nevent = 0;
for (i = 1:numel(startindex))
     datnow = dat(startindex(i):endindex(i)); timenow = timestamp(startindex(i):endindex(i));
     if strncmpi(peakmode, 'pos', 3)
        [peaknow, peakindex] = max(datnow); %%if either positive crossing
     elseif strncmpi(peakmode, 'neg', 3)
        [peaknow, peakindex] = min(datnow); %%if either negative crossing
     elseif strncmpi(peakmode, 'abs', 3)
        [peaknow, peakindex] = max(abs(datnow)); %%if either positive or negative crossing
     else
         peaknow = NaN; peakindex = NaN;
     end
     if (abs(peaknow) >= pth)
         [minindex, maxindex] = FindLocal(datnow); peakVal = abs(datnow(minindex)); %%%%!!!THis is different from ripples
         if (numel(find(peakVal>=pth))>= minPeakNum)
             nevent = nevent + 1;
             hvs.amp(nevent) = peaknow; hvs.peakT(nevent) = timenow(peakindex);
             hvs.startT(nevent) = timestamp(startindex(i)); hvs.endT(nevent) = timestamp(endindex(i));
             hvs.dur(nevent) = hvs.endT(nevent) - hvs.startT(nevent);
             hvs.freq(nevent) = numel(find(peakVal>=sth))/hvs.dur(nevent);
             hvs.eng(nevent) = mean(abs(datnow))* hvs.dur(nevent);
             hvs.minTime = [hvs.minTime; timenow(minindex)]; hvs.maxTime = [hvs.maxTime; timenow(maxindex)];
             hvs.minAmp = [hvs.minAmp; datnow(minindex)]; hvs.maxAmp = [hvs.maxAmp; datnow(maxindex)];
         end
     end
end

function hvs = adjustminmaxtimes(hvs, dat, timestamp, sth)
iip = find(hvs.minAmp <= -sth); iiq = find(hvs.maxAmp >= sth);
negmanV = mean(abs(hvs.minAmp(iip))); posmanV = mean(abs(hvs.maxAmp(iiq)));
%if (numel(iip) >= numel(iiq)) %%%if HV spikes are negative
if negmanV >= posmanV
    kkp = iip;
    hvs.minTime = hvs.minTime(kkp); hvs.minAmp = hvs.minAmp(kkp);
    maxAmp = NaN*ones(size(hvs.minAmp)); maxTime = NaN*ones(size(hvs.minTime));
    if (~isempty(kkp))
        for (ttp = 1:numel(kkp)-1)
            kkq = find( (hvs.maxTime>hvs.minTime(ttp)) & (hvs.maxTime<min([hvs.minTime(ttp+1) hvs.minTime(ttp)+0.1])) );
            if (~isempty(kkq))
                timenow = hvs.maxTime(kkq); ampfine = hvs.maxAmp(kkq);
                %[ttt, ddd] = min(timenow-hvs.minTime(ttp)); maxAmp(ttp) = ampfine(ddd); maxTime(ttp) = timenow(ddd); %%%%closest maxima
                [ddd, ttt] = max(ampfine); maxAmp(ttp) = ddd; maxTime(ttp) = timenow(ttt); %%%%highest point
            end
        end
        kkq = find( (hvs.maxTime>hvs.minTime(numel(kkp))) & (hvs.maxTime<hvs.minTime(numel(kkp))+0.1) );
        if (~isempty(kkq))
            timenow = hvs.maxTime(kkq); ampfine = hvs.maxAmp(kkq);
            %[ttt, ddd] = min(timenow-hvs.minTime(numel(kkp))); maxAmp(numel(kkp)) = ampfine(ddd); maxTime(numel(kkp)) = timenow(ddd); %%%%closest maxima
            [ddd, ttt] = max(ampfine); maxAmp(numel(kkp)) = ddd; maxTime(numel(kkp)) = timenow(ttt);
        end
    end
    hvs.maxAmp = maxAmp; hvs.maxTime = maxTime;
else %%%if HV spikes are positive
    kkp = iiq;
    hvs.maxTime = hvs.maxTime(kkp); hvs.maxAmp = hvs.maxAmp(kkp);
    minAmp = NaN*ones(size(hvs.maxAmp)); minTime = NaN*ones(size(hvs.maxTime));
    if (~isempty(kkp))
        for (ttp = 1:numel(kkp)-1)
            kkq = find( (hvs.minTime>hvs.maxTime(ttp)) & (hvs.minTime<min([hvs.maxTime(ttp+1) hvs.maxTime(ttp)+0.1])) );
            if (~isempty(kkq))
                timenow = hvs.minTime(kkq); ampfine = hvs.minAmp(kkq);
                %[ttt, ddd] = min(timenow-hvs.maxTime(ttp)); minAmp(ttp) = ampfine(ddd); minTime(ttp) = timenow(ddd); %%%%closest maxima
                [ddd, ttt] = min(ampfine); minAmp(ttp) = ddd; minTime(ttp) = timenow(ttt); %%%%highest point
            end
        end
        kkq = find( (hvs.minTime>hvs.maxTime(numel(kkp))) & (hvs.minTime<hvs.maxTime(numel(kkp))+0.1) );
        if (~isempty(kkq))
            timenow = hvs.minTime(kkq); ampfine = hvs.minAmp(kkq);
            %[ttt, ddd] = min(timenow-hvs.maxTime(numel(kkp))); minAmp(numel(kkp)) = ampfine(ddd); minTime(numel(kkp)) = timenow(ddd); %%%%closest maxima
            [ddd, ttt] = min(ampfine); minAmp(numel(kkp)) = ddd; minTime(numel(kkp)) = timenow(ttt);
        end
    end
    hvs.minAmp = minAmp; hvs.minTime = minTime; 
end
% function runstd = findrunstd(eeg, eegdata, finaldir)
% runstd = []; dat = [];
% %%%look for all the run events and select the data durign these run events
% iii = find( strcmp(eeg.general.finaldir, finaldir) & (strcmp(eeg.parm.band, 'hvs')) ); % &(strcmp(eeg.parm.sessType, 'linear')) );
% for (k = 1:numel(iii))
%     i = iii(k);
%     evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i};
%     tt = find( strcmp(evType, 'run'));
%     if (~isempty(tt))
%         [timestamp, datnow, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
%         datnow = datnow-mean(datnow); datok = [];
%         for (j = 1:numel(tt))
%             for (k = 1:numel(evTime{tt(j)}.start))
%                 datok = [datok; datnow( (timestamp>=evTime{tt(j)}.start(k)) & (timestamp<=evTime{tt(j)}.ent(k)) )];
%             end
%         end
%         dat = [dat; datok]; datnow = [];
%     end
% end
% if (~isempty(dat))
%     runstd = std(dat);
% end

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

function runstd = findrunsessstd(eeg, eegdata, finaldir, bandinfix, recarea)
runstd = []; dat = [];
%%%look for all the run events and select the data durign these run events
iii = find( strcmp(eeg.general.finaldir, finaldir) & strcmp(eeg.parm.band, bandinfix) & strcmp(eeg.general.recarea, recarea) );
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
[filepath, filename] = GetAllFile(posdir, '', '.pwr'); 
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
    iii = find(pwrnow >=40); state(iii) = 2*ones(size(iii)); %%%%running state: >50 fixed threshold in pixels per second
    %jjj = find(pwrnow <=30); state(jjj) = zeros(size(jjj)); %%%%stopp state: <30 fixed threshold in pixels per second
    [stateStart, stateEnd] = findstateepisode(state, 1, numel(state), 2);
    runstart = timenow(stateStart); runend = timenow(stateEnd);
    iii = find(runend-runstart>2); runstart = runstart(iii); runend = runend(iii); %%%select events longer than 2s
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

function maxgap = resetmaxgap(band, eegname)
%%% if use this function for non-hvs detection, need to reset groupspan
%%%% always set as slightly shorter than ~2 cycles
maxgap = 0.1; %default
switch band
    case 'gamma'
        if contains(eegname, 'EMG') %%% EMG filtered through slow gamma ~=40Hz
             maxgap = 0.05;             
        elseif contains(eegname, 'slow') %%% slow gamma ~=30Hz
             maxgap = 0.07;
        elseif contains(eegname, 'high') %%% high gamma ~=100Hz
             maxgap = 0.02;
        else   %%%%regular gamma ~=50Hz
             maxgap = 0.04;
        end
    case 'beta' %%%% 20Hz
         maxgap = 0.1;
    case 'ripple' %%%%~180Hz
         maxgap = 0.012;
    case 'theta' %%%% ~8Hz
         maxgap = 0.25;
end









