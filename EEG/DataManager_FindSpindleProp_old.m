function [eeg, eegdata] = DataManager_FindSpindleProp(eeg, eegdata, eegind, vv)

%%%%first add filtered spindle-band EEG traces
iii = find( strcmp(eeg.parm.band, 'spindle') );
SS = questdlg(['Found ', num2str(numel(iii)), ' spindle files. Filter broad band files?']);
ok = 1;
if (strcmp(SS, 'Yes'))
    %%%%%read spindle filer
    [MCroot, MCname, DAname, DEname] = CurrentVersion;
    filterfile = fullfile(MCroot, 'Filters', 'EEGspindleresample.mat'); %filter = EEG spindle, sampling freq needs to be ~200Hz
    s = load(filterfile); num = s.Num; den = s.Den; s = []; %%denumerator 
    %%%%%filter the broad-band EEG traces
    iii = find( strcmp(eeg.parm.band, 'broad') ); neeg = numel(eeg.general.eegfile);
    for (i = 1:numel(iii))
        eegfilenow = eeg.general.eegfile{iii(i)}; buffersize = eeg.parm.buffersize(iii(i)); infix = 'spindle';
        disp(['--------> filter file: ', eegfilenow]);
        EEG_ResampleFilterData(1, 10, num, den, eegfilenow, buffersize, infix, vv);
        [eeg, eegdata] = appendfile(eeg, eegdata, eegfilenow, infix, iii(i));
    end
    eegind = 1:numel(eeg.general.eegfile);
elseif (strcmp(SS, 'Cancel'))
    ok = 0; disp(['--------------> action cancelles']);
end
if ok
   [eeg, eegdata] = DoSpindleAnalysis(eeg, eegdata, eegind, vv);
end

function [eeg, eegdata] = DoSpindleAnalysis(eeg, eegdata, eegind, vv)
%%%detect spindle events and compute their parameters
neeg = numel(eeg.general.eegfile);
%%%the following variables are assigned
% if (~isfield(eegdata, 'spindle')) eegdata.spindle = []; end
% if (~isfield(eegdata.spindle, 'allspindles')) eegdata.spindle.allspindles = cell(1, neeg); end %%all computed results are here
%%%%displayed variables
if (~isfield(eeg, 'spindle')) eeg.spindle = []; end
if (~isfield(eeg.spindle, 'sessNum')) eeg.spindle.sessNum = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessRate')) eeg.spindle.sessRate = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanAmp')) eeg.spindle.sessMeanAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanAmp')) eeg.spindle.sessMeanFreq = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanDur')) eeg.spindle.sessMeanDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanEng')) eeg.spindle.sessMeanEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'runNum')) eeg.spindle.runNum = cell(1, neeg); end
if (~isfield(eeg.spindle, 'runRate')) eeg.spindle.runRate = cell(1, neeg); end
if (~isfield(eeg.spindle, 'runMeanAmp')) eeg.spindle.runMeanAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanAmp')) eeg.spindle.runMeanFreq = cell(1, neeg); end
if (~isfield(eeg.spindle, 'runMeanDur')) eeg.spindle.runMeanDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'runMeanEng')) eeg.spindle.runMeanEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'stopNum')) eeg.spindle.stopNum = cell(1, neeg); end
if (~isfield(eeg.spindle, 'stopRate')) eeg.spindle.stopRate = cell(1, neeg); end
if (~isfield(eeg.spindle, 'stopMeanAmp')) eeg.spindle.stopMeanAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanAmp')) eeg.spindle.stopMeanFreq = cell(1, neeg); end
if (~isfield(eeg.spindle, 'stopMeanDur')) eeg.spindle.stopMeanDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'stopMeanEng')) eeg.spindle.stopMeanEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'swsNum')) eeg.spindle.swsNum = cell(1, neeg); end
if (~isfield(eeg.spindle, 'swsRate')) eeg.spindle.swsRate = cell(1, neeg); end
if (~isfield(eeg.spindle, 'swsMeanAmp')) eeg.spindle.swsMeanAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanAmp')) eeg.spindle.swsMeanFreq = cell(1, neeg); end
if (~isfield(eeg.spindle, 'swsMeanDur')) eeg.spindle.swsMeanDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'swsMeanEng')) eeg.spindle.swsMeanEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'remNum')) eeg.spindle.remNum = cell(1, neeg); end
if (~isfield(eeg.spindle, 'remRate')) eeg.spindle.remRate = cell(1, neeg); end
if (~isfield(eeg.spindle, 'remMeanAmp')) eeg.spindle.remMeanAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanAmp')) eeg.spindle.remMeanFreq = cell(1, neeg); end
if (~isfield(eeg.spindle, 'remMeanDur')) eeg.spindle.remMeanDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'remMeanEng')) eeg.spindle.remMeanEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'minTime')) eeg.spindle.minTime = cell(1, neeg); end
if (~isfield(eeg.spindle, 'maxTime')) eeg.spindle.maxTime = cell(1, neeg); end
if (~isfield(eeg.spindle, 'minAmp')) eeg.spindle.minAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'maxAmp')) eeg.spindle.maxAmp = cell(1, neeg); end

if (~isfield(eeg.spindle, 'sessStartT')) eeg.spindle.sessStartT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessEndT')) eeg.spindle.sessEndT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessPeakT')) eeg.spindle.sessPeakT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessDur')) eeg.spindle.sessDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessAmp')) eeg.spindle.sessAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessEng')) eeg.spindle.sessEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'runStartT')) eeg.spindle.runStartT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'runEndT')) eeg.spindle.runEndT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'runPeakT')) eeg.spindle.runPeakT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'runDur')) eeg.spindle.runDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'runAmp')) eeg.spindle.runAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'runEng')) eeg.spindle.runEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'stopStartT')) eeg.spindle.stopStartT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'stopEndT')) eeg.spindle.stopEndT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'stopPeakT')) eeg.spindle.stopPeakT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'stopDur')) eeg.spindle.stopDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'stopAmp')) eeg.spindle.stopAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'stopEng')) eeg.spindle.stopEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'swsStartT')) eeg.spindle.swsStartT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'swsEndT')) eeg.spindle.swsEndT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'swsPeakT')) eeg.spindle.swsPeakT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'swsDur')) eeg.spindle.swsDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'swsAmp')) eeg.spindle.swsAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'swsEng')) eeg.spindle.swsEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'remStartT')) eeg.spindle.remStartT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'remEndT')) eeg.spindle.remEndT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'remPeakT')) eeg.spindle.remPeakT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'remDur')) eeg.spindle.remDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'remAmp')) eeg.spindle.remAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'remEng')) eeg.spindle.remEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'startThreshold')) eeg.spindle.startThreshold = cell(1, neeg); end
if (~isfield(eeg.spindle, 'peakThreshold')) eeg.spindle.peakThreshold = cell(1, neeg); end

%if (~isfield(eeg.spindle, 'swsMaxAmp')) eeg.spindle.swsMaxAmp = cell(1, neeg); end 
%if (~isfield(eeg.spindle, 'swsMaxDur')) eeg.spindle.swsMaxDur = cell(1, neeg); end  %top 10% spindles duration
%if (~isfield(eeg.spindle, 'swsMaxEng')) eeg.spindle.swsMaxEng = cell(1, neeg); end 
%if (~isfield(eeg.spindle, 'stopMaxAmp')) eeg.spindle.stopMaxAmp = cell(1, neeg); end 
%if (~isfield(eeg.spindle, 'stopMaxDur')) eeg.spindle.stopMaxDur = cell(1, neeg); end  %top 10% spindles duration
%if (~isfield(eeg.spindle, 'stopMaxEng')) eeg.spindle.stopMaxEng = cell(1, neeg); end 

for (iiik = 1:numel(eegind))
i = eegind(iiik);
if strcmp(eeg.parm.band{i}, 'spindle') %%only do this for filtered EEG traces
    disp(['--------> spindle analysis: ', eeg.general.eegfile{i}]);
    ripnorm = eeg.parm.spinNorm{i}; pth = eeg.parm.spinPeakThres(i); sth = eeg.parm.spinStartThres(i); 
    maxgap = eeg.parm.spinMaxGap(i)*eeg.general.freq{i}; %%%now mingap in datapoints
    mindur = eeg.parm.spinMinDur(i)*eeg.general.freq{i};
    maxdur = eeg.parm.spinMaxDur(i)*eeg.general.freq{i};
    %%%get eeg data
    [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
    dat = dat-mean(dat);
    %%%find parameters for spindle detection
    if (strcmp(ripnorm, 'self'))
        pth = pth*std(dat); sth = sth*std(dat); disp(['---------------> self std: ', num2str(std(dat))]);
    elseif (strcmp(ripnorm, 'run'))
        runstd = findrunstd(eeg, eegdata, eeg.general.finaldir{i}); 
        if (~isempty(runstd))
           pth = pth*runstd; sth = sth*runstd; disp(['---------------> run std: ', num2str(runstd)]);
        else
           msgbox('Warning: run eeg spindle data not found; use self std instead as threshold unit');
           pth = pth*std(dat); sth = sth*std(dat); disp(['---------------> self std: ', num2str(std(dat))]);
        end
    end
    eeg.spindle.startThreshold{i} = sth; eeg.spindle.peakThreshold{i} = pth;
    %%%detect spindles
    spin = findallspindles(timestamp, dat, pth, sth, maxgap, mindur, maxdur); %%spin.startT(i), .peakT, .endT, .dur, .amp
    %%%session variables
    if (~isempty(eeg.general.sesslength{i})) & (eeg.general.sesslength{i}>0)
    eeg.spindle.sessNum{i} = numel(spin.startT); eeg.spindle.sessRate{i} = numel(spin.startT)/eeg.general.sesslength{i};
    eeg.spindle.sessMeanAmp{i} = mean(spin.amp); eeg.spindle.sessMeanDur{i} = mean(spin.dur); eeg.spindle.sessMeanEng{i} = mean(spin.eng);
    eeg.spindle.sessStartT{i} = spin.startT; eeg.spindle.sessPeakT{i} = spin.peakT; eeg.spindle.sessEndT{i} = spin.endT; 
    eeg.spindle.sessAmp{i} = spin.amp; eeg.spindle.sessDur{i} = spin.dur; eeg.spindle.sessEng{i} = spin.eng;
    eeg.spindle.minTime{i} = spin.minTime; eeg.spindle.maxTime{i} = spin.maxTime;
    eeg.spindle.minAmp{i} = spin.minAmp; eeg.spindle.maxAmp{i} = spin.maxAmp;
    %%%event variables
    evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i};
    runid = []; stopid = [];  swsid = []; remid = []; runleng = 0; stopleng = 0;  swsleng = 0; remleng = 0; 
    for (j = 1:numel(evTime))
        st = evTime{j}.start; et = evTime{j}.ent; nev = numel(st); iii = []; leng = 0;
        for (ttt = 1:nev)
            ik = find( (spin.peakT>=st(ttt)) & (spin.peakT<et(ttt)) ); iii = union(iii, ik); ik = []; leng = leng + et(ttt)-st(ttt);
        end
        if (strcmp(evType{j}, 'run')) runid = union(runid, iii); runleng = runleng+ leng; end
        if (strcmp(evType{j}, 'stop')) stopid = union(stopid, iii); stopleng = stopleng + leng; end
        if (strcmp(evType{j}, 'sws')) swsid = union(swsid, iii); swsleng = swsleng+leng; end
        if (strcmp(evType{j}, 'rem')) remid = union(remid, iii); remleng = remleng + leng; end
    end    
    %%%%%%run variables
    if (runleng > 0)
       eeg.spindle.runNum{i} = numel(spin.startT(runid)); eeg.spindle.runRate{i} = numel(spin.startT(runid))/runleng;
       eeg.spindle.runMeanAmp{i} = mean(spin.amp(runid)); eeg.spindle.runMeanDur{i} = mean(spin.dur(runid));
       eeg.spindle.runMeanEng{i} = mean(spin.eng(runid));
       eeg.spindle.runStartT{i} = spin.startT(runid); eeg.spindle.runPeakT{i} = spin.peakT(runid); eeg.spindle.runEndT{i} = spin.endT(runid); 
       eeg.spindle.runAmp{i} = spin.amp(runid); eeg.spindle.runDur{i} = spin.dur(runid); eeg.spindle.runEng{i} = spin.eng;
    end
    if (stopleng > 0)
       eeg.spindle.stopNum{i} = numel(spin.startT(stopid)); eeg.spindle.stopRate{i} = numel(spin.startT(stopid))/stopleng;
       eeg.spindle.stopMeanAmp{i} = mean(spin.amp(stopid)); eeg.spindle.stopMeanDur{i} = mean(spin.dur(stopid));
       eeg.spindle.stopMeanEng{i} = mean(spin.eng(stopid));
       eeg.spindle.stopStartT{i} = spin.startT(stopid); eeg.spindle.stopPeakT{i} = spin.peakT(stopid); eeg.spindle.stopEndT{i} = spin.endT(stopid); 
       eeg.spindle.stopAmp{i} = spin.amp(stopid); eeg.spindle.stopDur{i} = spin.dur(stopid); eeg.spindle.stopEng{i} = spin.eng;
       %%do max spindles: top 10%
       ttt = prctile(eeg.spindle.stopAmp{i}, 10); iii = find(eeg.spindle.stopAmp{i}>= ttt);
       eeg.spindle.stopMaxAmp{i} = eeg.spindle.stopAmp{i}(iii); eeg.spindle.stopMaxDur{i} = eeg.spindle.stopDur{i}(iii); 
       eeg.spindle.stopMaxEng{i} = eeg.spindle.stopEng{i}(iii);
    end
    if (swsleng > 0)
       eeg.spindle.swsNum{i} = numel(spin.startT(swsid)); eeg.spindle.swsRate{i} = numel(spin.startT(swsid))/swsleng;
       eeg.spindle.swsMeanAmp{i} = mean(spin.amp(swsid)); eeg.spindle.swsMeanDur{i} = mean(spin.dur(swsid));
       eeg.spindle.swsMeanEng{i} = mean(spin.eng(swsid));
       eeg.spindle.swsStartT{i} = spin.startT(swsid); eeg.spindle.swsPeakT{i} = spin.peakT(swsid); eeg.spindle.swsEndT{i} = spin.endT(swsid); 
       eeg.spindle.swsAmp{i} = spin.amp(swsid); eeg.spindle.swsDur{i} = spin.dur(swsid); eeg.spindle.swsEng{i} = spin.eng;
       %%do max spindles: top 10%
       ttt = prctile(eeg.spindle.swsAmp{i}, 10); iii = find(eeg.spindle.swsAmp{i}>= ttt);
       eeg.spindle.swsMaxAmp{i} = eeg.spindle.swsAmp{i}(iii); eeg.spindle.swsMaxDur{i} = eeg.spindle.swsDur{i}(iii); 
       eeg.spindle.swsMaxEng{i} = eeg.spindle.swsEng{i}(iii);
    end  
    if (remleng > 0)
       eeg.spindle.remNum{i} = numel(spin.startT(remid)); eeg.spindle.remRate{i} = numel(spin.startT(remid))/remleng;
       eeg.spindle.remMeanAmp{i} = mean(spin.amp(remid)); eeg.spindle.remMeanDur{i} = mean(spin.dur(remid));
       eeg.spindle.remMeanEng{i} = mean(spin.eng(remid));
       eeg.spindle.remStartT{i} = spin.startT(remid); eeg.spindle.remPeakT{i} = spin.peakT(remid); eeg.spindle.remEndT{i} = spin.endT(remid); 
       eeg.spindle.remAmp{i} = spin.amp(remid); eeg.spindle.remDur{i} = spin.dur(remid); eeg.spindle.remEng{i} = spin.eng;
    end  
    end
end
end

function spin = findallspindles(timestamp, dat, pth, sth, maxgap, mindur, maxdur) %%spin.startT(i), .peakT, .endT, .dur, .amp
spin.startT = []; spin.peakT = []; spin.endT = []; spin.dur = []; spin.amp = []; spin.eng = []; 
spin.minTime = []; spin.maxTime = []; spin.minAmp = []; spin.maxAmp = [];
index = find(abs(dat) >= sth); %all threshold crossing for both positive crossing and negative crossing
if (~isempty(index))
    [startindex, endindex] = findeventindex(index, maxgap); %%start/end indices for start threshold crossing
    %%%%now do a filtering to select those with durations within [mindur maxdur]
    alldur = endindex-startindex+1; iii = find( (alldur>=mindur) & (alldur<=maxdur) );
    startindex = startindex(iii); endindex = endindex(iii);
    spin = findamp(dat, timestamp, startindex, endindex, pth);
    %disp([startindex(1:100)' endindex(1:100)']);
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

function spin = findamp(dat, timestamp, startindex, endindex, pth)
%%individual event indicated by [startindex, endindex], indices across startthreshold
%%peaktime, amplitude = maximum values within [startindex endindex] (minimum for slow wave detection)
%%starttime = index in [startindex-groupspan startindex] but cross the startthreshold
%%endtime = index in [endindex endindex+groupspan]but cross the startthreshold
%%area = integrate through [tarttime endtime]
spin.startT = []; spin.peakT = []; spin.endT = []; spin.dur = []; spin.amp = []; spin.eng = []; 
spin.minTime = []; spin.maxTime = []; spin.minAmp = []; spin.maxAmp = [];
nevent = 0;
for (i = 1:numel(startindex))
     datnow = dat(startindex(i):endindex(i)); timenow = timestamp(startindex(i):endindex(i));
     [abspeak, peakindex] = max(abs(datnow)); %%if either positive crossing or negative crossing
     if (abspeak >= pth)
         nevent = nevent + 1;
         spin.amp(nevent) = abspeak; spin.peakT(nevent) = timenow(peakindex);
         spin.startT(nevent) = timestamp(startindex(i)); spin.endT(nevent) = timestamp(endindex(i));
         spin.dur(nevent) = spin.endT(nevent) - spin.startT(nevent);
         spin.eng(nevent) = mean(abs(datnow))* spin.dur(nevent);
         [minindex, maxindex] = FindLocal(datnow); 
         spin.minTime = [spin.minTime; timenow(minindex)]; spin.maxTime = [spin.maxTime; timenow(maxindex)];
         spin.minAmp = [spin.minAmp; datnow(minindex)]; spin.maxAmp = [spin.maxAmp; datnow(maxindex)];
     end
end

function runstd = findrunstd(eeg, eegdata, finaldir)
runstd = []; dat = [];
%%%look for all the run events and select the data durign these run events
iii = find( strcmp(eeg.general.finaldir, finaldir) & (strcmp(eeg.parm.band, 'spindle')) ); % &(strcmp(eeg.parm.sessType, 'linear')) );
for (k = 1:numel(iii))
    i = iii(k);
    evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i};
    tt = find( strcmp(evType, 'run'));
    if (~isempty(tt))
        [timestamp, datnow, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
        datnow = datnow-mean(datnow); datok = [];
        for (j = 1:numel(tt))
            for (k = 1:numel(evTime{tt(j)}.start))
                datok = [datok; datnow( (timestamp>=evTime{tt(j)}.start(k)) & (timestamp<=evTime{tt(j)}.ent(k)) )];
            end
        end
        dat = [dat; datok]; datnow = [];
    end
end
if (~isempty(dat))
    runstd = std(dat);
end

function [eeg, eegdata] = appendfile(eeg, eegdata, eegfilenow, infix, index)
nfile = numel(eeg.general.eegfile); [pp, nn, ee] = fileparts(eegfilenow);
ID = strcat(nn, '_', infix); newname = fullfile(pp, strcat(ID, ee));
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




