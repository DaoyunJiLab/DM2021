function [eeg, eegdata] = DataManager_FindRippleProp(eeg, eegdata, eegind, vv)
%%%detect ripple events and compute their parameters
neeg = numel(eeg.general.eegfile);
%%%the following variables are assigned
% if (~isfield(eegdata, 'ripple')) eegdata.ripple = []; end
% if (~isfield(eegdata.ripple, 'allripples')) eegdata.ripple.allripples = cell(1, neeg); end %%all computed results are here
%%%%displayed variables
if (~isfield(eeg, 'ripple')) eeg.ripple = []; end
if (~isfield(eeg.ripple, 'sessNum')) eeg.ripple.sessNum = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessRate')) eeg.ripple.sessRate = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessMeanAmp')) eeg.ripple.sessMeanAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessMeanDur')) eeg.ripple.sessMeanDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessMeanEng')) eeg.ripple.sessMeanEng = cell(1, neeg); end

if (~isfield(eeg.ripple, 'runNum')) eeg.ripple.runNum = cell(1, neeg); end
if (~isfield(eeg.ripple, 'runRate')) eeg.ripple.runRate = cell(1, neeg); end
if (~isfield(eeg.ripple, 'runMeanAmp')) eeg.ripple.runMeanAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'runMeanDur')) eeg.ripple.runMeanDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'runMeanEng')) eeg.ripple.runMeanEng = cell(1, neeg); end

if (~isfield(eeg.ripple, 'stopNum')) eeg.ripple.stopNum = cell(1, neeg); end
if (~isfield(eeg.ripple, 'stopRate')) eeg.ripple.stopRate = cell(1, neeg); end
if (~isfield(eeg.ripple, 'stopMeanAmp')) eeg.ripple.stopMeanAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'stopMeanDur')) eeg.ripple.stopMeanDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'stopMeanEng')) eeg.ripple.stopMeanEng = cell(1, neeg); end

if (~isfield(eeg.ripple, 'swsNum')) eeg.ripple.swsNum = cell(1, neeg); end
if (~isfield(eeg.ripple, 'swsRate')) eeg.ripple.swsRate = cell(1, neeg); end
if (~isfield(eeg.ripple, 'swsMeanAmp')) eeg.ripple.swsMeanAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'swsMeanDur')) eeg.ripple.swsMeanDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'swsMeanEng')) eeg.ripple.swsMeanEng = cell(1, neeg); end

if (~isfield(eeg.ripple, 'remNum')) eeg.ripple.remNum = cell(1, neeg); end
if (~isfield(eeg.ripple, 'remRate')) eeg.ripple.remRate = cell(1, neeg); end
if (~isfield(eeg.ripple, 'remMeanAmp')) eeg.ripple.remMeanAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'remMeanDur')) eeg.ripple.remMeanDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'remMeanEng')) eeg.ripple.remMeanEng = cell(1, neeg); end

if (~isfield(eeg.ripple, 'sessStartT')) eeg.ripple.sessStartT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessEndT')) eeg.ripple.sessEndT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessPeakT')) eeg.ripple.sessPeakT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessDur')) eeg.ripple.sessDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessAmp')) eeg.ripple.sessAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'sessEng')) eeg.ripple.sessEng = cell(1, neeg); end

if (~isfield(eeg.ripple, 'runStartT')) eeg.ripple.runStartT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'runEndT')) eeg.ripple.runEndT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'runPeakT')) eeg.ripple.runPeakT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'runDur')) eeg.ripple.runDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'runAmp')) eeg.ripple.runAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'runEng')) eeg.ripple.runEng = cell(1, neeg); end

if (~isfield(eeg.ripple, 'stopStartT')) eeg.ripple.stopStartT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'stopEndT')) eeg.ripple.stopEndT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'stopPeakT')) eeg.ripple.stopPeakT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'stopDur')) eeg.ripple.stopDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'stopAmp')) eeg.ripple.stopAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'stopEng')) eeg.ripple.stopEng = cell(1, neeg); end

if (~isfield(eeg.ripple, 'swsStartT')) eeg.ripple.swsStartT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'swsEndT')) eeg.ripple.swsEndT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'swsPeakT')) eeg.ripple.swsPeakT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'swsDur')) eeg.ripple.swsDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'swsAmp')) eeg.ripple.swsAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'swsEng')) eeg.ripple.swsEng = cell(1, neeg); end

if (~isfield(eeg.ripple, 'remStartT')) eeg.ripple.remStartT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'remEndT')) eeg.ripple.remEndT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'remPeakT')) eeg.ripple.remPeakT = cell(1, neeg); end
if (~isfield(eeg.ripple, 'remDur')) eeg.ripple.remDur = cell(1, neeg); end
if (~isfield(eeg.ripple, 'remAmp')) eeg.ripple.remAmp = cell(1, neeg); end
if (~isfield(eeg.ripple, 'remEng')) eeg.ripple.remEng = cell(1, neeg); end

if (~isfield(eeg.ripple, 'sessMaxAmp')) eeg.ripple.sessMaxAmp = cell(1, neeg); end 
if (~isfield(eeg.ripple, 'sessMaxDur')) eeg.ripple.sessMaxDur = cell(1, neeg); end  %top 10% ripples duration
if (~isfield(eeg.ripple, 'sessMaxEng')) eeg.ripple.sessMaxEng = cell(1, neeg); end 
if (~isfield(eeg.ripple, 'swsMaxAmp')) eeg.ripple.swsMaxAmp = cell(1, neeg); end 
if (~isfield(eeg.ripple, 'swsMaxDur')) eeg.ripple.swsMaxDur = cell(1, neeg); end  %top 10% ripples duration
if (~isfield(eeg.ripple, 'swsMaxEng')) eeg.ripple.swsMaxEng = cell(1, neeg); end 
if (~isfield(eeg.ripple, 'stopMaxAmp')) eeg.ripple.stopMaxAmp = cell(1, neeg); end 
if (~isfield(eeg.ripple, 'stopMaxDur')) eeg.ripple.stopMaxDur = cell(1, neeg); end  %top 10% ripples duration
if (~isfield(eeg.ripple, 'stopMaxEng')) eeg.ripple.stopMaxEng = cell(1, neeg); end 

if (~isfield(eeg.ripple, 'startThreshold')) eeg.ripple.startThreshold = cell(1, neeg); end
if (~isfield(eeg.ripple, 'peakThreshold')) eeg.ripple.peakThreshold = cell(1, neeg); end

for (iiik = 1:numel(eegind))
i = eegind(iiik);
if strcmp(eeg.parm.band{i}, 'ripple') %%only do this for filtered EEG traces
    disp(['--------> ripple analysis: ', eeg.general.eegfile{i}]);
    ripnorm = eeg.parm.rippNorm{i}; pth = eeg.parm.rippPeakThres(i); sth = eeg.parm.rippStartThres(i); 
    maxgap = eeg.parm.rippMaxGap(i)*eeg.general.freq{i}; %%%now mingap in datapoints
    mindur = eeg.parm.rippMinDur(i)*eeg.general.freq{i};
    maxdur = eeg.parm.rippMaxDur(i)*eeg.general.freq{i};
    %%%get eeg data
    [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
    dat = dat-mean(dat);
    %%%find parameters for ripple detection
    if (strcmp(ripnorm, 'self'))
        pth = pth*std(dat); sth = sth*std(dat); disp(['---------------> self std: ', num2str(std(dat))]);
    elseif (strcmp(ripnorm, 'run'))
        runstd = findrunstd(eeg, eegdata, eeg.general.finaldir{i}); 
        if (~isempty(runstd))
           pth = pth*runstd; sth = sth*runstd; disp(['---------------> run std: ', num2str(runstd)]);
        else
           msgbox('Warning: run eeg ripple data not found; use self std instead as threshold unit');
           pth = pth*std(dat); sth = sth*std(dat); disp(['---------------> self std: ', num2str(std(dat))]);
        end
    end
    eeg.ripple.startThreshold{i} = sth; eeg.ripple.peakThreshold{i} = pth;
    %%%detect ripples
    ripp = findallripples(timestamp, dat, pth, sth, maxgap, mindur, maxdur); %%ripp.startT(i), .peakT, .endT, .dur, .amp
    %%%session variables
    if (~isempty(eeg.general.sesslength{i})) & (eeg.general.sesslength{i}>0)
    eeg.ripple.sessNum{i} = numel(ripp.startT); eeg.ripple.sessRate{i} = numel(ripp.startT)/eeg.general.sesslength{i};
    eeg.ripple.sessMeanAmp{i} = mean(ripp.amp); eeg.ripple.sessMeanDur{i} = mean(ripp.dur); eeg.ripple.sessMeanEng{i} = mean(ripp.eng);
    eeg.ripple.sessStartT{i} = ripp.startT; eeg.ripple.sessPeakT{i} = ripp.peakT; eeg.ripple.sessEndT{i} = ripp.endT; 
    eeg.ripple.sessAmp{i} = ripp.amp; eeg.ripple.sessDur{i} = ripp.dur; eeg.ripple.sessEng{i} = ripp.eng;
    %%do max ripples: top 10% power
       ttt = prctile(eeg.ripple.sessEng{i}, 90); iii = find(eeg.ripple.sessEng{i}>= ttt);
       eeg.ripple.sessMaxAmp{i} = eeg.ripple.sessAmp{i}(iii); eeg.ripple.sessMaxDur{i} = eeg.ripple.sessDur{i}(iii); 
       eeg.ripple.sessMaxEng{i} = eeg.ripple.sessEng{i}(iii);
    %%%event variables
    evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i};
    runid = []; stopid = [];  swsid = []; remid = []; runleng = 0; stopleng = 0;  swsleng = 0; remleng = 0; 
    for (j = 1:numel(evTime))
        st = evTime{j}.start; et = evTime{j}.ent; nev = numel(st); iii = []; leng = 0;
        for (ttt = 1:nev)
            ik = find( (ripp.peakT>=st(ttt)) & (ripp.peakT<et(ttt)) ); iii = union(iii, ik); ik = []; leng = leng + et(ttt)-st(ttt);
        end
        if (strcmp(evType{j}, 'run')) runid = union(runid, iii); runleng = runleng+ leng; end
        if (strcmp(evType{j}, 'stop')) stopid = union(stopid, iii); stopleng = stopleng + leng; end
        if (strcmp(evType{j}, 'sws')) swsid = union(swsid, iii); swsleng = swsleng+leng; end
        if (strcmp(evType{j}, 'rem')) remid = union(remid, iii); remleng = remleng + leng; end
    end    
    %%%%%%run variables
    if (runleng > 0)
       eeg.ripple.runNum{i} = numel(ripp.startT(runid)); eeg.ripple.runRate{i} = numel(ripp.startT(runid))/runleng;
       eeg.ripple.runMeanAmp{i} = mean(ripp.amp(runid)); eeg.ripple.runMeanDur{i} = mean(ripp.dur(runid));
       eeg.ripple.runMeanEng{i} = mean(ripp.eng(runid));
       eeg.ripple.runStartT{i} = ripp.startT(runid); eeg.ripple.runPeakT{i} = ripp.peakT(runid); eeg.ripple.runEndT{i} = ripp.endT(runid); 
       eeg.ripple.runAmp{i} = ripp.amp(runid); eeg.ripple.runDur{i} = ripp.dur(runid); eeg.ripple.runEng{i} = ripp.eng;
    end
    if (stopleng > 0)
       eeg.ripple.stopNum{i} = numel(ripp.startT(stopid)); eeg.ripple.stopRate{i} = numel(ripp.startT(stopid))/stopleng;
       eeg.ripple.stopMeanAmp{i} = mean(ripp.amp(stopid)); eeg.ripple.stopMeanDur{i} = mean(ripp.dur(stopid));
       eeg.ripple.stopMeanEng{i} = mean(ripp.eng(stopid));
       eeg.ripple.stopStartT{i} = ripp.startT(stopid); eeg.ripple.stopPeakT{i} = ripp.peakT(stopid); eeg.ripple.stopEndT{i} = ripp.endT(stopid); 
       eeg.ripple.stopAmp{i} = ripp.amp(stopid); eeg.ripple.stopDur{i} = ripp.dur(stopid); eeg.ripple.stopEng{i} = ripp.eng;
       %%do max ripples: top 10%
       ttt = prctile(eeg.ripple.stopEng{i}, 90); iii = find(eeg.ripple.stopEng{i}>= ttt);
       eeg.ripple.stopMaxAmp{i} = eeg.ripple.stopAmp{i}(iii); eeg.ripple.stopMaxDur{i} = eeg.ripple.stopDur{i}(iii); 
       eeg.ripple.stopMaxEng{i} = eeg.ripple.stopEng{i}(iii);
    end
    if (swsleng > 0)
       eeg.ripple.swsNum{i} = numel(ripp.startT(swsid)); eeg.ripple.swsRate{i} = numel(ripp.startT(swsid))/swsleng;
       eeg.ripple.swsMeanAmp{i} = mean(ripp.amp(swsid)); eeg.ripple.swsMeanDur{i} = mean(ripp.dur(swsid));
       eeg.ripple.swsMeanEng{i} = mean(ripp.eng(swsid));
       eeg.ripple.swsStartT{i} = ripp.startT(swsid); eeg.ripple.swsPeakT{i} = ripp.peakT(swsid); eeg.ripple.swsEndT{i} = ripp.endT(swsid); 
       eeg.ripple.swsAmp{i} = ripp.amp(swsid); eeg.ripple.swsDur{i} = ripp.dur(swsid); eeg.ripple.swsEng{i} = ripp.eng(swsid);
       %%do max ripples: top 10%
       ttt = prctile(eeg.ripple.swsEng{i}, 90); iii = find(eeg.ripple.swsEng{i}>= ttt);
       eeg.ripple.swsMaxAmp{i} = eeg.ripple.swsAmp{i}(iii); eeg.ripple.swsMaxDur{i} = eeg.ripple.swsDur{i}(iii); 
       eeg.ripple.swsMaxEng{i} = eeg.ripple.swsEng{i}(iii);
    end  
    if (remleng > 0)
       eeg.ripple.remNum{i} = numel(ripp.startT(remid)); eeg.ripple.remRate{i} = numel(ripp.startT(remid))/remleng;
       eeg.ripple.remMeanAmp{i} = mean(ripp.amp(remid)); eeg.ripple.remMeanDur{i} = mean(ripp.dur(remid));
       eeg.ripple.remMeanEng{i} = mean(ripp.eng(remid));
       eeg.ripple.remStartT{i} = ripp.startT(remid); eeg.ripple.remPeakT{i} = ripp.peakT(remid); eeg.ripple.remEndT{i} = ripp.endT(remid); 
       eeg.ripple.remAmp{i} = ripp.amp(remid); eeg.ripple.remDur{i} = ripp.dur(remid); eeg.ripple.remEng{i} = ripp.eng(remid);
    end  
    end
end
end

function ripp = findallripples(timestamp, dat, pth, sth, maxgap, mindur, maxdur) %%ripp.startT(i), .peakT, .endT, .dur, .amp
ripp.startT = []; ripp.peakT = []; ripp.endT = []; ripp.dur = []; ripp.amp = []; ripp.eng = [];
index = find(abs(dat) >= sth); %all threshold crossing for both positive crossing and negative crossing
if (~isempty(index))
    [startindex, endindex] = findeventindex(index, maxgap); %%start/end indices for start threshold crossing
    %%%%now do a filtering to select those with durations within [mindur maxdur]
    alldur = endindex-startindex+1; iii = find( (alldur>=mindur) & (alldur<=maxdur) );
    startindex = startindex(iii); endindex = endindex(iii);
    ripp = findamp(dat, timestamp, startindex, endindex, pth);
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

function ripp = findamp(dat, timestamp, startindex, endindex, pth)
%%individual event indicated by [startindex, endindex], indices across startthreshold
%%peaktime, amplitude = maximum values within [startindex endindex] (minimum for slow wave detection)
%%starttime = index in [startindex-groupspan startindex] but cross the startthreshold
%%endtime = index in [endindex endindex+groupspan]but cross the startthreshold
%%area = integrate through [tarttime endtime]
ripp.startT = []; ripp.peakT = []; ripp.endT = []; ripp.dur = []; ripp.amp = []; ripp.eng = [];
nevent = 0;
for (i = 1:numel(startindex))
     datnow = dat(startindex(i):endindex(i)); timenow = timestamp(startindex(i):endindex(i));
     [abspeak, peakindex] = max(abs(datnow)); %%if either positive crossing or negative crossing
     if (abspeak >= pth)
         nevent = nevent + 1;
         ripp.amp(nevent) = abspeak; ripp.peakT(nevent) = timenow(peakindex);
         ripp.startT(nevent) = timestamp(startindex(i)); ripp.endT(nevent) = timestamp(endindex(i));
         ripp.dur(nevent) = ripp.endT(nevent) - ripp.startT(nevent);
         ripp.eng(nevent) = mean(abs(datnow))* ripp.dur(nevent);
     end
end

function runstd = findrunstd(eeg, eegdata, finaldir)
runstd = []; dat = [];
%%%look for all the run events and select the data durign these run events
iii = find( strcmp(eeg.general.finaldir, finaldir) & (strcmp(eeg.parm.band, 'ripple')) ); % &(strcmp(eeg.parm.sessType, 'linear')) );
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



