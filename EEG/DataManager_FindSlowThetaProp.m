function [eeg, eegdata] = DataManager_FindSlowThetaProp(eeg, eegdata, eegind, vv)
ok = 1;
%%%%first add filtered ripple-band EEG traces
ttt = find( strcmp(eeg.parm.band(eegind), 'cluster0') ); 
if numel(ttt) == 0
   SS = questdlg('Found no cluster0 files. Continue?'); ok = 0;
   if (strcmp(SS, 'Yes')) ok=1; end
end
if ok
   [eeg, eegdata] = DoSlowThetaAnalysis(eeg, eegdata, eegind);
end

function [eeg, eegdata] = DoSlowThetaAnalysis(eeg, eegdata, eegind)
%%%detect slowtheta events and compute their parameters
neeg = numel(eeg.general.eegfile);
%%%the following variables are assigned
%%%%displayed variables
if (~isfield(eeg, 'slowtheta')) eeg.slowtheta = []; end
if (~isfield(eeg.slowtheta, 'startThreshold')) eeg.slowtheta.startThreshold = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'troughThreshold')) eeg.slowtheta.troughThreshold = cell(1, neeg); end

if (~isfield(eeg.slowtheta, 'sessNum')) eeg.slowtheta.sessNum = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessRate')) eeg.slowtheta.sessRate = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessTimeRatio')) eeg.slowtheta.sessTimeRatio = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessMeanAmp')) eeg.slowtheta.sessMeanAmp = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessMeanFreq')) eeg.slowtheta.sessMeanFreq = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessMeanDur')) eeg.slowtheta.sessMeanDur = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessMeanEng')) eeg.slowtheta.sessMeanEng = cell(1, neeg); end

if (~isfield(eeg.slowtheta, 'evtNum')) eeg.slowtheta.evtNum = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'evtRate')) eeg.slowtheta.evtRate = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'evtTimeRatio')) eeg.slowtheta.evtTimeRatio = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'evtMeanAmp')) eeg.slowtheta.evtMeanAmp = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'evtMeanFreq')) eeg.slowtheta.evtMeanFreq = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'evtMeanDur')) eeg.slowtheta.evtMeanDur = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'evtMeanEng')) eeg.slowtheta.evtMeanEng = cell(1, neeg); end

if (~isfield(eeg.slowtheta, 'sessStartT')) eeg.slowtheta.sessStartT = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessEndT')) eeg.slowtheta.sessEndT = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessPeakT')) eeg.slowtheta.sessPeakT = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessDur')) eeg.slowtheta.sessDur = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessAmp')) eeg.slowtheta.sessAmp = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessEng')) eeg.slowtheta.sessEng = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'sessFreq')) eeg.slowtheta.sessFreq = cell(1, neeg); end

if (~isfield(eeg.slowtheta, 'minTime')) eeg.slowtheta.minTime = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'maxTime')) eeg.slowtheta.maxTime = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'minAmp')) eeg.slowtheta.minAmp = cell(1, neeg); end
if (~isfield(eeg.slowtheta, 'maxAmp')) eeg.slowtheta.maxAmp = cell(1, neeg); end

for (iiik = 1:numel(eegind))
i = eegind(iiik);
if strcmp(eeg.parm.band{i}, 'cluster0') %%only do this for filtered EEG traces remove this now
    disp(['--------> detecting slowtheta events: ', eeg.general.eegfile{i}]);
    pth = eeg.parm.slowthetaTroughThres(i); sth = eeg.parm.slowthetaStartThres(i); 
    mindur = eeg.parm.slowthetaMinDur(i); %*eeg.general.freq{i};
    maxdur = eeg.parm.slowthetaMaxDur(i); mineventgap = eeg.parm.slowthetaMinEventGap(i);
    %disp(['---------------> groupspan set as <: ', num2str(maxgap), 's']);    
    minPeakNum = eeg.parm.slowthetaMinPeakNum(i); %%%*eeg.general.freq{i}; ---!!!This changed from ripple detection = min number of peaks passing the peak threshold
    %%%here are data to analyze
    if strcmp(eeg.parm.band{i}, 'cluster0')
        timestamp = eegdata.cluster0.timepoint{i}; dat = eegdata.cluster0.cnt{i}; 
        binsize = eeg.parm.cl0Binsize(i);
        mindur = round(mindur/binsize); mineventgap = round(mineventgap/binsize);%%%now maxgap in datapoints 
        maxdur = round(maxdur/binsize);
    else
        [timestamp, dat, ~, ~] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i));% timestamp/point in second, dat/cnt in mV
        mindur = mindur*eeg.general.freq{i}; %%%now maxgap in datapoints 
        mineventgap = mineventgap*eeg.general.freq{i}; %%%now mingap in datapoints 
        maxdur = maxdur*eeg.general.freq{i};
    end
    eeg.slowtheta.startThreshold{i} = sth; eeg.slowtheta.troughThreshold{i} = pth;
    %%%session variables
    if (~isempty(eeg.general.sesslength{i})) && (eeg.general.sesslength{i}>0)
        %%%detect slowthetas
        slowtheta = findallslowthetas(timestamp, dat, pth, sth, mineventgap, mindur, maxdur, minPeakNum); %%slowtheta.startT(i), .peakT, .endT, .dur, .amp
        %%%%%%%%%%% here need to filter thru the peak/trough times: choose the ones truly spike (negative/positive) peaks
        %%%%%%%%%%%       then match the positive/negative peaks to the negative/positive peaks
        slowtheta = adjustminmaxtimes(slowtheta, dat, timestamp, sth);
        %disp('-----ignore adjusting min max times');
        %%%%%%%%%%%%now assign session variables        
        eeg.slowtheta.sessNum{i} = numel(slowtheta.startT); eeg.slowtheta.sessRate{i} = numel(slowtheta.startT)/eeg.general.sesslength{i};
        eeg.slowtheta.sessTimeRatio{i} = sum(slowtheta.dur)/eeg.general.sesslength{i};
        eeg.slowtheta.sessMeanAmp{i} = mean(slowtheta.amp); eeg.slowtheta.sessMeanDur{i} = mean(slowtheta.dur); eeg.slowtheta.sessMeanEng{i} = mean(slowtheta.eng);
        eeg.slowtheta.sessStartT{i} = slowtheta.startT; eeg.slowtheta.sessPeakT{i} = slowtheta.peakT; eeg.slowtheta.sessEndT{i} = slowtheta.endT;
        eeg.slowtheta.sessAmp{i} = slowtheta.amp; eeg.slowtheta.sessDur{i} = slowtheta.dur; eeg.slowtheta.sessEng{i} = slowtheta.eng;
        eeg.slowtheta.sessFreq{i} = slowtheta.freq;
        eeg.slowtheta.minTime{i} = slowtheta.minTime; eeg.slowtheta.minAmp{i} = slowtheta.minAmp;
        eeg.slowtheta.maxTime{i} = slowtheta.maxTime; eeg.slowtheta.maxAmp{i} = slowtheta.maxAmp;
        eeg.slowtheta.sessMeanFreq{i} = mean(slowtheta.freq);
        
        %%%event variables
        evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i};
        for (j = 1:numel(evTime))
            eeg.slowtheta.evtNum{i}(j) = 0; eeg.slowtheta.evtRate{i}(j) = 0; eeg.slowtheta.evtMeanAmp{i}(j) = NaN;
            eeg.slowtheta.evtMeanFreq{i}(j) = NaN; eeg.slowtheta.evtMeanDur{i}(j) = NaN; eeg.slowtheta.evtMeanEng{i}(j) = NaN;
            st = evTime{j}.start; et = evTime{j}.ent; nev = numel(st); evtid = []; leng = 0;
            for (ttt = 1:nev)
                ik = find( (slowtheta.peakT>=st(ttt)) & (slowtheta.peakT<et(ttt)) ); evtid = union(evtid, ik); ik = []; leng = leng + et(ttt)-st(ttt);
            end
            if (leng > 0)
                eeg.slowtheta.evtNum{i}(j) = numel(slowtheta.startT(evtid)); eeg.slowtheta.evtRate{i}(j) = numel(slowtheta.startT(evtid))/leng;
                eeg.slowtheta.evtMeanAmp{i}(j) = mean(slowtheta.amp(evtid)); eeg.slowtheta.evtMeanDur{i}(j) = mean(slowtheta.dur(evtid));
                eeg.slowtheta.evtTimeRatio{i}(j) = sum(slowtheta.dur(evtid))/leng;
                eeg.slowtheta.evtMeanEng{i}(j) = mean(slowtheta.eng(evtid)); eeg.slowtheta.evtMeanFreq{i}(j) = mean(slowtheta.freq(evtid));
            end
        end
    end
else %%%if not a cluster0 file, get rid of those parm assignments only for slowtheta files
     eeg.parm.slowthetaNorm{i} = []; eeg.parm.slowthetaPeakThres(i) = NaN; eeg.parm.slowthetaStartThres(i) = NaN; 
     eeg.parm.slowthetaMaxGap(i) = NaN; eeg.parm.slowthetaMinDur(i) = NaN; eeg.parm.slowthetaMaxDur(i) = NaN;
end
end

function slowtheta = findallslowthetas(timestamp, dat, pth, sth, mineventgap, mindur, maxdur, minPeakNum) %%slowtheta.startT(i), .peakT, .endT, .dur, .amp
slowtheta.startT = []; slowtheta.peakT = []; slowtheta.endT = []; slowtheta.dur = []; slowtheta.amp = []; slowtheta.eng = []; 
slowtheta.minTime = []; slowtheta.maxTime = []; slowtheta.minAmp = []; slowtheta.maxAmp = []; slowtheta.freq = [];
 
index = find(dat <= sth); %all threshold negative crossing ---Use this one because MUA agrees with max negative peaks more than the max positive peaks
if (~isempty(index))
    %[startindex, endindex] = findeventindex(index, maxgap); %%start/end indices for start threshold crossing
    [startindex, endindex] = findeventindex(index, 1); %%start/end indices for start threshold crossing
    %%%%now do a filtering to select those with durations within [mindur maxdur]
    alldur = endindex-startindex+1; iii = find( (alldur>=mindur) & (alldur<=maxdur) );
    startindex = startindex(iii); endindex = endindex(iii); %disp('****duraton filtering: '); disp(numel(startindex));
    if (~isempty(startindex))
       % disp(['****threshold filtering: ' num2str(numel(startindex)) '% ' num2str(mean(endindex-startindex))]);
       %%%%now join events if event gaps are too small 
       [startindex, ttt] = sort(startindex); endindex = endindex(ttt);
       gapnow = startindex(2:numel(startindex)) - endindex(1:numel(startindex)-1);
       smallgapind = find( gapnow < mineventgap ); largegapind = setdiff( [1:numel(gapnow)], smallgapind );
       newstart = startindex(largegapind+1); newend = endindex(largegapind);
       addstart = []; addend = [];
       if (~isempty(smallgapind))
          [sind, eind] = findeventindex(smallgapind, 1); %%start/end indices for start threshold crossing
          addstart = startindex(sind); addend = endindex(eind+1); 
       end
       startindex = union(union(startindex(1), newstart), addstart);
       endindex = union(union(endindex(numel(endindex)), newend), addend);
       startindex = sort(startindex); endindex = sort(endindex); 
       %disp(['****gap filtering: ' num2str(numel(startindex)) '% ' num2str(mean(endindex-startindex))]);
       %%%%%%%compute parameters in the events
       slowtheta = findamp(dat, timestamp, startindex, endindex, pth, minPeakNum, sth); %%filter thrlugh pth and compute
       %disp([startindex(1:100)' endindex(1:100)']);
       %disp(['****peak amp/number filtering: ' num2str(numel(slowtheta.startT))]);
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

function slowtheta = findamp(dat, timestamp, startindex, endindex, pth, minPeakNum, sth)
%%individual event indicated by [startindex, endindex], indices across startthreshold
%%peaktime, amplitude = maximum values within [startindex endindex] (minimum for slow wave detection)
%%starttime = index in [startindex-groupspan startindex] but cross the startthreshold
%%endtime = index in [endindex endindex+groupspan]but cross the startthreshold
%%area = integrate through [tarttime endtime]
slowtheta.startT = []; slowtheta.peakT = []; slowtheta.endT = []; slowtheta.dur = []; slowtheta.amp = []; slowtheta.eng = []; 
slowtheta.minTime = []; slowtheta.maxTime = []; slowtheta.minAmp = []; slowtheta.maxAmp = []; slowtheta.freq = [];
nevent = 0;
for (i = 1:numel(startindex))
     datnow = dat(startindex(i):endindex(i)); timenow = timestamp(startindex(i):endindex(i));
     [peaknow, peakindex] = min(datnow); %%if either negative crossing
     if (peaknow <= pth)
         [minindex, maxindex] = FindLocal(datnow); peakVal = datnow(minindex); %%%%!!!THis is different from ripples
         %disp(['this is ok:' num2str(numel(minindex))]);
         if (numel(find(peakVal<=pth))>= minPeakNum)
             nevent = nevent + 1;
             slowtheta.amp(nevent) = peaknow; slowtheta.peakT(nevent) = timenow(peakindex);
             slowtheta.startT(nevent) = timestamp(startindex(i)); slowtheta.endT(nevent) = timestamp(endindex(i));
             slowtheta.dur(nevent) = slowtheta.endT(nevent) - slowtheta.startT(nevent);
             slowtheta.freq(nevent) = numel(find(peakVal<=sth))/slowtheta.dur(nevent);
             slowtheta.eng(nevent) = mean(abs(datnow))* slowtheta.dur(nevent);
             slowtheta.minTime = [slowtheta.minTime timenow(minindex)]; slowtheta.maxTime = [slowtheta.maxTime timenow(maxindex)];
             slowtheta.minAmp = [slowtheta.minAmp datnow(minindex)]; slowtheta.maxAmp = [slowtheta.maxAmp datnow(maxindex)];
         end
     end
end

function slowtheta = adjustminmaxtimes(slowtheta, dat, timestamp, sth)
iip = find(slowtheta.minAmp <= sth); iiq = find(slowtheta.maxAmp > sth);
negmanV = mean(abs(slowtheta.minAmp(iip))); posmanV = mean(abs(slowtheta.maxAmp(iiq)));
%if (numel(iip) >= numel(iiq)) %%%if HV spikes are negative
if negmanV >= posmanV
    kkp = iip;
    slowtheta.minTime = slowtheta.minTime(kkp); slowtheta.minAmp = slowtheta.minAmp(kkp);
    maxAmp = NaN*ones(size(slowtheta.minAmp)); maxTime = NaN*ones(size(slowtheta.minTime));
    if (~isempty(kkp))
        for (ttp = 1:numel(kkp)-1)
            kkq = find( (slowtheta.maxTime>slowtheta.minTime(ttp)) & (slowtheta.maxTime<min([slowtheta.minTime(ttp+1) slowtheta.minTime(ttp)+0.1])) );
            if (~isempty(kkq))
                timenow = slowtheta.maxTime(kkq); ampfine = slowtheta.maxAmp(kkq);
                %[ttt, ddd] = min(timenow-slowtheta.minTime(ttp)); maxAmp(ttp) = ampfine(ddd); maxTime(ttp) = timenow(ddd); %%%%closest maxima
                [ddd, ttt] = max(ampfine); maxAmp(ttp) = ddd; maxTime(ttp) = timenow(ttt); %%%%highest point
            end
        end
        kkq = find( (slowtheta.maxTime>slowtheta.minTime(numel(kkp))) & (slowtheta.maxTime<slowtheta.minTime(numel(kkp))+0.1) );
        if (~isempty(kkq))
            timenow = slowtheta.maxTime(kkq); ampfine = slowtheta.maxAmp(kkq);
            %[ttt, ddd] = min(timenow-slowtheta.minTime(numel(kkp))); maxAmp(numel(kkp)) = ampfine(ddd); maxTime(numel(kkp)) = timenow(ddd); %%%%closest maxima
            [ddd, ttt] = max(ampfine); maxAmp(numel(kkp)) = ddd; maxTime(numel(kkp)) = timenow(ttt);
        end
    end
    slowtheta.maxAmp = maxAmp; slowtheta.maxTime = maxTime;
else %%%if HV spikes are positive
    kkp = iiq;
    slowtheta.maxTime = slowtheta.maxTime(kkp); slowtheta.maxAmp = slowtheta.maxAmp(kkp);
    minAmp = NaN*ones(size(slowtheta.maxAmp)); minTime = NaN*ones(size(slowtheta.maxTime));
    if (~isempty(kkp))
        for (ttp = 1:numel(kkp)-1)
            kkq = find( (slowtheta.minTime>slowtheta.maxTime(ttp)) & (slowtheta.minTime<min([slowtheta.maxTime(ttp+1) slowtheta.maxTime(ttp)+0.1])) );
            if (~isempty(kkq))
                timenow = slowtheta.minTime(kkq); ampfine = slowtheta.minAmp(kkq);
                %[ttt, ddd] = min(timenow-slowtheta.maxTime(ttp)); minAmp(ttp) = ampfine(ddd); minTime(ttp) = timenow(ddd); %%%%closest maxima
                [ddd, ttt] = min(ampfine); minAmp(ttp) = ddd; minTime(ttp) = timenow(ttt); %%%%highest point
            end
        end
        kkq = find( (slowtheta.minTime>slowtheta.maxTime(numel(kkp))) & (slowtheta.minTime<slowtheta.maxTime(numel(kkp))+0.1) );
        if (~isempty(kkq))
            timenow = slowtheta.minTime(kkq); ampfine = slowtheta.minAmp(kkq);
            %[ttt, ddd] = min(timenow-slowtheta.maxTime(numel(kkp))); minAmp(numel(kkp)) = ampfine(ddd); minTime(numel(kkp)) = timenow(ddd); %%%%closest maxima
            [ddd, ttt] = min(ampfine); minAmp(numel(kkp)) = ddd; minTime(numel(kkp)) = timenow(ttt);
        end
    end
    slowtheta.minAmp = minAmp; slowtheta.minTime = minTime; 
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
if (numel(mindex) < 1)
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









