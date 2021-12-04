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
if (~isfield(eeg.spindle, 'startThreshold')) eeg.spindle.startThreshold = cell(1, neeg); end
if (~isfield(eeg.spindle, 'peakThreshold')) eeg.spindle.peakThreshold = cell(1, neeg); end

if (~isfield(eeg.spindle, 'sessNum')) eeg.spindle.sessNum = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessRate')) eeg.spindle.sessRate = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanAmp')) eeg.spindle.sessMeanAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanFreq')) eeg.spindle.sessMeanFreq = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanDur')) eeg.spindle.sessMeanDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessMeanEng')) eeg.spindle.sessMeanEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'evtNum')) eeg.spindle.evtNum = cell(1, neeg); end
if (~isfield(eeg.spindle, 'evtRate')) eeg.spindle.evtRate = cell(1, neeg); end
if (~isfield(eeg.spindle, 'evtMeanAmp')) eeg.spindle.evtMeanAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'evtMeanFreq')) eeg.spindle.evtMeanFreq = cell(1, neeg); end
if (~isfield(eeg.spindle, 'evtMeanDur')) eeg.spindle.evtMeanDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'evtMeanEng')) eeg.spindle.evtMeanEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'sessStartT')) eeg.spindle.sessStartT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessEndT')) eeg.spindle.sessEndT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessPeakT')) eeg.spindle.sessPeakT = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessDur')) eeg.spindle.sessDur = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessAmp')) eeg.spindle.sessAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'sessEng')) eeg.spindle.sessEng = cell(1, neeg); end

if (~isfield(eeg.spindle, 'minTime')) eeg.spindle.minTime = cell(1, neeg); end
if (~isfield(eeg.spindle, 'maxTime')) eeg.spindle.maxTime = cell(1, neeg); end
if (~isfield(eeg.spindle, 'minAmp')) eeg.spindle.minAmp = cell(1, neeg); end
if (~isfield(eeg.spindle, 'maxAmp')) eeg.spindle.maxAmp = cell(1, neeg); end

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
    for (j = 1:numel(evTime))
        eeg.spindle.evtNum{i}(j) = 0; eeg.spindle.evtRate{i}(j) = 0; eeg.spindle.evtMeanAmp{i}(j) = NaN; 
        eeg.spindle.evtMeanFreq{i}(j) = NaN; eeg.spindle.evtMeanDur{i}(j) = NaN; eeg.spindle.evtMeanEng{i}(j) = NaN;
        st = evTime{j}.start; et = evTime{j}.ent; nev = numel(st); evtid = []; leng = 0;
        for (ttt = 1:nev)
            ik = find( (spin.peakT>=st(ttt)) & (spin.peakT<et(ttt)) ); evtid = union(evtid, ik); ik = []; leng = leng + et(ttt)-st(ttt);
        end
        if (leng > 0)
            eeg.spindle.evtNum{i}(j) = numel(spin.startT(evtid)); eeg.spindle.evtRate{i}(j) = numel(spin.startT(evtid))/leng;
            eeg.spindle.evtMeanAmp{i}(j) = mean(spin.amp(evtid)); eeg.spindle.evtMeanDur{i}(j) = mean(spin.dur(evtid));
            eeg.spindle.evtMeanEng{i}(j) = mean(spin.eng(evtid));
        end
    end    
    end
else %%%if not a spindle file, get rid of those parm assignments only for spindle files
    eeg.parm.spinNorm{i} = []; eeg.parm.spinPeakThres(i) = NaN; eeg.parm.spinStartThres(i) = NaN; 
    eeg.parm.spinMaxGap(i) = NaN; eeg.parm.spinMinDur(i) = NaN; eeg.parm.spinMaxDur(i) = NaN;
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




