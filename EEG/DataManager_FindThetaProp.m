function [eeg, eegdata] = DataManager_FindThetaProp(eeg, eegdata, eegind, vv)
%%%detect ripple events and compute their parameters
neeg = numel(eeg.general.eegfile);
%%%the following variables are assigned
% if (~isfield(eegdata, 'ripple')) eegdata.ripple = []; end
% if (~isfield(eegdata.ripple, 'allripples')) eegdata.ripple.allripples = cell(1, neeg); end %%all computed results are here
%%%%displayed variables
if (~isfield(eeg, 'theta')) eeg.theta = []; end
if (~isfield(eeg.theta, 'sessNum')) eeg.theta.sessNum = cell(1, neeg); end
if (~isfield(eeg.theta, 'sessMeanAmp')) eeg.theta.sessMeanAmp = cell(1, neeg); end
if (~isfield(eeg.theta, 'sessMeanFreq')) eeg.theta.sessMeanFreq = cell(1, neeg); end

if (~isfield(eeg.theta, 'evtNum')) eeg.theta.evtNum = cell(1, neeg); end
if (~isfield(eeg.theta, 'evtMeanAmp')) eeg.theta.evtMeanAmp = cell(1, neeg); end
if (~isfield(eeg.theta, 'evtMeanFreq')) eeg.theta.evtMeanFreq = cell(1, neeg); end

if (~isfield(eeg.theta, 'maxTime')) eeg.theta.maxTime = cell(1, neeg); end
if (~isfield(eeg.theta, 'minTime')) eeg.theta.minTime = cell(1, neeg); end
if (~isfield(eeg.theta, 'maxAmp')) eeg.theta.maxAmp = cell(1, neeg); end
if (~isfield(eeg.theta, 'minAmp')) eeg.theta.minAmp = cell(1, neeg); end

for (iiik = 1:numel(eegind))
i = eegind(iiik);
if ~isempty(strfind(eeg.general.eegfile{i}, 'theta_smooth')) %%only do this for filtered and smoothed theta EEG traces
    disp(['--------> theta analysis: ', eeg.general.eegfile{i}]);
    %%%get eeg data
    [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
    dat = dat-mean(dat);
    %%%detect maximum and minimum timepoints
    theta = findthetapoints(timestamp, dat); %theta.minT, .maxT, .minA, .maxA, .int %%%all these have same lengths (matched)
    eeg.theta.maxTime{i} = theta.maxT; eeg.theta.minTime{i} = theta.minT; eeg.theta.maxAmp{i} = theta.maxA; eeg.theta.minAmp{i} = theta.minA; 
    
    %%%session variables
    eeg.theta.sessNum{i} = numel(theta.maxT);
    eeg.theta.sessMeanAmp{i} = mean(abs([theta.minA; theta.maxA]));
    intnow = abs(diff(sort([theta.minT theta.maxT])));
    goodint = intnow( (intnow>= 0.05) & (intnow<=0.5) ); 
    if (~isempty(goodint)) eeg.theta.sessMeanFreq{i} = 1/mean(goodint); end
    
    %%%event variables
    evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i};
    for (j = 1:numel(evTime))
        eeg.theta.evtNum{i}(j) = 0; eeg.theta.evtMeanAmp{i}(j) = NaN; eeg.theta.evtMeanFreq{i}(j) = NaN;
        st = evTime{j}.start; et = evTime{j}.ent; nev = numel(st); evtid = [];
        for (ttt = 1:nev)
            ik = find( (theta.maxT>=st(ttt)) & (theta.maxT<et(ttt)) ); evtid = union(evtid, ik); ik = []; 
        end
        if (~isempty(evtid))
            eeg.theta.evtNum{i}(j) = numel(evtid); 
            eeg.theta.evtMeanAmp{i}(j) = mean(abs([theta.minA(evtid); theta.maxA(evtid)]));
            intnow = abs(diff(sort([theta.minT(evtid) theta.maxT(evtid)]))); goodint = intnow( (intnow>= 0.05) & (intnow<=0.5) ); 
            if (~isempty(goodint)) eeg.theta.evtMeanFreq{i}(j) = 1/mean(goodint); end
        end
    end    
end
end

function theta = findthetapoints(timestamp, dat)
theta.minT = []; theta.maxT = []; theta.minA = []; theta.maxA = []; %theta.int = []; 
[minindex, maxindex] = FindLocal(dat); %%%find minima and maxima indices in dat and timestamp
%[minindex, maxindex] = findmatch(timestamp, minindex, maxindex); %%%match each minindex with each maxindex
theta.maxT = timestamp(maxindex); theta.minT = timestamp(minindex); %%theta peak times in second: reversed
theta.maxA = dat(maxindex); theta.minA = dat(minindex); %theta.int = abs(theta.maxT-theta.minT); 

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

%function [minindex, maxindex] = findmatch(timestamp, mminindex, mmaxindex)
% numel(mminindex)
% numel(mmaxindex)
% disp([mminindex(1:10)' mmaxindex(1:10)']);
% mintime = timestamp(mminindex); maxtime = timestamp(mmaxindex); size(mintime)
% alltime = [mintime; maxtime]; allind = [mminindex; mmaxindex]; size(maxtime)
% minsign = -1*ones(size(mintime)); maxsign = ones(size(maxtime)); allsign = [minsign; maxsign];
% [timeok, ind] = sort(alltime); signok = allsign(ind); allind = allind(ind);%%%now singular times are sorted. if neighboring signs are the same, eliminat the middle ones
% signdiff = diff(signok); smallind = find(signdiff == 2); bigind = find(signdiff == -2);
% finalnum = min([numel(smallind) numel(bigind)]);
% minindex = allind(smallind(1:finalnum)); maxindex = allind(bigind(1:finalnum));
% numel(minindex)
% numel(maxindex)
% disp([minindex(1:10) maxindex(1:10)]); 


