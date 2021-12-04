function [eeg, eegdata] = DataManager_FindEEGParm(eeg, eegdata)
%%set up all parameters for later on calculations
%%parm.timeunit = timeunit for eeg timestamp (1)

neeg = numel(eeg.general.eegfile);
eeg.parm.timeunit = 0.0001*ones(1, neeg);   %timeunit = 0.0001s - don't change
eeg.parm.buffersize = 512*ones(1, neeg); %Default buffersize = 512 datapoints - don't change
eeg.parm.band = cell(1, neeg); eeg.parm.sessType = cell(1,neeg); eeg.parm.eventtype = cell(1, neeg);
eeg.parm.specMinFreq = NaN*ones(1, neeg); eeg.parm.specMaxFreq = NaN*ones(1, neeg);
eeg.parm.specWinSize = NaN*ones(1, neeg); eeg.parm.specWinShift = NaN*ones(1,neeg); 

for (i = 1:neeg)
    efile = eeg.general.eegfile{i}; sessnow = eeg.general.sessname{i};
    eeg.parm.sessType{i} = findsesstype(sessnow);
    [btype, pnow] = findbandtype(efile);
    eeg.parm.band{i} = btype; eeg.parm.specMinFreq(i) = pnow.specMinFreq; eeg.parm.specMaxFreq(i) = pnow.specMaxFreq;
    eeg.parm.specFreqStep(i) = pnow.specFreqStep;
    eeg.parm.specWinSize(i) = pnow.specWinSize; eeg.parm.specWinShift(i) = pnow.specWinShift; 
    evn = eeg.general.eventname{i};
    for (j = 1:numel(evn))
        eeg.parm.eventtype{i}{j} = findeventtype(evn{j});
    end
end

function sessType = findsesstype(sessnow)
sessType = 'others';
if (~isempty(strfind(sessnow, 'open'))) | (~isempty(strfind(sessnow, 'field'))) | (~isempty(strfind(sessnow, 'OF'))) ...
        | (~isempty(strfind(sessnow, 'OP'))) | (~isempty(strfind(sessnow, 'platform'))) | (~isempty(strfind(sessnow, 'circ')))...
        | (~isempty(strfind(sessnow, 'square'))) | (~isempty(strfind(sessnow, 'Open')))...
        | (~isempty(strfind(sessnow, 'QS'))) | (~isempty(strfind(sessnow, 'QW'))) | (~isempty(strfind(sessnow, 'QN')))...
        | (~isempty(strfind(sessnow, 'QE'))) | (~isempty(strfind(sessnow, 'QX'))) | (~isempty(strfind(lower(sessnow), 'box')))
    sessType = 'open';
elseif (~isempty(strfind(sessnow, 'track'))) | (~isempty(strfind(sessnow, 'linear'))) | (~isempty(strfind(sessnow, 'LT')))...
        | (~isempty(strfind(sessnow, 'Linear'))) | (~isempty(strfind(sessnow, 'Track'))) | (~isempty(strfind(sessnow, 'TRK')))...
        | (~isempty(strfind(lower(sessnow), 'run'))) | (~isempty(strfind(lower(sessnow), 'obs'))) 
    sessType = 'linear';
elseif (~isempty(strfind(sessnow, 'sleep'))) | (~isempty(strfind(sessnow, 'Sleep'))) | (~isempty(strfind(sessnow, 'rest')))...
        | (~isempty(strfind(sessnow, 'Rest'))) | (~isempty(strfind(sessnow, 'SLP')))
    sessType = 'sleep';
end

function [btype, pnow] = findbandtype(efile)
btype = []; pnow.specMinFreq = NaN; pnow.specMaxFreq = NaN; pnow.specWinSize = NaN; pnow.specWinShift = NaN;
if (~isempty(strfind(efile, 'ripple')))
    btype = 'ripple';
elseif (~isempty(strfind(efile, 'delta')))
    btype = 'delta';
elseif (~isempty(strfind(efile, 'theta')))
    btype = 'theta';
elseif (~isempty(strfind(efile, 'sws')))
    btype = 'sws';
elseif (~isempty(strfind(efile, 'SWS')))
    btype = 'sws';
elseif (~isempty(strfind(efile, 'spindle')))
    btype = 'spindle';
elseif (~isempty(strfind(efile, 'hvs')))
    btype = 'hvs';
elseif (~isempty(strfind(efile, 'gamma')))
    btype = 'gamma';
elseif (~isempty(strfind(efile, 'beta')))
    btype = 'beta';
elseif (~isempty(strfind(efile, 'EMG')))
    btype = 'emg';
elseif (~isempty(strfind(efile, 'cluster0')))
    btype = 'cluster0';
else
    btype = 'broad'; 
end
pnow.specMinFreq = 0.5; pnow.specMaxFreq = 400; pnow.specFreqStep = 0.5; pnow.specWinSize = 4; pnow.specWinShift = 2;

function evType = findeventtype(evn)
evType = 'others'; 
if ( (~isempty(strfind(evn, 'cw'))) | (~isempty(strfind(evn, 'acw'))) | (~isempty(strfind(evn, 'CW'))) | (~isempty(strfind(evn, 'ACW')))...
        | (~isempty(strfind(evn, 'leftright'))) | (~isempty(strfind(evn, 'rightleft'))) | (~isempty(strfind(evn, 'Leftright')))...
        | (~isempty(strfind(evn, 'Rightleft'))) | (~isempty(strfind(evn, 'ncw'))) | (~isempty(strfind(evn, 'NCW')))...
        | (~isempty(strfind(evn, 'RL'))) | (~isempty(strfind(evn, 'LR'))) | (~isempty(strfind(lower(evn), 'run')))...
        | (~isempty(strfind(evn, 'lefttoright'))) | (~isempty(strfind(evn, 'righttoleft'))) )...
    & ( (isempty(strfind(evn, 'ripple'))) & (isempty(strfind(evn, 'spindle'))) & (isempty(strfind(evn, 'full'))) ...
        & (isempty(strfind(evn, 'session'))) & (isempty(strfind(evn, 'allep'))) )
   evType = 'run';
elseif (~isempty(strfind(evn, 'stop'))) | (~isempty(strfind(evn, 'Stop'))) | (~isempty(strfind(evn, 'eat'))) | (~isempty(strfind(evn, 'Eat')))...
        | (~isempty(strfind(evn, 'leftleft'))) | (~isempty(strfind(evn, 'rightright'))) | (~isempty(strfind(evn, 'Leftleft')))...
        | (~isempty(strfind(evn, 'Rightright')))
    evType = 'stop';
% elseif (~isempty(strfind(evn, 'sleep'))) | (~isempty(strfind(evn, 'Sleep')))...
%         | (~isempty(strfind(evn, 'rest'))) | (~isempty(strfind(evn, 'Rest')))
%     evType = 'sleep';
  
elseif (~isempty(strfind(evn, 'REM'))) | (~isempty(strfind(evn, 'rem')))
    evType = 'rem';
elseif (~isempty(strfind(evn, 'SWS'))) | (~isempty(strfind(evn, 'sws')))
    evType = 'sws';    
end
if (~isempty(strfind(evn, 'ripp')))
    evType = 'ripple';
end
if (~isempty(strfind(evn, 'theta')))
    evType = 'theta';
end
if (~isempty(strfind(lower(evn), 'stop'))) | (~isempty(strfind(evn, 'leftleft'))) | (~isempty(strfind(evn, 'rightright'))) %%%stop keyword override run keyword
    evType = 'stop';
end


