function [pinfo, data] = DataManager_FindParm(pinfo, data)
%%set up all parameters for later on calculations
%%parm.timeunit = timeunit for spike timestamp (1)
%%parm.timebinsize = time bin size in second for calculating rate, etc
%%parm.space1dbin = spatial bin size for 1d linearized track
%%parm.space2dbin
%%parm.maxgap = maximum gap in pixel between fields used to determine field boundary
%%parm.zerothreshold = percentage of maximum rate in field used to treat low rate firing as zero random firing

nspike = numel(pinfo.general.animalname);
pinfo.parm.timeunit = ones(1, nspike);   %timeunit = 1s
pinfo.parm.fs = 32556 * ones(1, nspike);     %sampling frequency per channel
pinfo.parm.spikebasepoint = 0 * ones(1, nspike); %number of points to compuate the baseline of a spike waveform 

pinfo.parm.maxratetimebin = ones(1, nspike);   %time binsize for computing max rate = 1s
pinfo.parm.spikeburstint = 0.01*ones(1,nspike); %spike interval for computing spike burst index (<=10ms considered a burst) 
pinfo.parm.spikerefracint = 0.002*ones(1,nspike); %spike interval for computing spike refraction index (<=2ms considered a spike within refraction periods) 

pinfo.parm.sessType = cell(1, nspike);
pinfo.parm.eventtype = cell(1,nspike);
%%specify parameters for particular analysis (samples below):

%%select which diode position to compute position-related information
% for (i = 1:nspike) %default diode position used to extract position-related analysis: 'front', 'back', or 'middle'
%     pinfo.parm.diodepos{i} = 'green';
% end
 

% pinfo.parm.space1dbin = 100 * ones(1, nspike);   %1d grid number: linearized track spatial binsize = 100 bins for the one path
% pinfo.parm.spacefullbin = 200 * ones(1, nspike);   %1d grid number: linearized track spatial binsize = 400 bins for the whole track
% pinfo.parm.space2dbin = 32 * ones(1, nspike);   %2d track spatial binsize = 32x32 gridss for the whole track
% pinfo.parm.maxgap = 5 * ones(1, nspike); %if (zero) gap >= 5 space1dbin (grid), split to two fields
% pinfo.parm.basepercent = 0.3 * ones(1, nspike); %baseline firng rate = average of 30% grids with lowest firing rate
% pinfo.parm.zerothreshold = 0.05 * ones(1, nspike); %if firing rate < 0.05*maxrate, treat it as zero
% pinfo.parm.minfieldsize = 3 * ones(1, nspike); %if a linear field size <= 3 space1dbin, discard it
% pinfo.parm.minpeakrate = 1 * ones(1, nspike); %if a linear field with peak firing rate <= 1 hz, discard it
 
% %%good run criteria: since each animal/day varies slightly, modify for some days
% pinfo.parm.epmaxlength = 10 * ones(1, nspike); %if a run episode (ent-start) longer than 10 second, reject it as bad run
% pinfo.parm.fullepmaxlength = 50 * ones(1, nspike); %if a fullrun episode (ent-start) longer than 60 second, reject it as bad run

for (i = 1:nspike)
    sessnow = pinfo.general.sessionname{i}; 
    for (j = 1:numel(sessnow))
        pinfo.parm.sessType{i}{j} = findsesstype(sessnow{j});
    end
    evn = pinfo.general.eventname{i};
    for (j = 1:numel(evn))
        pinfo.parm.eventtype{i}{j} = findeventtype(evn{j});
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

function evType = findeventtype(evn)
evType = 'others'; 
if ( (~isempty(strfind(evn, 'cw'))) | (~isempty(strfind(evn, 'acw'))) | (~isempty(strfind(evn, 'CW'))) | (~isempty(strfind(evn, 'ACW')))...
        | (~isempty(strfind(evn, 'leftright'))) | (~isempty(strfind(evn, 'rightleft'))) | (~isempty(strfind(evn, 'Leftright')))...
        | (~isempty(strfind(evn, 'Rightleft'))) | (~isempty(strfind(evn, 'ncw'))) | (~isempty(strfind(evn, 'NCW')))...
        | (~isempty(strfind(evn, 'RL'))) | (~isempty(strfind(evn, 'LR'))) | (~isempty(strfind(lower(evn), 'run')))...
        | (~isempty(strfind(evn, 'lefttoright'))) | (~isempty(strfind(evn, 'righttoleft'))) ...
        | (~isempty(strfind(evn, 'CenterRight'))) | (~isempty(strfind(evn, 'RightCenter'))) ...
        | (~isempty(strfind(evn, 'CenterLeftt'))) | (~isempty(strfind(evn, 'LeftCenter'))) )...
    & ( (isempty(strfind(evn, 'ripple'))) & (isempty(strfind(evn, 'spindle'))) & (isempty(strfind(evn, 'full'))) ...
        & (isempty(strfind(evn, 'session'))) & (isempty(strfind(evn, 'allep'))) & (isempty(strfind(evn, 'Wrong')))  )
% if ( (~isempty(strfind(evn, 'cw'))) | (~isempty(strfind(evn, 'acw'))) | (~isempty(strfind(evn, 'CW'))) | (~isempty(strfind(evn, 'ACW')))...
%         | (~isempty(strfind(evn, 'leftright'))) | (~isempty(strfind(evn, 'rightleft'))) | (~isempty(strfind(evn, 'Leftright')))...
%         | (~isempty(strfind(evn, 'Rightleft'))) | (~isempty(strfind(evn, 'ncw'))) | (~isempty(strfind(evn, 'NCW')))...
%         | (~isempty(strfind(evn, 'RL'))) | (~isempty(strfind(evn, 'LR'))) | (~isempty(strfind(lower(evn), 'run')))...
%         | (~isempty(strfind(evn, 'lefttoright'))) | (~isempty(strfind(evn, 'righttoleft'))) )...
%     & ( (isempty(strfind(evn, 'ripple'))) & (isempty(strfind(evn, 'spindle'))) & (isempty(strfind(evn, 'full'))) ...
%         & (isempty(strfind(evn, 'session'))) & (isempty(strfind(evn, 'allep'))) )
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
if contains(lower(evn), 'pbe')
    evType = 'ripple';
end
if contains(lower(evn), 'box') && contains(lower(evn), 'sesspeakt')
    evType = 'ripple';
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