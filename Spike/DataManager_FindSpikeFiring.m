function [pinfo,data] = DataManager_FindSpikeFiring(pinfo,data, cellind, vv)
%%find spike firing information by looking at spike data
%%fields assigned here:
if (~isempty(cellind))
  nspike = numel(pinfo.general.animalname);
  if (~isfield(pinfo, 'firing')) pinfo.firing = []; end
  if (~isfield(pinfo.firing, 'sessspikeN')) pinfo.firing.sessspikeN = cell(1, nspike); end %number of spikes for all task sessions
  if (~isfield(pinfo.firing, 'sessmaxrate')) pinfo.firing.sessmaxrate = cell(1, nspike); end %maximum firing rate during task (averaged across all task sessions), bined 1s
  if (~isfield(pinfo.firing, 'sessmeanrate')) pinfo.firing.sessmeanrate = cell(1, nspike); end %overall firing rate during task
  if (~isfield(pinfo.firing, 'sessBI')) pinfo.firing.sessBI = cell(1, nspike); end %bursting index during task
  if (~isfield(pinfo.firing, 'sessRI')) pinfo.firing.sessRI = cell(1, nspike); end %refraction index during task
  
%   if (~isfield(pinfo.firing, 'taskspikeN')) pinfo.firing.taskspikeN = cell(1, nspike); end %number of spikes for all task sessions
%   if (~isfield(pinfo.firing, 'taskmaxrate')) pinfo.firing.taskmaxrate = cell(1, nspike); end %maximum firing rate during task (averaged across all task sessions), bined 1s
%   if (~isfield(pinfo.firing, 'taskmeanrate')) pinfo.firing.taskmeanrate = cell(1, nspike); end %overall firing rate during task
%   if (~isfield(pinfo.firing, 'taskBI')) pinfo.firing.taskBI = cell(1, nspike); end %bursting index during task
%   if (~isfield(pinfo.firing, 'taskRI')) pinfo.firing.taskRI = cell(1, nspike); end %refraction index during task
%   
%   if (~isfield(pinfo.firing, 'sleepspikeN')) pinfo.firing.sleepspikeN = cell(1, nspike); end %number of spikes for all sleep sessions
%   if (~isfield(pinfo.firing, 'sleepmaxrate')) pinfo.firing.sleepmaxrate = cell(1, nspike); end %maximum firing rate during sleep (averaged across all sleep sessions), bined 1s
%   if (~isfield(pinfo.firing, 'sleepmeanrate')) pinfo.firing.sleepmeanrate = cell(1, nspike); end %overall firing rate during sleep
%   if (~isfield(pinfo.firing, 'sleepBI')) pinfo.firing.sleepBI = cell(1, nspike); end %bursting index during sleep
%   if (~isfield(pinfo.firing, 'sleepRI')) pinfo.firing.sleepRI = cell(1, nspike); end %refraction index during sleep
  
  if (~isfield(pinfo.firing, 'evtspikeN')) pinfo.firing.evtspikeN = cell(1, nspike); end %number of spikes for all run events
  if (~isfield(pinfo.firing, 'evtmaxrate')) pinfo.firing.evtmaxrate = cell(1, nspike); end %maximum firing rate during run (averaged across all run sessions), bined 1s
  if (~isfield(pinfo.firing, 'evtmeanrate')) pinfo.firing.evtmeanrate = cell(1, nspike); end %overall firing rate during run
  if (~isfield(pinfo.firing, 'evtBI')) pinfo.firing.evtBI = cell(1, nspike); end %bursting index during run
  if (~isfield(pinfo.firing, 'evtRI')) pinfo.firing.evtRI = cell(1, nspike); end %refraction index during run
  
%   if (~isfield(pinfo.firing, 'stopspikeN')) pinfo.firing.stopspikeN = cell(1, nspike); end %number of spikes for all run events
%   if (~isfield(pinfo.firing, 'stopmaxrate')) pinfo.firing.stopmaxrate = cell(1, nspike); end %maximum firing rate during run (averaged across all run sessions), bined 1s
%   if (~isfield(pinfo.firing, 'stopmeanrate')) pinfo.firing.stopmeanrate = cell(1, nspike); end %overall firing rate during run
%   if (~isfield(pinfo.firing, 'stopBI')) pinfo.firing.stopBI = cell(1, nspike); end %bursting index during run
%   if (~isfield(pinfo.firing, 'stopRI')) pinfo.firing.stopRI = cell(1, nspike); end %refraction index during run
%   
%   if (~isfield(pinfo.firing, 'SWSspikeN')) pinfo.firing.SWSspikeN = cell(1, nspike); end %number of spikes for all SWS events
%   if (~isfield(pinfo.firing, 'SWSmaxrate')) pinfo.firing.SWSmaxrate = cell(1, nspike); end %maximum firing rate during SWS (averaged across all SWS sessions), bined 1s
%   if (~isfield(pinfo.firing, 'SWSmeanrate')) pinfo.firing.SWSmeanrate = cell(1, nspike); end %overall firing rate during SWS
%   if (~isfield(pinfo.firing, 'SWSBI')) pinfo.firing.SWSBI = cell(1, nspike); end %bursting index during SWS 
%   if (~isfield(pinfo.firing, 'SWSRI')) pinfo.firing.SWSRI = cell(1, nspike); end %refraction index during SWS
%   
%   if (~isfield(pinfo.firing, 'REMspikeN')) pinfo.firing.REMspikeN = cell(1, nspike); end %number of spikes for all REM events
%   if (~isfield(pinfo.firing, 'REMmaxrate')) pinfo.firing.REMmaxrate = cell(1, nspike); end %maximum firing rate during REM sleep, bined 1s
%   if (~isfield(pinfo.firing, 'REMmeanrate')) pinfo.firing.REMmeanrate = cell(1, nspike); end %overall firing rate during REM
%   if (~isfield(pinfo.firing, 'REMBI')) pinfo.firing.REMBI = cell(1, nspike); end %bursting index during REM sleep
%   if (~isfield(pinfo.firing, 'REMRI')) pinfo.firing.REMRI = cell(1, nspike); end %refraction index during REM
end  
for (jjjjk = 1:numel(cellind))
    i = cellind(jjjjk);
    disp(strcat('-----> compute firing rates and burst indices -- ', pinfo.general.parmfile{i}));
    timeunit = pinfo.parm.timeunit(i); spiketime = data.spike.spiketime{i}*timeunit;
    binsize = pinfo.parm.maxratetimebin(i); BIint = pinfo.parm.spikeburstint(i); RIint = pinfo.parm.spikerefracint(i);
    
    sessions = pinfo.general.sessionname{i}; nsess = numel(sessions);
    for (tt = 1:nsess)
        sT = pinfo.general.sessionstartT{i}(tt); eT = pinfo.general.sessionendT{i}(tt); 
        [spikeN, meanrate, maxrate, BI, RI] = ComputeRateNumBI(spiketime, sT, eT, binsize, BIint, RIint);
        pinfo.firing.sessspikeN{i}(tt) = spikeN; pinfo.firing.sessmeanrate{i}(tt) = meanrate; pinfo.firing.sessmaxrate{i}(tt) = maxrate;
        pinfo.firing.sessBI{i}(tt) = BI; pinfo.firing.sessRI{i}(tt) = RI;
    end
    
    %%%%firing in events 
    evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i};
    for (j = 1:numel(evTime))
        starttime = evTime{j}.start; endtime = evTime{j}.ent;
        [spikeN, meanrate, maxrate, BI, RI] = ComputeRateNumBI(spiketime, starttime, endtime, binsize, BIint, RIint);
        pinfo.firing.evtmaxrate{i}(j) = maxrate; pinfo.firing.evtmeanrate{i}(j) = meanrate; 
        pinfo.firing.evtBI{i}(j) = BI; pinfo.firing.evtRI{i}(j) = RI; pinfo.firing.evtspikeN{i}(j) = spikeN;
    end
end

function [spikeN, meanrate, maxrate, BI, RI] = ComputeRateNumBI(spiketime, starttime, endtime, binsize, BIint, RIint)
spikeN = NaN; meanrate = NaN; maxrate = NaN; BI = NaN; RI = NaN;
if (~isempty(starttime))
    ep.start = starttime; ep.ent = endtime;
    [spikenow, epid] = SpikeEventFilter(spiketime, ep);
    spikeN = numel(spikenow); meanrate = spikeN/sum(endtime-starttime);
  if spikeN > 0
    %%%maxrate & BI
    maxrr = zeros(numel(starttime),1); burstspike = 0; totalspike = 0; refracspike = 0;
    for (i = 1:numel(starttime))
        index = find( (spikenow>=starttime(i)) & (spikenow<=endtime(i)) ); spikeok = sort(spikenow(index)); timebin = [];
        if (~isempty(index))
           nbin = ceil( (endtime(i)-starttime(i))/binsize );
           timebin = starttime(i) + ([1:nbin] - 1) * binsize;
           count = histc(spikeok', timebin);
           if ~isempty(count) maxrr(i) = max(count)/binsize; end
        end
        if (numel(index)>=1)
            spikeint = diff(spikeok); 
            burstspike = burstspike + numel(find(spikeint<=BIint)); totalspike = totalspike + numel(index);
            refracspike = refracspike + numel(find(spikeint<=RIint));
        end
    end
    maxrate = max(maxrr); 
    if (totalspike > 0)
        BI = burstspike/totalspike; RI = refracspike/totalspike; 
    end
  end
end
