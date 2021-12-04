function DataManager_ReComputeDatabase_NoStops
%%% (pinfo,data,cellind,behav,bhdata,vv)

hf = gcbf; hgroup = getappdata(hf, 'hgroup');hfield = getappdata(hf, 'hfield');
plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
fname = get(hf, 'Name'); ftt = strfind(fname, '__'); currentfilename = fname(ftt+2:numel(fname));
[~, ~, ee] = fileparts(currentfilename);
[writefilename, okk] = getoutputfile(ee, currentfilename, ow);


%do it for all cellind
if okk && strcmp(ee, '.spikedb')
    if (plotparm.linkbehav == 0)
        disp(['--------> no behavioral (non-stop) database linked']); okk = 0;
    else
        pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
        behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
        [pinfo,data] = ReComputeSpikeDBNow(pinfo,data,behav,bhdata);
    end
end
if okk && strcmp(ee, '.eegdb')
    if (plotparm.linkbehav == 0)
        disp(['--------> no behavioral data linked']); okk = 0;
    else
        disp('--------> EEG database cannont be recomputed at this time; create new event files instead.'); okk = 0;
        %pinfo = getappdata(hf, 'eeg'); data = getappdata(hf, 'eegdata');
        %behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
        %[pinfo,data] = ReComputeEEGDBNow(pinfo,data,behav,bhdata);
    end
end
if okk && strcmp(ee, '.behavdb')
    pinfo = getappdata(hf, 'behav'); data = getappdata(hf, 'bhdata');
    if cc    
        [pinfo,data] = ReComputeBehavDBNow(pinfo,data);
    else
        pinfo = assignparameter(pinfo);
    end
end
if okk 
    if (ow)
        iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
        fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
        newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
    else
        hmain = DataManager_DataManager_Callback; plotparm.linkspike = 1; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
    end
    setappdata(hmain,'plotparm', plotparm);
    set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
    if strcmp(ee, '.spikedb')
        s = whos('data'); Gb1 = s.bytes/(1024^3);
        s = whos('pinfo'); Gb2 = s.bytes/(1024^3);
        if Gb1 + Gb2 < 2
           save(writefilename, 'pinfo', 'data'); %disp('Small one');
        else
           save(writefilename, 'pinfo', 'data', '-mat', '-v7.3'); %disp('Big one');
        end
        setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data);
        DataManager_PlotSpikeDatabase(hmain, pinfo, data);
    elseif strcmp(ee, '.behavdb')
        behav = pinfo; bhdata = data; 
        save(writefilename, 'behav', 'bhdata');
        setappdata(hmain, 'behav', pinfo); setappdata(hmain, 'bhdata', data);
        DataManager_PlotSpikeDatabase(hmain, behav, bhdata, 'Sessions', 'sessID', '.behavdb');
    %elseif strcmp(ee, '.eegdb')
    %    eeg = pinfo; eegdata = data; save(writefilename, 'eeg', 'eegdata', '-mat', '-v7.3');
    %    setappdata(hmain, 'eeg', pinfo); setappdata(hmain, 'eegdata', data);
    %    DataManager_PlotSpikeDatabase(hmain, eeg, eegdata, 'EEG File', 'eegID', '.eegdb');
    end
    pinfo = []; data = [];
end
disp('**********************');

function pinfo = assignparameter(pinfo)
nspike = numel(pinfo.general.datedir); 
%%%%%assign parameters for computing place field dynamics
pp = {'[Max min] 1D speed (cm/s)'; 'Maximum 2D speed (cm/s)'; 'Max gap (s)'; 'Minimum stop duration (s)';...
    'Smooth speed?'; 'Speed smoothing sigma (s)'}; 
def = {'10 -1000'; '15'; '0.2'; '0.5'; 'Yes'; '0.2'};
III=inputdlg(pp, 'Parameters for identifying stopping points', 6, def, 'on');
if (~isempty(III))
   if (~isfield(pinfo.parm, 'max1DSpeed')) pinfo.parm.max1DSpeed = zeros(1, nspike); end
   if (~isfield(pinfo.parm, 'min1DSpeed')) pinfo.parm.min1DSpeed = zeros(1, nspike); end
   if (~isfield(pinfo.parm, 'max2DSpeed')) pinfo.parm.max2DSpeed = zeros(1, nspike); end
   if (~isfield(pinfo.parm, 'maxgap')) pinfo.parm.maxgap = zeros(1, nspike); end
   if (~isfield(pinfo.parm, 'minStopDur')) pinfo.parm.minStopDur = zeros(1, nspike); end 
   if (~isfield(pinfo.parm, 'smoothSpeed')) pinfo.parm.smoothSpeed = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'smoothSpeedSigma')) pinfo.parm.smoothSpeedSigma = zeros(1, nspike); end
   aa = str2num(III{1});
   for (i = 1:nspike) pinfo.parm.max1DSpeed(i) = aa(1); pinfo.parm.min1DSpeed(i) = aa(2); end
   for (i = 1:nspike) pinfo.parm.max2DSpeed(i) = str2num(III{2}); end
   for (i = 1:nspike) pinfo.parm.maxgap(i) = str2num(III{3}); end
   for (i = 1:nspike) pinfo.parm.minStopDur(i) = str2num(III{4}); end
   for (i = 1:nspike) pinfo.parm.smoothSpeed{i} = III{5}; end
   for (i = 1:nspike) pinfo.parm.smoothSpeedSigma(i) = str2num(III{6}); end
end

%%%%%%%%%%%%%%  *********** Filter out stopping periods on open platform and linear track *******  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [behav,bhdata] = ReComputeBehavDBNow(behav,bhdata)
%%%%%%%%%Regarding events: event (start end) times cannot be changed, but an event can be empty if all stops within
%%%%%%%%%
%%%%Identifying stopping points
%%%%Filter out stopping points; eliminate them from timestamps 
%%%%redo 2D and 1D analysis
nsess = numel(behav.general.datedir);
if (~isfield(behav.behavior, 'sessActDur')) behav.behavior.sessActDur = cell(1, nsess); end %
if (~isfield(behav.behavior, 'sessStopEvtNum')) behav.behavior.sessStopEvtNum = cell(1, nsess); end %%number of low speed laps during each run session
if (~isfield(behav.behavior, 'sessStoptotDur')) behav.behavior.sessStoptotDur = cell(1, nsess); end
if (~isfield(behav.behavior, 'sessStopTimePercent')) behav.behavior.sessStopTimePercent = cell(1, nsess); end
if (~isfield(behav.behavior, 'evtStopEvtNum')) behav.behavior.evtStopEvtNum = cell(1, nsess); end %%number of low speed laps during each run session
if (~isfield(behav.behavior, 'evtStopTotDur')) behav.behavior.evtStopTotDur = cell(1, nsess); end
if (~isfield(behav.behavior, 'evtStopEvtPerLap')) behav.behavior.evtStopEvtPerLap = cell(1, nsess); end
if (~isfield(behav.behavior, 'evtStopDurPerLap')) behav.behavior.evtStopDurPerLap = cell(1, nsess); end
if (~isfield(behav.behavior, 'evtStopTimePercent')) behav.behavior.evtStopTimePercent = cell(1, nsess); end
if (~isfield(behav.behavior, 'sessStopStart')) behav.behavior.sessStopStart = cell(1, nsess); end
if (~isfield(behav.behavior, 'sessStopEnd')) behav.behavior.sessStopEnd = cell(1, nsess); end
if (~isfield(behav.behavior, 'evtStopStart')) behav.behavior.evtStopStart = cell(1, nsess); end
if (~isfield(behav.behavior, 'evtStopEnd')) behav.behavior.evtStopEnd = cell(1, nsess); end
if (~isfield(bhdata.event, 'sessStopevTimes')) bhdata.event.sessStopevTimes = cell(1, nsess); end
if (~isfield(bhdata.event, 'evtStopevTimes')) bhdata.event.evtStopevTimes = cell(1, nsess); end

for (i = 1:nsess)
    max1Dspeed = behav.parm.max1DSpeed(i); min1Dspeed = -10000; 
    if isfield(behav.parm, 'min1DSpeed') min1Dspeed = behav.parm.min1DSpeed(i); end 
    max2Dspeed = behav.parm.max2DSpeed(i); minstopdur = behav.parm.minStopDur(i);
    smoothspeed = behav.parm.smoothSpeed{i}; speedsigma = behav.parm.smoothSpeedSigma(i); framerate = behav.parm.framerate(i);
    maxgap = behav.parm.maxgap(i);
    %if (~strcmp(behav.parm.sessType{i}, 'sleep'))
       disp(strcat('-----> filter out stopping periods and recompute position data properties ---', behav.general.sessID{i}));      
       %%%% session stopping periods: 2D speed here is the smoothed instantaneous velocity
       allvel = bhdata.sess.AllVel{i}; 
       [startind, endind] = identifystops(allvel, max2Dspeed, 0, minstopdur, smoothspeed, speedsigma, framerate, maxgap);
       postimestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i); npos = numel(postimestamp);
       stopev.start = postimestamp(startind); stopev.ent = postimestamp(endind); nstop = numel(startind);
       bhdata.event.sessStopevTimes{i} = stopev; behav.behavior.sessStopEvtNum{i} = nstop; 
       behav.behavior.sessStopStart{i} = stopev.start; behav.behavior.sessStopEnd{i} = stopev.ent;
       if (nstop>0) 
          behav.behavior.sessStoptotDur{i} = sum(stopev.ent-stopev.start);
       else
          behav.behavior.sessStoptotDur{i} = 0;
       end
        behav.behavior.sessStopTimePercent{i} = behav.behavior.sessStoptotDur{i} / behav.behavior.sessDur{i};
       %%%%%%%%%%%get rid of stop pointss from actual postimestamp and position data
       sessbadind = [];
       for (j = 1:nstop)
           sessbadind = union(sessbadind, (startind(j):endind(j)));
       end
       sessposind = setdiff( (1:npos-1), sessbadind); 
       %%%%[bhdata.pos.postimestamp{i}, ij] = sort(posnow);
       bhdata.pos.postimestamp{i} = bhdata.pos.postimestamp{i}(sessposind); 
       for (ik = 1:numel(bhdata.pos.XX{i}))
           bhdata.pos.XX{i}{ik} = bhdata.pos.XX{i}{ik}(sessposind); bhdata.pos.YY{i}{ik} = bhdata.pos.YY{i}{ik}(sessposind);
       end
       %%%%recompute 2D properties
       behav.behavior.sessActDur{i} = behav.behavior.sessDur{i} - behav.behavior.sessStoptotDur{i}; 
       postimestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i); 
       stime = behav.general.sessstartT{i}; etime = behav.general.sessendT{i}; 
       segtime = behav.parm.sessSegTime(i); nseg = floor((etime-stime)/segtime);
       Pmarker = behav.parm.sessPmarker{i}; Fmarker = behav.parm.sessFmarker{i}; Bmarker = behav.parm.sessBmarker{i}; 
       allposmarker = behav.general.posMarker{i}; 
       posXX = []; posYY = []; ik = find(strcmp(allposmarker, Pmarker)); 
       if (numel(ik) == 1) posXX = bhdata.pos.XX{i}{ik}*behav.parm.pixelXSize(i); posYY = bhdata.pos.YY{i}{ik}*behav.parm.pixelYSize(i); end
       fXX = []; fYY = []; ik = find(strcmp(allposmarker, Fmarker)); 
       if (numel(ik) == 1) fXX = bhdata.pos.XX{i}{ik}*behav.parm.pixelXSize(i); fYY = bhdata.pos.YY{i}{ik}*behav.parm.pixelYSize(i); end
       bXX = []; bYY = []; ik = find(strcmp(allposmarker, Bmarker)); 
       if (numel(ik) == 1) bXX = bhdata.pos.XX{i}{ik}*behav.parm.pixelXSize(i); bYY = bhdata.pos.YY{i}{ik}*behav.parm.pixelYSize(i); end
       %%compute session parameters
       %AllVel = ComputeV(postimestamp, posXX, posYY, 1); AllHeadDir = ComputeDir(postimestamp, fXX, bXX, fYY, bYY);
       if ~isempty(posXX)
           bhdata.sess.AllVel{i} = allvel(sessposind); bhdata.sess.AllHeadDir{i} = bhdata.sess.AllHeadDir{i}(sessposind); 
              aa = allvel(sessposind); aa = aa(~isnan(aa));
           behav.behavior.sessMeanSpeed{i} = mean(aa); behav.behavior.sessMedSpeed{i} = median(aa);
       else
           bhdata.sess.AllVel{i} = []; bhdata.sess.AllHeadDir{i} = []; behav.behavior.sessMeanSpeed{i} = NaN;
       end
       xybin = bhdata.sess.gridXYbin{i}; 
       occutime = OccupancyPoint(xybin{1}, xybin{2}, posXX, posYY,0)/framerate;
       bhdata.sess.gridOccuptime{i} = occutime;
       for (j = 1:nseg)
            stnow = stime + (j-1)*segtime; etnow = stnow + segtime; if (etnow>etime) etnow = etime; end
            tii = find( (postimestamp>=stnow)&(postimestamp<etnow) );
            if ~isempty(posXX)
                pxnow = posXX(tii); pynow = posYY(tii);
            else
                pxnow = []; pynow = [];
            end
            occutime = OccupancyPoint(xybin{1}, xybin{2}, pxnow, pynow, 0)/framerate;
            bhdata.sess.gridSegOccuptime{i}{j} = occutime; 
       end
    
       %%%%%event stopping periods
       if (strcmp(behav.parm.sessType{i}, 'linear'))
           allevstopstart = []; allevstopent = [];
           evname = behav.general.eventname{i}; evType = behav.parm.eventType{i}; evTimes = bhdata.event.eventtimes{i};
           %lowVtime = 0; lowevTimes = []; lowVLapInd = [];
       for (j = 1:numel(evname))
         %lowevTimes{j} = []; lowVLapInd{j} = [];
         behav.behavior.evtStopEvtNum{i}(j) = NaN; behav.behavior.evtStopTotDur{i}(j) = NaN; bhdata.event.evtStopevTimes{i}{j} = [];
         behav.behavior.evtStopEvtPerLap{i}(j) = NaN; behav.behavior.evtStopDurPerLap{i}(j) = NaN; behav.behavior.evtStopTimePercent{i}(j) = NaN;
         if (strcmp(evType{j}, 'run')) && contains(behav.parm.eventPosltr{i}{j}, '.ltr')  %if this is run event & can be linearized
             nlap = numel(evTimes{j}.start); 
             nstop =0; stopT = 0; stopevnow.start = []; stopevnow.ent = [];
             xbin = bhdata.event.Xbin{i}{j}; 
             for (k = 1:nlap) %%%find stopping points lap by lap
                 lappostime = bhdata.event.LapAllPostimestamp{i}{j}{k}; npos = numel(lappostime); %%%%time unit is already in second 
                 lapx = bhdata.event.LapAllX{i}{j}{k};
                 %allvel = bhdata.event.LapAllVel{i}{j}{k}; 
                 %[startind, endind] = identifystopsfromstopev(lappostime, stopev);
                 %lapspeed = diff(lapx') ./ diff(lappostime); lapspeed = [0; lapspeed]; %%%colum vector with same dimension as lappostime
                 lapspeed = bhdata.event.LapAll1DSpeed{i}{j}{k}; lapspeed = lapspeed';
                 [startind, endind] = identifystops(lapspeed, max1Dspeed, min1Dspeed, minstopdur, smoothspeed, speedsigma, framerate, maxgap);
                 if (numel(startind)>0) 
                     nstop = nstop + numel(startind); 
                     stopevnow.start = [stopevnow.start; lappostime(startind)]; stopevnow.ent = [stopevnow.ent; lappostime(endind)];
                     stopTnow = sum(lappostime(endind)-lappostime(startind));
                     stopT = stopT + stopTnow;
                     %%%%%%%%%%% get rid of stop pointss from postimestamp and position data in event variables
                     evbadind = [];
                     for (m = 1:numel(startind))
                         evbadind = union(evbadind, (startind(m):endind(m)));
                     end
                     %%% be very careful here: there could be session stopping points in event laps too
                     sessbadind = getseestoppoints(behav.behavior.sessStopStart{i}, behav.behavior.sessStopEnd{i}, lappostime);
                     allbadind = sort(union(evbadind, sessbadind));
                     posind = setdiff( (1:npos-1), allbadind); allstopTnow = numel(allbadind)/framerate;  
                     %%%%re-asign lap parameters
                     bhdata.event.LapTrajleng{i}{j}(k) = NaN; %not computed here - too cumbersome to do
                     bhdata.event.LapAllPostimestamp{i}{j}{k} = lappostime(posind);
                     bhdata.event.LapAllVel{i}{j}{k} = bhdata.event.LapAllVel{i}{j}{k}(posind);
                     bhdata.event.LapAllHeadDir{i}{j}{k}= bhdata.event.LapAllHeadDir{i}{j}{k}(posind);
                     bhdata.event.LapDur{i}{j}(k)=bhdata.event.LapDur{i}{j}(k) - allstopTnow;
                     bhdata.event.LapAllX{i}{j}{k}=lapx(posind);
                     bhdata.event.LapAllY{i}{j}{k}=bhdata.event.LapAllY{i}{j}{k}(posind);
                     bhdata.event.LapMeanSpeed{i}{j}(k)= (max(xbin)-min(xbin))/((evTimes{j}.ent(k)-evTimes{j}.start(k))-allstopTnow);
                       aa = bhdata.event.LapAll1DSpeed{i}{j}{k}(posind);
                     bhdata.event.LapAll1DSpeed{i}{j}{k} = aa; bhdata.event.LapMed1DSpeed{i}{j}(k) = median(aa(~isnan(aa)));
                     %%%%re-do occupancy time
                     bhdata.event.LapOccuptime{i}{j}{k}= (histc(lapx(posind), xbin))'/framerate;
                     %%%%% lap trajectory length not done! bhdata.event.LapTrajleng{i}{j}(k) = findtrajleng(, ynow);
                 end
             end
             behav.behavior.evtStopEvtNum{i}(j) = nstop; behav.behavior.evtStopTotDur{i}(j) = stopT; 
             behav.behavior.evtStopEvtPerLap{i}(j) = nstop/nlap; behav.behavior.evtStopDurPerLap{i}(j) = stopT/nlap;
             behav.behavior.evtStopTimePercent{i}(j) = stopT/sum(evTimes{j}.ent - evTimes{j}.start)*100;
             bhdata.event.evtStopevTimes{i}{j} = stopevnow;
             allevstopstart = [allevstopstart; stopevnow.start]; allevstopent = [allevstopent; stopevnow.ent];
             alloccu = zeros(numel(xbin),1);
             for (k = 1:nlap)
                 alloccu = alloccu + bhdata.event.LapOccuptime{i}{j}{k};
             end
             bhdata.event.Occuptime{i}{j} = alloccu; 
             %%%%%%re-set event variables
             behav.behavior.evtTotDur{i}(j) = sum(bhdata.event.LapDur{i}{j});
               aa = bhdata.event.LapDur{i}{j};  
             behav.behavior.evtMeanDur{i}(j) = mean(aa); behav.behavior.evtMedDur{i}(j) = median(aa);
               aa = addtogether(bhdata.event.LapAllVel{i}{j}); aa = aa(~isnan(aa));
             behav.behavior.evtMeanSpeed{i}(j) = mean(aa); behav.behavior.evtMedSpeed{i}(j) = median(aa);
               aa = addtogether(bhdata.event.LapAll1DSpeed{i}{j}); aa = aa(~isnan(aa));
             behav.behavior.evtMean1DSpeed{i}(j) = mean(aa); behav.behavior.evtMed1DSpeed{i}(j) = median(aa);
             %%%%% below is not computed - too cumbersome to do
             behav.behavior.evtMeanTrajleng{i}(j) = NaN; behav.behavior.evtMedTrajleng{i}(j) = NaN;             
         end
       end  
             [allevstopstart, iii] = sort(allevstopstart); allevstopent = allevstopent(iii);
             behav.behavior.evtStopStart{i} = allevstopstart; behav.behavior.evtStopEnd{i} = allevstopent;
       end
    %end
end
function badind = getseestoppoints(stopstart, stopend, lappostime)
badind = [];
for i =1:numel(stopstart)
    iii = find( (lappostime>=stopstart(i)) & (lappostime<=stopend(i)) );
    badind = union(badind, iii);
end

function [pinfo,data] = ReComputeSpikeDBNow(pinfo,data,behav,bhdata)
%%%%%%Regarding events: (start end) times cannot be changed!! But may contain no spikes after spike elimination
%%%%%%change data.spike.spikteime: get rid of spikes during stopping periods (both 2D session stopping events and 1D event stopping events)
%%%%%%recompute firing properties; keep wave properties intact
nspike = numel(pinfo.general.datedir);
if (~isfield(pinfo.parm, 'max1DSpeed')) pinfo.parm.max1DSpeed = zeros(1, nspike); end
if (~isfield(pinfo.parm, 'max2DSpeed')) pinfo.parm.max2DSpeed = zeros(1, nspike); end
if (~isfield(pinfo.parm, 'maxGap')) pinfo.parm.maxgap = zeros(1, nspike); end
if (~isfield(pinfo.parm, 'minStopDur')) pinfo.parm.minStopDur = zeros(1, nspike); end 
if (~isfield(pinfo.parm, 'smoothSpeed')) pinfo.parm.smoothSpeed = cell(1, nspike); end
if (~isfield(pinfo.parm, 'smoothSpeedSigma')) pinfo.parm.smoothSpeedSigma = zeros(1, nspike); end
for (i = 1:nspike)
    disp(strcat('-----> filter out stopping periods and recompute spike properties ---', pinfo.general.clname{i}));
    finaldirnow = pinfo.general.finaldir{i};
    spiketime = data.spike.spiketime{i}*pinfo.parm.timeunit(i);
    binsize = pinfo.parm.maxratetimebin(i); BIint = pinfo.parm.spikeburstint(i); RIint = pinfo.parm.spikerefracint(i);
    sessions = pinfo.general.sessionname{i}; nsess = numel(sessions);
    for (tt = 1:nsess)
        posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, sessions{tt}) );
        if numel(posid) == 1
           if isfield(behav.parm, 'max1DSpeed') pinfo.parm.max1DSpeed(i) = behav.parm.max1DSpeed(posid); end
           if isfield(behav.parm, 'max2DSpeed') pinfo.parm.max2DSpeed(i) = behav.parm.max2DSpeed(posid); end
           if isfield(behav.parm, 'minStopDur') pinfo.parm.minStopDur(i) = behav.parm.minStopDur(posid); end
           if isfield(behav.parm, 'smoothSpeed') pinfo.parm.smoothSpeed{i} = behav.parm.smoothSpeed{posid}; end
           if isfield(behav.parm, 'smoothSpeedSigma') pinfo.parm.smoothSpeedSigma(i) = behav.parm.smoothSpeedSigma(posid); end 
           if isfield(behav.parm, 'maxgap') pinfo.parm.maxgap(i) = behav.parm.maxgap(posid); end
           break
        end
    end
    sessstopev = cell(size(sessions));
    for (tt = 1:nsess)
        %disp(pinfo.parm.sessType{i}{tt});
      if (~strcmp(pinfo.parm.sessType{i}{tt}, 'sleep')) %%%do not eliminate spikes in sleep sessions
        posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, sessions{tt}) );
        %disp(['----> found posid: ', num2str(posid)]);
        sessstopev{tt}.ent = []; sessstopev{tt}.start = []; 
        if numel(posid) == 1
            %disp('---> inside now');
           if isfield(bhdata.event, 'sessStopevTimes')
              sessstopev{tt} = bhdata.event.sessStopevTimes{posid};
           elseif isfield(bhdata.event, 'sessPatchevTimes')
              sessstopev{tt} = bhdata.event.sessPatchevTimes{posid};
           end
           spiketime = eliminatespikes(spiketime, sessstopev{tt}.start, sessstopev{tt}.ent); %%%eliminate spikes in 2D session stopping periods
           sT = pinfo.general.sessionstartT{i}(tt); eT = pinfo.general.sessionendT{i}(tt); 
           [sT, eT] = filterouteventtimesnow(sT, eT, sessstopev{tt});
           [spikeN, meanrate, maxrate, BI, RI] = ComputeRateNumBI(spiketime, sT, eT, binsize, BIint, RIint);
           pinfo.firing.sessspikeN{i}(tt) = spikeN; pinfo.firing.sessmeanrate{i}(tt) = meanrate; pinfo.firing.sessmaxrate{i}(tt) = maxrate;
           pinfo.firing.sessBI{i}(tt) = BI; pinfo.firing.sessRI{i}(tt) = RI;
        end
      end
    end
    evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i}; evName = pinfo.general.eventname{i}; 
    for (j = 1:numel(evName))
        if (strcmp(evType{j}, 'run')) 
            %%%%locate event position data
            [evSess, sessid] = identifysession(evTime{j}, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
            posid = []; evid = []; nev = numel(evTime{j}.start);
            if (~isempty(evSess))
                posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess) );
            end
            if numel(posid) == 1
                if (isfield(behav.general, 'eventname'))
                    evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
                else
                    evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
                end
            end
            if (numel(posid)~=1)||(numel(evid)~=1)
               disp(['-------------> stopping periods not filtered: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
            elseif contains(behav.parm.eventPosltr{posid}{evid}, '.ltr')  %if this is run event & can be linearized
               stopev = bhdata.event.evtStopevTimes{posid}{evid}; 
               spiketime = eliminatespikes(spiketime, stopev.start, stopev.ent);%%%%eliminate spikes in event stopping periods
               sT = evTime{j}.start; eT = evTime{j}.ent; 
               [sT, eT] = filterouteventtimesnow(sT, eT, stopev); [sT, eT] = filterouteventtimesnow(sT, eT, sessstopev{sessid});
               [spikeN, meanrate, maxrate, BI, RI] = ComputeRateNumBI(spiketime, sT, eT, binsize, BIint, RIint);
               %%%%%% The above needs to be careful: spiketime already
               %%%%%% eliminated those in both session and event stopev's
               %%%%%% timepoints need to eliminate both too. 
               pinfo.firing.evtmaxrate{i}(j) = maxrate; pinfo.firing.evtmeanrate{i}(j) = meanrate; 
               pinfo.firing.evtBI{i}(j) = BI; pinfo.firing.evtRI{i}(j) = RI; pinfo.firing.evtspikeN{i}(j) = spikeN;
            end
        end
    end
    data.spike.spiketime{i} = spiketime / pinfo.parm.timeunit(i);
end

function aa = addtogether(A)
aa = []; 
for (i=1:numel(A))
    S = size(A{i}); if S(1)>1 A{i} = A{i}'; end
    aa = [aa A{i}];
end
  
function [evSess, sessid] = identifysession(evTime, sessionname, startT, endT)
evSess = []; sessid = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii}; sessid = iii;
end

function spiketime = eliminatespikes(spiketime, st, et)
nspike = numel(spiketime); ind = [];
for (i = 1:numel(st))
    ii = find( (spiketime>=st(i)) & (spiketime<=et(i)) );
    ind = union(ind, ii);
end
jjj = setdiff( (1:nspike), ind );
spiketime = spiketime(jjj);

function [writefilename, okk] = getoutputfile(ext, currentfilename, ow)
okk = 1; writefilename = [];
if (ow == 0)
   [fname, pname] = uiputfile(fullfile(cd, ext), 'Write the new database to:');
   if (numel(fname)>1)
      writefilename = fullfile(pname, fname);
   else
      okk = 0;
   end
else
   writefilename = currentfilename;      
end

function [startind, endind] = identifystops(oldvel, maxspeed, minspeed, minstopdur, smoothspeed, sigma, framerate, maxgap)
startind = []; endind = []; vel = oldvel(2:numel(oldvel)); %%%get rid of first point
%%% also get rid of rediculus speed point: >= 50*std
aaa = abs(vel(~isnan(vel))); %%%%vel already smoothed and all >= 0
if ~isempty(aaa) 
    ss = std(aaa); md = median(aaa);
    ind = find(abs(vel)>md+50*ss); vel(ind) = md*ones(size(ind)); 
end
if (~isempty(vel))
if strncmpi(smoothspeed, 'yes', 1)
   sigma = round(sigma*framerate); %%%sigma now if number of frames
   XX = [-5*sigma:1:5*sigma]; wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
   vel = conv(vel, wind); vel = vel( ((nw-1)/2+1):(numel(vel)-(nw-1)/2) );
end
maxgap = round(maxgap*framerate);  mindur =  round(minstopdur*framerate);
index = find( (vel<=maxspeed) & (vel>=minspeed) ); %all threshold crossing: vel cound be negative, but only for 1D speed; for 2D, vel always >=0
if (~isempty(index))
    [startindnow, endindnow] = findeventindex(index, maxgap); %%start/end indices for start threshold crossing 
    iii = find( (endindnow-startindnow)>= mindur );
    startind = startindnow(iii); endind = endindnow(iii);
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
       nindex = numel(spanindex); startindex = ones(1,nindex+1); endindex = ones(1,nindex+1);
       startindex(1) = index(1);
       endindex(1) = index(spanindex(1));
       for (i = 2:nindex)
           startindex(i) = index(spanindex(i-1)+1);
           endindex(i) = index(spanindex(i));
       end
       startindex(nindex+1) = index(spanindex(nindex)+1);
       endindex(nindex+1) = index(numel(index));
   end
end

% iii = find(vel<maxspeed); %%%%%%%all low speed time points
% if (~isempty(iii))
%    jjj = diff(iii); gapind = find(jjj>1); ngap = numel(gapind);
%    gapstartind = zeros(1, ngap+1); gapendind = zeros(1, ngap+1); 
%    gapstartind(1)= iii(1); 
%    if (ngap > 0)
%       gapendind(1) = iii(gapind(1))+1;
%       for (i = 2:ngap)
%        gapstartind(i) = iii(gapind(i-1)+1); gapendind(i) = iii(gapind(i))+1;
%       end
%       gapstartind(ngap+1) = iii(gapind(ngap)+1); gapendind(ngap+1) = iii(numel(iii))+1;
%    else
%       gapendind(1) = iii(numel(iii))+1;
%    end
%    gapleng = (gapendind - gapstartind)/framerate;
%    jjj = find(gapleng >= minstopdur); 
%    startind = gapstartind(jjj); endind = gapendind(jjj);
% end
% end

function [spikeN, meanrate, maxrate, BI, RI] = ComputeRateNumBI(spiketime, starttime, endtime, binsize, BIint, RIint)
spikeN = NaN; meanrate = NaN; maxrate = NaN; BI = NaN; RI = NaN;
if (~isempty(starttime))
    ep.start = starttime; ep.ent = endtime;
    [spikenow, epid] = SpikeEventFilter(spiketime, ep);
    spikeN = numel(spikenow); meanrate = spikeN/sum(endtime-starttime);
    %%%maxrate & BI
    maxrr = zeros(numel(starttime),1); burstspike = 0; totalspike = 0; refracspike = 0;
    for (i = 1:numel(starttime))
        index = find( (spikenow>=starttime(i)) & (spikenow<=endtime(i)) ); spikeok = sort(spikenow(index)); timebin = [];
        if (~isempty(index))
           nbin = ceil( (endtime(i)-starttime(i))/binsize );
           timebin = starttime(i) + ([1:nbin] - 1) * binsize;
           count = histc(spikeok', timebin); maxrr(i) = max(count)/binsize;
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

function [ev1now, ev2now] = filterouteventtimesnow(ev1now, ev2now, evTimes)%%%%This is too stringent, need to find overlaps
timeunit = 0.005; %time resolution
%%%%% Finding timepoints in [ev1now ev2now], but not in evTimes.start/ent (is an infinitely complex problem)
%%% Here use a comprehensive counting approach: judge point by point
evTstart = evTimes.start; evTent = evTimes.ent;
%%%%now find the overlap between ev1/2now and evTstart/ent
nev = numel(ev1now); nT = numel(evTstart); 
%disp(['filter start end time ' num2str(evTstart(1)) '* end time ' num2str(evTent(numel(evTent)))]);
if (~isempty(ev1now)) && (~isempty(evTstart))
   st = ev1now(1); et = ev2now(nev); allT = st:timeunit:et; 
   s1=findpoints(allT, ev1now, ev2now); s2=findpoints(allT, evTstart, evTent);
   sel = s1&(~s2); [sind, eind] = parseindex(sel);
   ev1now = allT(sind); ev2now = allT(eind); 
end
function sel = findpoints(allT, sT, eT)
sel = zeros(1,numel(allT)); iii = []; 
for (j = 1:numel(sT))
     inow = find( ((allT>=sT(j)) & (allT<=eT(j))) );
     iii = union(iii, inow);
end
sel(iii) = ones(1,numel(iii));
function [sind, eind] = parseindex(sel)
sind = []; eind = [];
if numel(sel)>1
   if isempty(find(sel==0))
       sind = 1; eind = numel(sel);
   else
       D = diff(sel); sind = (find(D==1))+1; eind = find(D==-1);
       if sel(1)==1 sind = [1 sind]; end
       if sel(numel(sel))==1 eind = [eind numel(sel)]; end
   end
end
if isempty(sind) || isempty(eind)
    disp('-----------> warning: event starts or ends not found; no events generated after deleting stopping times')
    sind = []; eind = [];    
elseif (numel(sind) ~= numel(eind)) || (~isempty(find(sind-eind>0)))
    disp('-----------> warning: event starts and ends disordered; no events generated after deleting stopping times')
    sind = []; eind = [];
end
