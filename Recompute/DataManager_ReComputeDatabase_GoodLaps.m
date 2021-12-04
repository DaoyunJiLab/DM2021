function DataManager_ReComputeDatabase_GoodLaps
%%% Eliminate bad laps in linearized positon data
%%% It is not recommended to eliminated too many "bad" laps, as long as stopping periods are eliminated

hf = gcbf; hgroup = getappdata(hf, 'hgroup');hfield = getappdata(hf, 'hfield');
plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
fname = get(hf, 'Name'); ftt = strfind(fname, '__'); currentfilename = fname(ftt+2:numel(fname));
[pp, nn, ee] = fileparts(currentfilename);
[writefilename, okk] = getoutputfile(ee, currentfilename, ow);
%do it for all cellind
if okk && strcmp(ee, '.spikedb')
    if (plotparm.linkbehav == 0)
        disp(['--------> no behavioral data linked']); okk = 0;
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
        disp('--------> EEG database cannot be re-computed at this time; create new event file instead.'); okk = 0;
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
        behav = pinfo; bhdata = data; save(writefilename, 'behav', 'bhdata');
        setappdata(hmain, 'behav', pinfo); setappdata(hmain, 'bhdata', data);
        DataManager_PlotSpikeDatabase(hmain, behav, bhdata, 'Sessions', 'sessID', '.behavdb');
    %elseif strcmp(ee, '.eegdb')
    %    eeg = pinfo; eegdata = data; save(writefilename, 'eeg', 'eegdata');
    %    setappdata(hmain, 'eeg', pinfo); setappdata(hmain, 'eegdata', data);
    %    DataManager_PlotSpikeDatabase(hmain, eeg, eegdata, 'EEG File', 'eegID', '.eegdb');
    end
    pinfo = []; data = [];
end
disp('**********************');

function pinfo = assignparameter(pinfo)
nspike = numel(pinfo.general.datedir); 
%%%%%assign parameters for computing place field dynamics
pp = {'Min median speed (cm/s):'; 'Max lap duration (s)'}; 
def = {'5'; '50'}; %%%default should be leneint '10', '40'
III=inputdlg(pp, 'Parameters for identifying good run laps', 1, def, 'on');
if (~isempty(III))
   if (~isfield(pinfo.parm, 'minLaMedSpeed')) pinfo.parm.minLapMedSpeed = zeros(1, nspike); end
   if (~isfield(pinfo.parm, 'maxLapDur')) pinfo.parm.maxLapDur = zeros(1, nspike); end
   for (i = 1:nspike) pinfo.parm.minLapMedSpeed(i) = str2num(III{1}); pinfo.parm.maxLapDur(i) = str2num(III{2}); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [behav,bhdata] = ReComputeBehavDBNow(behav,bhdata)
%%%%Do nothing to session variables; 
%%%%Filter out low-speed laps; eliminate them from computed event variables
%%%%(actual raw positon data not affected)
nsess = numel(behav.general.datedir);
if (~isfield(behav.behavior, 'evtLowVLapNum')) behav.behavior.evtLowVLapNum = cell(1, nsess); end %%number of low speed laps during each run session
if (~isfield(behav.behavior, 'evtLowVLapPercent')) behav.behavior.evtLowVLapPercent = cell(1, nsess); end %%Percentage of low speed laps among all laps
if (~isfield(behav.behavior, 'evtLowVLapTotDur')) behav.behavior.evtLowVLapTotDur = cell(1, nsess); end
if (~isfield(bhdata.event, 'lowVLapTimes')) bhdata.event.lowVLapTimes = cell(1, nsess); end
if (~isfield(bhdata.event, 'lowVLapInd')) bhdata.event.lowVLapInd = cell(1, nsess); end

for (i = 1:nsess)
if (strcmp(behav.parm.sessType{i}, 'linear'))
    disp(strcat('-----> filter out lowV laps and recompute position data properties ---', behav.general.sessID{i}));
    MinLapMedSpeed = behav.parm.minLapMedSpeed(i); MaxLapDur = behav.parm.maxLapDur(i);
    postimestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i); npos = numel(postimestamp);
    evname = behav.general.eventname{i}; evType = behav.parm.eventType{i}; evTimes = bhdata.event.eventtimes{i};
    lowevTimes = []; lowVLapInd = [];%badind = [];
    for (j = 1:numel(evname))
         lowevTimes{j} = []; lowVLapInd{j} = [];
         if (strcmp(evType{j}, 'run')) && contains(behav.parm.eventPosltr{i}{j}, '.ltr')  %if this is run event & can be linearized
             lapspeednow = bhdata.event.LapMed1DSpeed{i}{j};
             lapdurnow = bhdata.event.eventtimes{i}{j}.ent - bhdata.event.eventtimes{i}{j}.start; %%bhdata.event.LapDur{i}{j}; this shrinks after stopping removal; use original instead
             if ~isempty(find((size(lapspeednow) ~= size(lapdurnow)))) lapdurnow = lapdurnow'; end %%%if their dimensions not match (row vs column), flip it
             iii = find( (lapspeednow>=MinLapMedSpeed) & (lapdurnow<=MaxLapDur) ); %%%%good laps
             jjj = find( (lapspeednow<MinLapMedSpeed) | (lapdurnow>MaxLapDur) ); 
             lowVLapInd{j} = jjj; %%%%bad laps
             %%%%re-asign lap parameters
             bhdata.event.eventtimes{i}{j}.start = bhdata.event.eventtimes{i}{j}.start(iii);
             bhdata.event.eventtimes{i}{j}.ent = bhdata.event.eventtimes{i}{j}.ent(iii);
             bhdata.event.LapAllPostimestamp{i}{j} = bhdata.event.LapAllPostimestamp{i}{j}(iii);
             bhdata.event.LapAllVel{i}{j} = bhdata.event.LapAllVel{i}{j}(iii);
             bhdata.event.LapAllHeadDir{i}{j}= bhdata.event.LapAllHeadDir{i}{j}(iii);
             bhdata.event.LapDur{i}{j}=bhdata.event.LapDur{i}{j}(iii);
             bhdata.event.LapTrajleng{i}{j}=bhdata.event.LapTrajleng{i}{j}(iii);
             bhdata.event.LapOccuptime{i}{j}=bhdata.event.LapOccuptime{i}{j}(iii);
             bhdata.event.LapAllX{i}{j}=bhdata.event.LapAllX{i}{j}(iii);
             bhdata.event.LapAllY{i}{j}=bhdata.event.LapAllY{i}{j}(iii);
             bhdata.event.LapMeanSpeed{i}{j}=bhdata.event.LapMeanSpeed{i}{j}(iii);
             bhdata.event.LapAll1DSpeed{i}{j} = bhdata.event.LapAll1DSpeed{i}{j}(iii); 
             bhdata.event.LapMed1DSpeed{i}{j} = bhdata.event.LapMed1DSpeed{i}{j}(iii);
             behav.behavior.evtNum{i}(j) = numel(iii); 
             
             %%%%%%re-set event variables
             behav.behavior.evtTotDur{i}(j) = sum(bhdata.event.LapDur{i}{j});
               aa = bhdata.event.LapDur{i}{j};  
             behav.behavior.evtMeanDur{i}(j) = mean(aa); behav.behavior.evtMedDur{i}(j) = median(aa);
               aa = addtogether(bhdata.event.LapAllVel{i}{j}); aa = aa(~isnan(aa));
             behav.behavior.evtMeanSpeed{i}(j) = mean(aa); behav.behavior.evtMedSpeed{i}(j) = median(aa);
               aa = addtogether(bhdata.event.LapAll1DSpeed{i}{j}); aa = aa(~isnan(aa));
             behav.behavior.evtMean1DSpeed{i}(j) = mean(aa); behav.behavior.evtMed1DSpeed{i}(j) = median(aa);
               aa = bhdata.event.LapTrajleng{i}{j}; aa = aa(~isnan(aa));
             behav.behavior.evtMeanTrajleng{i}(j) = NaN; behav.behavior.evtMedTrajleng{i}(j) = NaN;
             if ~isempty(aa) behav.behavior.evtMeanTrajleng{i}(j) = mean(aa); behav.behavior.evtMedTrajleng{i}(j) = median(aa); end
            
             behav.behavior.evtLowVLapNum{i}(j) = numel(jjj); behav.behavior.evtLowVLapPercent{i}(j) = 100*numel(jjj)/(numel(iii)+numel(jjj));
             if ~isempty(jjj)
                behav.behavior.evtLowVLapTotDur{i}(j)= sum(evTimes{j}.ent(jjj) - evTimes{j}.start(jjj));
             else
                behav.behavior.evtLowVLapTotDur{i}(j) = 0;
             end
             %%%re-asign total run occup time
             xbin = bhdata.event.Xbin{i}{j};
             alloccu = zeros(numel(xbin),1);
             for (tt = 1:numel(bhdata.event.eventtimes{i}{j}.start))
                  alloccu = alloccu + bhdata.event.LapOccuptime{i}{j}{tt};
             end
             bhdata.event.Occuptime{i}{j} = alloccu;
             lowevTimes{j}.start = evTimes{j}.start(jjj); lowevTimes{j}.ent = evTimes{j}.ent(jjj); 
% %              %%%%timestamps to be eliminated - why not do this? - See below
%              for (tt = 1:numel(jjj))
%                  indnow = find( (postimestamp>=evTimes{j}.start(jjj(tt))) & (postimestamp<=evTimes{j}.ent(jjj(tt))) );
%                  badind = union(badind, indnow); 
%              end
         end
    end
    %%%reset event-related session variables
    bhdata.event.lowVLapInd{i} = lowVLapInd;
    bhdata.event.lowVLapTimes{i} = lowevTimes;
    
%     %%%eliminated bad lap timestamps and position points and redo 2D
%     %%%properties -----Do NOT DO this!-----TO KEEP CONSISTENT WITH
%     %%%No-Stop processing: eliminate 1D bad events only within linearized
%     %%%position data
%     bhdata.event.lowVLapInd{i} = lowVLapInd;
%     posind = setdiff( (1:npos), badind); bhdata.pos.postimestamp{i} = bhdata.pos.postimestamp{i}(posind); 
%     %%%%[bhdata.pos.postimestamp{i}, ij] = sort(posnow);
%     for (ik = 1:numel(bhdata.pos.XX{i}))
%         bhdata.pos.XX{i}{ik} = bhdata.pos.XX{i}{ik}(posind); bhdata.pos.YY{i}{ik} = bhdata.pos.YY{i}{ik}(posind);
%     end
%     %%%%recompute 2D properties
%     postimestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i); framerate = behav.parm.framerate(i); 
%     stime = behav.general.sessstartT{i}; etime = behav.general.sessendT{i}; 
%     segtime = behav.parm.sessSegTime(i); nseg = floor((etime-stime)/segtime);
%     Pmarker = behav.parm.sessPmarker{i}; Fmarker = behav.parm.sessFmarker{i}; Bmarker = behav.parm.sessBmarker{i}; 
%     allposmarker = behav.general.posMarker{i}; 
%     posXX = []; posYY = []; ik = find(strcmp(allposmarker, Pmarker)); 
%     if (numel(ik) == 1) posXX = bhdata.pos.XX{i}{ik}; posYY = bhdata.pos.YY{i}{ik}; end
%     fXX = []; fYY = []; ik = find(strcmp(allposmarker, Fmarker)); 
%     if (numel(ik) == 1) fXX = bhdata.pos.XX{i}{ik}; fYY = bhdata.pos.YY{i}{ik}; end
%     bXX = []; bYY = []; ik = find(strcmp(allposmarker, Bmarker)); 
%     if (numel(ik) == 1) bXX = bhdata.pos.XX{i}{ik}; bYY = bhdata.pos.YY{i}{ik}; end
%     %%compute session parameters
%     AllVel = ComputeV(postimestamp, posXX, posYY, 1); AllHeadDir = ComputeDir(postimestamp, fXX, bXX, fYY, bYY);
%     bhdata.sess.AllVel{i} = AllVel; bhdata.sess.AllHeadDir{i} = AllHeadDir; 
%     xybin = bhdata.sess.gridXYbin{i}; 
%     occutime = OccupancyPoint(xybin{1}, xybin{2}, posXX, posYY,0)/framerate;
%     bhdata.sess.gridOccuptime{i} = occutime;
%     for (j = 1:nseg)
%             stnow = stime + (j-1)*segtime; etnow = stnow + segtime; if (etnow>etime) etnow = etime; end
%             tii = find( (postimestamp>=stnow)&(postimestamp<etnow) );
%             occutime = OccupancyPoint(xybin{1}, xybin{2}, posXX(tii), posYY(tii), 0)/framerate;
%             bhdata.sess.gridSegOccuptime{i}{j} = occutime; 
%     end
    
end
end

function [pinfo,data] = ReComputeSpikeDBNow(pinfo,data,behav,bhdata)
%%%%%%update data.events.eventtimes
%%%%%%change data.spike.spikteime: get rid of spikes during event bad laps
%%%%%%recompute firing properties
%%%%%keep wave properties intact
nspike = numel(pinfo.general.datedir);
if (~isfield(pinfo.parm, 'minLapSpeed')) pinfo.parm.minLapSpeed = zeros(1, nspike); end
for (i = 1:nspike)
    disp(strcat('-----> filter out lowV laps and recompute spike properties ---', pinfo.general.clname{i}));
    finaldirnow = pinfo.general.finaldir{i};
    spiketime = data.spike.spiketime{i}*pinfo.parm.timeunit(i);
    binsize = pinfo.parm.maxratetimebin(i); BIint = pinfo.parm.spikeburstint(i); RIint = pinfo.parm.spikerefracint(i);
    evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i}; evName = pinfo.general.eventname{i}; 
    for (j = 1:numel(evName))
        if (strcmp(evType{j}, 'run')) 
            %%%%locate event position data
            [evSess, sessID] = identifysession(evTime{j}, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
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
               disp(['-------------> low V laps not filtered: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
            elseif contains(behav.parm.eventPosltr{posid}{evid}, '.ltr')  %if this is run event & can be linearized
               if isfield(bhdata.event, 'sessStopevTimes')
                  sessstopev = bhdata.event.sessStopevTimes{posid};
               elseif isfield(bhdata.event, 'sessPatchevTimes')
                  sessstopev = bhdata.event.sessPatchevTimes{posid};
               end
               stopev = bhdata.event.evtStopevTimes{posid}{evid}; 
               pinfo.parm.minLapMedSpeed(i) = behav.parm.minLapMedSpeed(posid); 
               pinfo.parm.maxLapDur(i) = behav.parm.maxLapDur(posid);
               jjj = bhdata.event.lowVLapInd{posid}{evid}; iii = setdiff( (1:nev), jjj );
               data.events.eventtimes{i}{j}.start = evTime{j}.start(iii); data.events.eventtimes{i}{j}.ent = evTime{j}.ent(iii); 
               spiketime = eliminatespikes(spiketime, evTime{j}.start(jjj), evTime{j}.ent(jjj));
               sT = evTime{j}.start(iii); eT = evTime{j}.ent(iii); 
               [sT, eT] = filterouteventtimesnow(sT, eT, stopev); [sT, eT] = filterouteventtimesnow(sT, eT, sessstopev);
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
%[pinfo,data] = DataManager_FindSpikeFiring(pinfo,data, 1:nspike, 0);%%%%This is wrong when run events were chopped (e.g. excluding stopping times)

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

function aa = addtogether(A)
aa = []; 
for (i=1:numel(A))
    S = size(A{i}); if S(1)>1 A{i} = A{i}'; end
    aa = [aa A{i}];
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
