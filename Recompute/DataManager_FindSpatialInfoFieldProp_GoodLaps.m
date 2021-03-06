function DataManager_FindSpatialInfoFieldProp_GoodLaps
%%% (pinfo,data,cellind,behav,bhdata,vv)
hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); tagmark = get(gcbo, 'Tag');
hgroup = getappdata(hf, 'hgroup');hfield = getappdata(hf, 'hfield');
ok = 1; plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
[writefilename, okk] = getoutputfile(hf, ow);

%get selected cellind
groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
if (~isempty(cellind)) && okk
    if cc
       if (plotparm.linkbehav == 0);
          disp(['--------> no behavioral data linked']);
       else
          behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
          [pinfo,data] = DoItNow(pinfo,data,behav,bhdata,  cellind, vv);
       end
    else
       pinfo = assignparameter(pinfo, cellind);
    end
    %%%%%%%%save and plot the new database
    if (~isempty(pinfo.general.parmfile))
        if (ok)
            save(writefilename, 'pinfo', 'data');
            if (ow)
               iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
               fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
               newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
            else
               hmain = DataManager_DataManager_Callback; plotparm.linkspike = 1; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
            end
            setappdata(hmain,'plotparm', plotparm);
            DataManager_PlotSpikeDatabase(hmain, pinfo, data);
            set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
            setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); pinfo = []; data = [];
        end
    else
        disp('-------------> no cells in the database!');
    end
else
    disp('--------------> no groups selected, groups do not contain any cells');
end
disp('**********************');

function pinfo = assignparameter(pinfo, cellind)
nspike = numel(pinfo.general.parmfile); ncell = numel(cellind);
%%%%%assign parameters for computing place field dynamics
pp = {'Minimum mean speed (pixels/s)'; 'Minimum number of laps'}; 
def = {'50'; '4'};
III=inputdlg(pp, 'Parameters for identifying good run laps', 2, def, 'on');
if (~isempty(III))
   if (~isfield(pinfo.parm, 'minLapSpeed')) pinfo.parm.minLapSpeed = zeros(1, nspike); end
   if (~isfield(pinfo.parm, 'minLapNum')) pinfo.parm.minLapNum = zeros(1, nspike); end
   for (i = 1:ncell) pinfo.parm.minLapSpeed(cellind(i)) = str2num(III{1}); end
   for (i = 1:ncell) pinfo.parm.minLapNum(cellind(i)) = str2num(III{2}); end
end

%%%%%%%%Changes from the standard PF calculation:
%%%%%%%%%%%%%% 1. Filter out laps with obvious stopping points on linear track
%%%%%%%%%%%%%%  *********** 2. Filter out stopping periods on open platform and linear track ******* Not done yet, need to recompute behavdb 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Compute the spatial info and place field properties:
%%%     1. For every session, compute rate maps, then compute spatial info, and 2D place field parameters
%%%     2. For every run events, compute 1D space plots, then compute spatial info, and 1d field parameters

%require the following parameters (not all listed)
%   behav.parm.s1Dbinsize, behav.parm.s2Dbinsize, behav.parm.sessPmarker,
%   pinfo.parm.fBaseRate, pinfo.parm.fMinGap, pinfo.parm.MinLapSpeed, pinfo.parm.MinFieldSpeed
%   pinfo.parm.sessType, pinfo.parm.eventtype

%require the following variables (not all listed)
%   pinfo.general.sessionname, pinfo.general.eventname, 

%requre the following data (not all listed)
%   data.spike.spiketime; data.event.eventtimes
%   bhdata.pos.postimestamp; bhdata.pos.XX/YY
%   bhdata.sess.gridXYbin; bhdata.sess.gridOccuptime; 
%   bhdata.event.Xbin, bhdata.event.Xjoint, bhdata.event.Occuptime, 
%   bhdata.event.LapAllPostimestamp, bhdata.event.AllX

function [pinfo,data] = DoItNow(pinfo,data,behav,bhdata, cellind, vv)
%variable to assign
if (~isempty(cellind)) && (~isempty(behav))
   nspike = numel(pinfo.general.parmfile);
   %if (~isfield(pinfo.parm, 'MinLapMeanSpeed')) pinfo.parm.MinLapMeanSpeed = MinLapMeanSpeed*ones(1, nspike); end%%%minimum mean speed of the lap
   %if (~isfield(pinfo.parm, 'MinLapNum')) pinfo.parm.MinLapNum = MinLapNum*ones(1, nspike); end %%%minimum number of laps (if good laps smaller than, do not analyze)
   %pinfo.parm.MinInstantSpeed = 20; %%%minimum instantaneous speed --- not used yet
   if (~isfield(pinfo, 'GRfield')) pinfo.GRfield = []; end
   %%%evt field properties
   if (~isfield(pinfo.GRfield, 'runMeanRate')) pinfo.GRfield.runMeanRate = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'runBurstI')) pinfo.GRfield.runBurstI = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'runSptlInfo')) pinfo.GRfield.runSptlInfo = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'runSptlInfotime')) pinfo.GRfield.runSptlInfotime = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'runSparsty')) pinfo.GRfield.runSparsty = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'runActvArea')) pinfo.GRfield.runActvArea = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'run1DEvtPairs')) pinfo.GRfield.run1DEvtPairs = cell(1, nspike); end 
   if (~isfield(pinfo.GRfield, 'run1DStablty')) pinfo.GRfield.run1DStablty = cell(1, nspike); end 
   if (~isfield(pinfo.GRfield, 'run1DStabP')) pinfo.GRfield.run1DStabP = cell(1, nspike); end 
   if (~isfield(pinfo.GRfield, 'run1DTrajPairs')) pinfo.GRfield.run1DTrajPairs = cell(1, nspike); end  %directionality 
   if (~isfield(pinfo.GRfield, 'run1DDirCrr')) pinfo.GRfield.run1DDirCrr = cell(1, nspike); end  %directionality: simple correlation 
   if (~isfield(pinfo.GRfield, 'run1DDirOverlap')) pinfo.GRfield.run1DDirOverlap = cell(1, nspike); end  %directionality: overlap
   if (~isfield(pinfo.GRfield, 'run1DDirLocShift')) pinfo.GRfield.run1DDirLocShift = cell(1, nspike); end  %directionality: location shift 
   if (~isfield(pinfo.GRfield, 'run1DDirRateDiffI')) pinfo.GRfield.run1DDirRateDiffI = cell(1, nspike); end  %rate difference index for each field 
   
   if (~isfield(pinfo.GRfield, 'PF1DNfield')) pinfo.GRfield.PF1DNfield = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1Devt')) pinfo.GRfield.PF1Devt = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DBoundStart')) pinfo.GRfield.PF1DBoundStart = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DBoundEnd')) pinfo.GRfield.PF1DBoundEnd = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DLocPeakX')) pinfo.GRfield.PF1DLocPeakX = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DLocComX')) pinfo.GRfield.PF1DLocComX = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DLocStartX')) pinfo.GRfield.PF1DLocStartX = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DLocEndX')) pinfo.GRfield.PF1DLocEndX = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DLength')) pinfo.GRfield.PF1DLength = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DSize')) pinfo.GRfield.PF1DSize = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DSkewness')) pinfo.GRfield.PF1DSkewness = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DKurtosis')) pinfo.GRfield.PF1DKurtosis = cell(1, nspike); end
   if (~isfield(pinfo.GRfield, 'PF1DInMeanrate')) pinfo.GRfield.PF1DInMeanrate = cell(1, nspike); end 
   if (~isfield(pinfo.GRfield, 'PF1DInPeakrate')) pinfo.GRfield.PF1DInPeakrate = cell(1, nspike); end 
   if (~isfield(pinfo.GRfield, 'PF1DBaseRate')) pinfo.GRfield.PF1DBaseRate = cell(1, nspike); end 
      
   %%%data variables to assign
   if (~isfield(data, 'GRfield')) data.GRfield = []; end
   if (~isfield(data.GRfield, 'evt1DRateMaps')) data.GRfield.evt1DRateMaps = cell(1, nspike); end
end

for (jjjk = 1:numel(cellind))
    i = cellind(jjjk); MinLapMeanSpeed = pinfo.parm.minLapSpeed(i); MinLapNum = pinfo.parm.minLapNum(i);
    threshold1Drate = pinfo.parm.f1DThresholdRate(i); max1Dgap = pinfo.parm.f1DMaxGap(i); min1Dpeakrate = pinfo.parm.f1DMinPeakRate(i); 
    timeunit = pinfo.parm.timeunit(i); base1DrateN = pinfo.parm.f1DBaseRateGridN(i); 
    threshold2Drate = pinfo.parm.f2DThresholdRate(i); max2Dgap = pinfo.parm.f2DMaxGap(i); min2Dpeakrate = pinfo.parm.f2DMinPeakRate(i); 
    base2DrateN = pinfo.parm.f2DBaseRateGridN(i); 
    smParm.d2sigma = pinfo.parm.fSmooth2DSigma(i); smParm.d2Nsig = pinfo.parm.fSmooth2DNSigma(i); smParm.d2sm = pinfo.parm.fSmooth2D{i};
    smParm.d1sigma = pinfo.parm.fSmooth1DSigma(i); smParm.d1Nsig = pinfo.parm.fSmooth1DNSigma(i); smParm.d1sm = pinfo.parm.fSmooth1D{i};
    smParm.BIint = pinfo.parm.spikeburstint(i);
    RC(1) = pinfo.parm.fRotate2DCenX(i); RC(2) = pinfo.parm.fRotate2DCenY(i); theta = pinfo.parm.fRotate2DAngle(i); 
    dirNpix = pinfo.parm.fDirShift(i); rotornot = pinfo.parm.fRotate{i};
    disp(strcat('-----> compute field properties ---', pinfo.general.parmfile{i}));
        
    %%%%%%%%compute event (linear track) 1D rate maps
    evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i}; 
    finaldirnow = pinfo.general.finaldir{i};
    D1ratenow = []; evSess = []; xbin = []; occutime = []; runMeanRate = []; runBurstI = [];
    for (j = 1:numel(evTime))
        runMeanRate(j) = NaN; runBurstI(j) = NaN;
        D1ratenow{j} = []; evSess{j} = []; xbin{j} = []; occutime{j} = []; n1Dact(j) = NaN;
        if (strcmp(evType{j}, 'run')) 
            [spiketime, epid] = SpikeEventFilter(timeunit*data.spike.spiketime{i}, evTime{j});
            %%%%locate event position data
            evSess{j} = identifysession(evTime{j}, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
            posid = []; evid = [];
            if (~isempty(evSess{j}))
                posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess{j}) );
            end
            if numel(posid == 1)
                if (isfield(behav.general, 'eventname'))
                    evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
                else
                    evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
                end
            end
            if (numel(posid)~=1)|(numel(evid)~=1)
               disp(['-------------> field not computed: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
            else
               xbin{j} = bhdata.event.Xbin{posid}{evid}; lapoccuptime = bhdata.event.LapOccuptime{posid}{evid}; %%%%occutime{j} = bhdata.event.Occuptime{posid}{evid}; 
               lappostime = bhdata.event.LapAllPostimestamp{posid}{evid}; lapmeanspeed = bhdata.event.LapMeanSpeed{posid}{evid};
               lapx = bhdata.event.LapAllX{posid}{evid};
               [D1ratenow{j}, n1Dact(j), occutime{j}, runMeanRate(j), runBurstI(j)] = ...
                   getlinearmap(spiketime, evTime{j}, lappostime, lapx, xbin{j}, lapoccuptime, lapmeanspeed, MinLapMeanSpeed, MinLapNum, smParm);
            end
        end
    end
    data.GRfield.evt1DRateMaps{i} = D1ratenow; pinfo.GRfield.runMeanRate{i} = runMeanRate; pinfo.GRfield.runBurstI{i} = runBurstI;
    
    %%%%%%%%%%compute 1D rate map properties
    nnnev = numel(evTime); pinfo.GRfield.runSptlInfo{i} = NaN*ones(1,nnnev); pinfo.GRfield.runSptlInfotime{i} = NaN*ones(1,nnnev);
    pinfo.GRfield.runSparsty{i} = NaN*ones(1,nnnev); pinfo.GRfield.runActvArea{i} = NaN*ones(1,nnnev);
    for (j = 1:numel(evTime))
        if (strcmp(evType{j}, 'run')) 
            [spinfo, spinfotime, sparsty] = findratemapprop(D1ratenow{j}, occutime{j});
            pinfo.GRfield.runSptlInfo{i}(j) = spinfo; pinfo.GRfield.runSptlInfotime{i}(j) = spinfotime; 
            pinfo.GRfield.runSparsty{i}(j) = sparsty; pinfo.GRfield.runActvArea{i}(j) = n1Dact(j);
        end
    end
    
    %%%%compute evt 1D map stability and field directionality
    [dirpairs, sessevpairs, dirpairnames, sessevpairnames] = pairupevts(evName, evSess, evType);  
    pinfo.GRfield.run1DTrajPairs{i} = dirpairnames; 
    for (tt = 1:numel(dirpairs)) %%%directionality: shift rate maps by # bins, then compute the overlap
        rate1 = D1ratenow{dirpairs{tt}(1)}; rate2 = D1ratenow{dirpairs{tt}(2)};
        xx1 = xbin{dirpairs{tt}(1)}; xx2 = xbin{dirpairs{tt}(2)}; centernow = (min(xx2)+max(xx2))/2;
        crr = correlateratemaps(rate1, rate2, xx1, xx2, centernow, 1); %assume back/forth traj linearization reversed 
                                                                                  %centernow, 1); if x mirror
        [overlap, locshift] = computedirectionoverlap(rate1, rate2, xx1, xx2, dirNpix);
        pinfo.GRfield.run1DDirCrr{i}(tt)= crr; pinfo.GRfield.run1DDirOverlap{i}(tt) = overlap; pinfo.GRfield.run1DDirLocShift{i}(tt) = locshift; 
        if (max(rate1) + max(rate2))>0
            pinfo.GRfield.run1DDirRateDiffI{i}(tt) = abs(max(rate1)-max(rate2))/abs(max(rate1)+max(rate2));
        else
            pinfo.GRfield.run1DDirRateDiffI{i}(tt) = NaN;
        end
    end
    pinfo.GRfield.run1DEvtPairs{i} = sessevpairnames; %%%%stability
    for (tt = 1:numel(sessevpairs)) %%%stability
        rate1 = []; xx1 = []; rate2 = []; xx2 = []; maa = 0; mbb = 0;
        for (k = 1:numel(sessevpairs{tt}{1}))
            %disp([size(D1ratenow{sessevpairs{tt}{1}(k)}) 0 size(xbin{sessevpairs{tt}{1}(k)})]);
            rate1 = [rate1; D1ratenow{sessevpairs{tt}{1}(k)}]; xx1 = [xx1 maa+xbin{sessevpairs{tt}{1}(k)}]; maa = max(xx1);
        end
        for (k = 1:numel(sessevpairs{tt}{2}))
            rate2 = [rate2; D1ratenow{sessevpairs{tt}{2}(k)}]; xx2 = [xx2 mbb+xbin{sessevpairs{tt}{2}(k)}]; mbb = max(xx2);
        end
        [crr, pp] = correlateratemaps(rate1, rate2, xx1, xx2, 0, 0);
        pinfo.GRfield.run1DStablty{i}(tt) = crr; pinfo.GRfield.run1DStabP{i}(tt) = pp;
    end
  
    %%%%compute 1D place fields
    nfield = 0; baserate = NaN*ones(numel(evTime),1);
    for (tt = 1:numel(evTime)) 
        nx = numel(xbin{tt});
        if (nx > 1)
        baserate(tt) = computebaserate(D1ratenow{tt}, base1DrateN);
        pf = find1Dfieldprop(D1ratenow{tt}, xbin{tt}, threshold1Drate, max1Dgap, min1Dpeakrate, baserate(tt));
        for (j = 1:numel(pf))
             nfield = nfield + 1;
             pinfo.GRfield.PF1DBoundStart{i}(nfield) = pf(j).boundstart; pinfo.GRfield.PF1DBoundEnd{i}(nfield) = pf(j).boundend;
             pinfo.GRfield.PF1DLength{i}(nfield) = pf(j).leng; pinfo.GRfield.PF1DSize{i}(nfield) = pf(j).size; 
             pinfo.GRfield.PF1DInMeanrate{i}(nfield) = pf(j).inmeanrate; pinfo.GRfield.PF1DInPeakrate{i}(nfield) = pf(j).inpeakrate; 
             pinfo.GRfield.PF1DSkewness{i}(nfield) = pf(j).skew; pinfo.GRfield.PF1DKurtosis{i}(nfield) = pf(j).kurtosis;
             pinfo.GRfield.PF1DLocPeakX{i}(nfield) = pf(j).peakX; pinfo.GRfield.PF1DLocComX{i}(nfield) = pf(j).comX;
             pinfo.GRfield.PF1DLocStartX{i}(nfield) = pf(j).startX; pinfo.GRfield.PF1DLocEndX{i}(nfield) = pf(j).endX;
             pinfo.GRfield.PF1Devt{i}{nfield} = evName{tt};
        end
        end
    end
    pinfo.GRfield.PF1DNfield{i} = nfield; pinfo.GRfield.PF1DBaseRate{i} = baserate; disp(['-----------> number of 1D fields identified: ', num2str(nfield)]);
end
%%%%%%%%%%%%%%%END OF THE MAIN ROUTINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pf = find1Dfieldprop(rateingrid, xbin, zerothreshold, maxgap, minpeakrate, baselinerate)
pf = [];
%%algorithm: ubtract from baseline firing rate (for ctx), then examine maximum peak layer by layer 
%%break spike firing area into a number of place fields
%%work only with lineared track
%%input: rateingrid - spike firing rate in each grid on lineared track grid
%%       zerothreshold - a threshold (in percentage of average firing rate)
%%                    below which the area is not in its place field
%%       maxgap - maximum gap (rate below zero threshold) in a number of pixels within a place field
%%output: numfield - number of place fields
numfield = 0; binsize = xbin(2)-xbin(1); maxgap = round(maxgap/binsize); %%maxgap now in number of bins
gridnum = numel(rateingrid); gridleftindex = 1:gridnum; %%%this is the leftover indices in grid for new field identification
rateingrid = rateingrid - baselinerate; %iii = find(rateingrid<0); if (~isempty(iii)) rateingrid(iii) = zeros(1, numel(iii)); end
[peakratenow, peakindexnow] = max(rateingrid);
while (peakratenow >= minpeakrate) 
    startindex = []; endindex = []; threstartindex = []; threendindex = []; %%%field locations in rateingrid
    %%%%%determine the left bound of the current field
    zeropoint = 0;
    if (peakindexnow == 1)
        startindex = 1; threstartindex = 1;
    else
        for (k = peakindexnow-1:-1:1)
            if (rateingrid(k) < peakratenow * zerothreshold) % | (isnan(rateingrid(k))) 
               zeropoint = zeropoint + 1; 
            else
               zeropoint = 0;
            end
            if (zeropoint > maxgap)  
               startindex = k + floor(zeropoint/2); threstartindex = k + zeropoint; break;
            elseif ((~ismember(k, gridleftindex))|(k==1)) & (zeropoint >= floor(maxgap/2))
               startindex = k + zeropoint-floor(maxgap/2); threstartindex = k + zeropoint; break;
            elseif ((~ismember(k, gridleftindex))|(k==1)) & (zeropoint < floor(maxgap/2))
               startindex = k; threstartindex = k + zeropoint; break;   
            end
        end
    end
    %%%%determine the right bound of the current field
    zeropoint = 0;
    if (peakindexnow == gridnum)
        endindex = gridnum; threendindex = gridnum;
    else
        for (k = peakindexnow+1:gridnum)
            if (rateingrid(k) < peakratenow * zerothreshold) % | (isnan(rateingrid(k))) 
               zeropoint = zeropoint + 1; 
            else
               zeropoint = 0;
            end
            if (zeropoint >= maxgap)  
               endindex = k - floor(zeropoint/2); threendindex = k - zeropoint; break;
            elseif ((~ismember(k, gridleftindex))|(k==gridnum)) & (zeropoint >= floor(maxgap/2))
               endindex = k - zeropoint+floor(maxgap/2); threendindex = k - zeropoint; break;
            elseif ((~ismember(k, gridleftindex))|(k==gridnum)) & (zeropoint < floor(maxgap/2))
               endindex = k; threendindex = k - zeropoint; break;   
            end
        end
    end
    if (threendindex-threstartindex > 2)
        rrrnow = rateingrid(threstartindex:threendindex); rrrnow = rrrnow(~isnan(rrrnow));
        [~,ttt] = max(rrrnow); aaa = 0;
        if (~isempty(ttt)) aaa = abs(ttt-1)/numel(rrrnow); end
        %if ~isempty(ttt) && (ttt>2) && (ttt < numel(rrrnow)-1)
        if (aaa>0.2) && (aaa<0.8) %%%this is to check whether the fields is severely truncated
           numfield = numfield + 1;
           %%%%%%characteristics of the current field
           pf(numfield).peakX = xbin(peakindexnow); pf(numfield).startX = xbin(threstartindex); 
           pf(numfield).endX = xbin(threendindex);
           pf(numfield).boundstart = xbin(startindex); pf(numfield).boundend = xbin(endindex);
           pf(numfield).leng = (threendindex-threstartindex+1)*binsize; pf(numfield).size = binsize*sum(rateingrid(threstartindex:threendindex));
           pf(numfield).inmeanrate = baselinerate + mean(rateingrid(threstartindex:threendindex));
           pf(numfield).inpeakrate = peakratenow + baselinerate;
           %%%the center location index and skewness: ATTENTION: rateingrid is a (unnormalized) density function, but not samples
           locmean = NaN; skew = NaN; kurt = NaN;
           sumrate = sum(rateingrid(threstartindex:threendindex));
           if (sumrate >0)
               allloc = xbin(threstartindex:threendindex);  %%field location points
               densfun = rateingrid(threstartindex:threendindex)/sumrate;  %now normalized propability density function
               locmean = allloc * densfun; %sum(allloc .* densfun); 
               locvar = sqrt( ((allloc-locmean).^2) * densfun ); %sqrt( sum( ((allloc-locmean).^2) .* densfun ) );
               if (locvar > 0)
                  skew = ((allloc-locmean).^3) * densfun  / (locvar^3);
                  X = (allloc-locmean)/locvar; kurt = (X.^4)*densfun - 3;
                  %skew = sum( ((allloc-locmean).^3) .* densfun ) / (locvar^3);
                  %X = (allloc-locmean)/locvar; kurt = sum((X.^4).*densfun) - 3;
               end
           end
           pf(numfield).comX = locmean; pf(numfield).skew = skew; pf(numfield).kurtosis = kurt;
        end
    end
    %%%%%%leave leftovers for computinng the next field
    gridleftindex = setdiff(gridleftindex, startindex:endindex); %%%leftover grid indices
    ratenow = rateingrid(gridleftindex); [peakratenow, peaknewindex] = max(ratenow); peakindexnow = gridleftindex(peaknewindex);
end

function [overlap, locshift] = computedirectionoverlap(rate1, rate2, xx1, xx2, dirNpix)
overlap = NaN; locshift = NaN; %%%%assume xx1 and xx2 are physically reverse
if (numel(xx1)>1) && (numel(xx2)>1)
binsize = xx1(2)-xx1(1); dirNbin = round(dirNpix/binsize);
% size(rate1)
% size(rate2)
% disp([rate1 rate2]);
%%%%re-align the 1D grids
xx2 = xx2 - (max(xx2)+min(xx2))/2; xx1 = xx1 - (min(xx1)+max(xx1))/2; xx2 = -xx2;
minx = max([min(xx1) min(xx2)]); maxx = min([max(xx1) max(xx2)]);
gg = minx:binsize:maxx; ind1 = zeros(1, numel(gg)); ind2 = zeros(1, numel(gg)); %%%%re-grid the space
    for (j = 1:numel(gg))
        [mmm, kkk] = min(abs(xx1-gg(j))); ind1(j) = kkk;
        [mmm, kkk] = min(abs(xx2-gg(j))); ind2(j) = kkk;
    end
rate1 = rate1(ind1); rate2 = rate2(ind2);
% size(rate1)
% size(rate2)
% disp([rate1 rate2]);
%%%%%normalize
ind = (rate1>-10) & (rate2>-10); rate1 = rate1(ind); rate2 = rate2(ind);
if ((max(rate1)>=1) & (max(rate2)==0)) | ((max(rate2)>=1) & (max(rate1)==0))
    overlap = 0;
elseif (max(rate1) >= 1) & (max(rate2) >= 1)
    rate1 = rate1/mean(rate1); rate2 = rate2/mean(rate2);
    [x,y] = size(rate1); if (x>y) rate1 = rate1'; end
    [x,y] = size(rate2); if (x>y) rate2 = rate2'; end
    loc = -dirNbin:dirNbin; ov = zeros(size(loc)); 
    for (s = 1:numel(loc))
        ind = 1:numel(rate1); kk = find( (ind+s>=1) & (ind+s<=numel(rate2)) ); 
        r1now = rate1(ind(kk)); r2now = rate2(ind(kk)+s);
        if (numel(kk)>0)
            %disp(['-----------> number of overlapping bins: ', num2str(numel(kk))]);
            ov(s) = sum( min([r1now; r2now]) )/numel(kk);
        end
    end
    [overlap, tt] = max(ov); locshift = loc(tt)*binsize; 
end
end

function [crr,pp] = correlateratemaps(rate1, rate2, XY1, XY2, RC, theta)
%%%%here RC and theta are the rotation center and angle: RC[x y](coordinates in pixels), 
%%%% theta (angles in degree): 0, 90, -90, 180, 1 (x mirror), -1 (y mirror)  
%%%% assume XYgrid1 and XYgrid2 (could be 1D or 2D) have the same bin sizes
crr = NaN; pp = NaN;
if (numel(XY1)>1) && (numel(XY2)>1)
%%%%first transform the rotated rate maps for 2D maps; for 1D maps: only do the x mirror (only choice 1 is valid... wait... look at directionaly algorithm first)
%%%%for 2D maps: transform XY grid coordinates, then choose common grids and their rates
%%%%%%%%%%determine 1D or 2D rate maps
if (~iscell(XY1)) %%if a 1D rate map: add the cell dimension and RC input replaced by traj2 center
    nd = 1; binsize = XY2(2) - XY2(1);
else
    nd = 2; binsize = XY2{1}(2) - XY2{1}(1);
end
%%%%%first, this is the first map rates and coordinates (x,y), & the second map without rotation/mirror, shifting, etc.
if (nd == 1)
    r1 = rate1; x1 = XY1; r2 = rate2; x2 = XY2; RC = (min(XY2) + max(XY2))/2; 
elseif (nd == 2)
    nx1 = numel(XY1{1}); ny1 = numel(XY1{2});
    r1 = reshape(rate1, nx1*ny1, 1); x1 = zeros(nx1*ny1,1); y1 = zeros(nx1*ny1,1); 
    for (i = 1:nx1)
        for (j = 1:ny1)
            x1((i-1)*ny1+j) = XY1{1}(i); y1((i-1)*ny1+j) = XY1{2}(j);
        end
    end
    nx2 = numel(XY2{1}); ny2 = numel(XY2{2});
    r2 = reshape(rate2, nx2*ny2, 1); x2 = zeros(nx2*ny2,1); y2 = zeros(nx2*ny2,1); 
    for (i = 1:nx2)
        for (j = 1:ny2)
            x2((i-1)*ny2+j) = XY2{1}(i); y2((i-1)*ny2+j) = XY2{2}(j);
        end
    end
end
%%%%%%second: work on the second coordinates (rotated, mirrored, shifted, etc)
%%%%rotation or mirror transformation: shift cartesian origin to RC, then rotate in polar coordinates, then shift cartesian orgin back
%%%%%%%%%%%%shift origin
x2 = x2 - RC(1); if (nd == 2) y2 = y2-RC(2); end
%%%%%%%%%%%%rotate or mirror
if (abs(theta) > 2) %if rotate 2D maps
    [ang, rr] = cart2pol(x2, y2); ang = ang + theta*2*pi/360; [x2, y2] = pol2cart(ang, rr); %coordinates for all rate2 points
elseif (theta == 1) %if mirror X
    x2 = -x2;
elseif (theta == -1) %if mirror y: only available for 2D maps
    if (nd ==2) y2 = -y2; end
end
%%%%%%%%%%%%shift origin back
if (nd ==1) x2 = x2 + RC; end
if (nd == 2) x2 = x2+RC(1); y2 = y2+RC(2); end
%%%%match up (r1, x1, y1) with (r2, x2, y2): coordinate (x y) need to be within binsize 
iii = find(r1>-10); r1 = r1(iii); x1 = x1(iii); if (nd ==2) y1 = y1(iii); end
jjj = find(r2>-10); r2 = r2(jjj); x2 = x2(jjj); if (nd ==2) y2 = y2(jjj); end
nmatch = 0; rr1 = []; rr2 = []; %%%%have to work point-by-point
if (nd ==1)
    for (i = 1:numel(x1))
        [mmm, iii] = min( abs(x2-x1(i)) );
        if (mmm <= binsize)
            nmatch = nmatch + 1; rr1(nmatch) = r1(i); rr2(nmatch) = r2(iii);
        end
    end
elseif (nd == 2)
    for (i = 1:numel(x1))
        [mmm, iii] = min( abs(x2-x1(i)) + abs(y2-y1(i)) );
        if (mmm <= 2*binsize)
            nmatch = nmatch + 1; rr1(nmatch) = r1(i); rr2(nmatch) = r2(iii);
        end
    end
end
if (nmatch>=10) %%%%if overlap more than 10 bins
    m1 = mean(rr1); m2 = mean(rr2); s1 = std(rr1); s2 = std(rr2);
    if (s1~=0)&&(s2~=0)
        [R,P] = corrcoef(rr1',rr2'); crr = R(1,2); pp = P(1,2);
        %rr1 = (rr1-m1)/s1; rr2 = (rr2-m2)/s2; crr = rr1*rr2'/nmatch; 
    end
end
end


function [dirpairs, sessevpairs, dirpairnames, sessevpairnames] = pairupevts(evName, evSess, evType)
%%%dirpairs{npair}[1,2], dirpairnames{npair}{1,2} = traj pairs' indices and names for computing directionality 
%%%sessevpairs{npair}{1,2}, sessevpairnames{npair}{1,2} = two groups of events' indices and names for computing stability 
dirpairs = []; sessevpairs = []; dirpairnames = []; sessevpairnames = [];
ndirpair = 0; nevpair = 0; allsess = unique(evSess(find(strcmp(evType, 'run')))); 
for (i = 1:numel(allsess))
    evind{i} = find(strcmp(evSess, allsess{i})); sessEvName{i} = evName(evind{i}); %%all events in a session
    for (tt = 1:numel(sessEvName{i}))
        for (ss = tt+1:numel(sessEvName{i}))
            ndirpair = ndirpair + 1; dirpairnames{ndirpair} = strcat(sessEvName{i}{tt}, '_', sessEvName{i}{ss});
            dirpairs{ndirpair} = [evind{i}(tt) evind{i}(ss)];
        end
    end
end
for (i = 1:numel(allsess))
    for (j = i+1:numel(allsess))
        [pairnow, pairnames, ok] = ifmatches(sessEvName{i}, evind{i}, sessEvName{j}, evind{j});
        if (ok)
            nevpair = nevpair + 1;
            sessevpairs{nevpair} = pairnow; sessevpairnames{nevpair} = pairnames;
        end
    end
end

function [pairnow, pairnames, ok] = ifmatches(sessEvName1, evind1, sessEvName2, evind2)
pairnow{1} = [] ; pairnow{2} = []; pairnames{1} = []; pairnames{2} = []; nelm = 0; ok = 1;
%if (numel(evind1) ~= numel(evind2))
%    ok = 0;
%else
    [all2ndstr, all2ndtraj] = strtok(sessEvName2, '_');
    for (i = 1:numel(evind1))
        [str, trajnow] = strtok(sessEvName1{i}, '_'); 
        jj = find(strcmp(all2ndtraj, trajnow));
        if (numel(jj) == 1)
            nelm = nelm + 1; pairnow{1}(nelm) = evind1(i); pairnow{2}(nelm) = evind2(jj);
            pairnames{1} = strcat(pairnames{1},'+', sessEvName1{i}); pairnames{2} = strcat(pairnames{2},'+', sessEvName2{jj});
        end
    end
%end
if (isempty(pairnow{1})) ok = 0; end
    
function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end

function [D1rate, nact, occutime, meanrate, burstI] = getlinearmap(spiketime, evTime, lappostime, lapx, xbin, lapoccuptime,...
        lapmeanspeed, MinLapMeanSpeed, MinLapNum, smParm)
nx = numel(xbin); D1rate = NaN*ones(nx,1); nlap = numel(lappostime); count = zeros(nx,1); occutime = zeros(nx,1); spikex = []; nact = NaN;
meanrate = NaN; burstI = NaN; alltime = 0; spikeT = [];
if (nx > 1) && (numel(find(lapmeanspeed>=MinLapMeanSpeed))>=MinLapNum);
    binsize = xbin(2) - xbin(1); sigma = round(smParm.d1sigma/binsize); %%sigma now in number of bins
%%%%%count spikes in xbin lap bu lap
for (i = 1:nlap)
    if (lapmeanspeed(i)>=MinLapMeanSpeed)
       counti = zeros(nx,1);
       spikenow = sort( spiketime( (spiketime>=evTime.start(i)) & (spiketime<=evTime.ent(i)) ) ); spikeT = [spikeT; spikenow];
       spikex = NaN*ones(size(spikenow)); [laptimenow, iii] = sort(lappostime{i}); allxnow = lapx{i}(iii);
       lastpoint = 1;  %this only for saving time
       for (j = 1:numel(spikenow))
         for (k = lastpoint:numel(laptimenow)-1) %find corresponding time in position data
             if (laptimenow(k) <= spikenow(j)) && (laptimenow(k+1) > spikenow(j)) 
                 spikex(j) = allxnow(k); lastpoint = k; break; 
             end
         end
       end
       for (j = 1:nx-1) counti(j) = numel(find((spikex>=xbin(j)) & (spikex<xbin(j+1)))); end
       counti(nx) = numel(find(spikex>=xbin(nx)));
       count = count + counti;
       occutime = occutime + lapoccuptime{i}; alltime = alltime + sum(lapoccuptime{i});
    end
end
nnnT = numel(spikeT);
if (alltime > 0) meanrate = nnnT/alltime; end
if (nnnT>1) burstI = findBurstI(spikeT, smParm.BIint); end
%%%%%compute firing rate
for (i = 1:nx)
     if (occutime(i) ~= 0) D1rate(i) = count(i)/occutime(i); end
end
nact = numel(find(count>0))/numel(find(occutime>0));
%%now do a 1d smoothing
if (strcmpi(smParm.d1sm, 'yes'))
    D1rate = OneDSmooth_new(D1rate, occutime, sigma, smParm.d1Nsig);
    %XX = [-smParm.d1Nsig*sigma:1:smParm.d1Nsig*sigma]; 
    %wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
    %D1rate = conv(D1rate, wind); D1rate = D1rate((nw-1)/2+1:numel(D1rate)-(nw-1)/2);
end
end

function [spinfo, spinfotime, sparsty] = findratemapprop(rate, occu)
%%%%%%another variable: mutual information is not computed here
spinfo = NaN; spinfotime = NaN; sparsty = NaN;  %size(rate)
%size(occu)
[ny, nx] = size(rate); rate = reshape(rate, nx*ny, 1); occu = reshape(occu, nx*ny, 1); mr = max(rate);
if (sum(occu) > 0)
   occu = occu / sum(occu); jjj = find(rate>-1); meanrate = sum(occu(jjj).*rate(jjj));
   spinfo = 0; spinfotime = 0; sparsty = 0;  
   if (meanrate > 0)
       rateok = rate / meanrate; iii = find(rateok > 0); 
       spinfo = sum(occu(iii).*rateok(iii).*log2(rateok(iii))); 
       spinfotime = spinfo*meanrate; sparsty = 1/sum(occu(iii).*rateok(iii).*rateok(iii));
   end
end

function baserate = computebaserate(ratenow, baserategridN)
baserate = 0;
ratenow = sort(ratenow(ratenow>-10)); nn = round(baserategridN*numel(ratenow));
if (nn > 0) baserate = ratenow(nn); end
%if (nn > 0) baserate = mean(ratenow(1:nn)); end

function [writefilename, okk] = getoutputfile(hf, ow)
okk = 1; writefilename = [];
if (ow == 0)
   [fname, pname] = uiputfile(fullfile(cd, '*.spikedb'), 'Write the new spike database to:');
   if (numel(fname)>1)
      writefilename = fullfile(pname, fname);
   else
      okk = 0;
   end
else
   %input = questdlg('The current database will be altered and overwritten. Are you sure?', 'Overwrite?', 'Yes');
   %if (strcmp(input, 'Yes'))
      fname = get(hf, 'Name'); ftt = strfind(fname, '__'); writefilename = fname(ftt+2:numel(fname));
   %else
   %   okk = 0;
   %end
end

function BI = findBurstI(spikeT, BIint)
BI = NaN; totalspike = numel(spikeT);
if (totalspike>=10) 
    spikeint = diff(sort(spikeT)); 
    burstspike = numel(find(spikeint<=BIint));
    BI = burstspike/totalspike;
end
