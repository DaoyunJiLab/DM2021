function [pinfo,data] = DataManager_FindSpatialInfoFieldProp(pinfo,data,behav,bhdata,cellind,vv)
%%%Compute the spatial info and place field properties:
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

%variable to assign
if (~isempty(cellind))
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo, 'field')) pinfo.field = []; end
   %%%session field properties
   if (~isfield(pinfo.field, 'sessSptlInfo')) pinfo.field.sessSptlInfo = cell(1, nspike); end
   if (~isfield(pinfo.field, 'sessSptlInfotime')) pinfo.field.sessSptlInfotime = cell(1, nspike); end
   if (~isfield(pinfo.field, 'sessSparsty')) pinfo.field.sessSparsty = cell(1, nspike); end
   if (~isfield(pinfo.field, 'sessActvArea')) pinfo.field.sessActvArea = cell(1, nspike); end
   if (~isfield(pinfo.field, 'sessPairs')) pinfo.field.sessPairs = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'sessPairMinRate')) pinfo.field.sessPairMinRate = cell(1, nspike); end
   if (~isfield(pinfo.field, 'sessPairMaxRate')) pinfo.field.sessPairMaxRate = cell(1, nspike); end
   if (~isfield(pinfo.field, 'sess2DStablty')) pinfo.field.sess2DStablty = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'sess2DStabP')) pinfo.field.sess2DStabP = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'sess2DRotStablty')) pinfo.field.sess2DRotStablty = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'sess2DRotStabP')) pinfo.field.sess2DRotStabP = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'sess2DSlideShufStablty')) pinfo.field.sess2DSlideShufStablty = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'sess2DSlideShufRotStablty')) pinfo.field.sess2DSlideShufRotStablty = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'sess2DSlideShufStabltyZ')) pinfo.field.sess2DSlideShufStabltyZ = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'sess2DSlideShufRotStabltyZ')) pinfo.field.sess2DSlideShufRotStabltyZ = cell(1, nspike); end 
   
   if (~isfield(pinfo.field, 'PF2DNfield')) pinfo.field.PF2DNfield = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PFsessMeanRate')) pinfo.field.PFsessMeanRate = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF2Dsess')) pinfo.field.PF2Dsess = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF2DSize')) pinfo.field.PF2DSize = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF2DArea')) pinfo.field.PF2DArea = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF2DpeakX')) pinfo.field.PF2DpeakX = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF2DpeakY')) pinfo.field.PF2DpeakY = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF2DInMeanrate')) pinfo.field.PF2DInMeanrate = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'PF2DInPeakrate')) pinfo.field.PF2DInPeakrate = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'PF2DBaseRate')) pinfo.field.PF2DBaseRate = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'PF2DSessMeanRate')) pinfo.field.PF2DSessMeanRate = cell(1, nspike); end 
   
   %%%evt field properties
   if (~isfield(pinfo.field, 'runSptlInfo')) pinfo.field.runSptlInfo = cell(1, nspike); end
   if (~isfield(pinfo.field, 'runSptlInfotime')) pinfo.field.runSptlInfotime = cell(1, nspike); end
   if (~isfield(pinfo.field, 'runSparsty')) pinfo.field.runSparsty = cell(1, nspike); end   
   if (~isfield(pinfo.field, 'runActvArea')) pinfo.field.runActvArea = cell(1, nspike); end
   if (~isfield(pinfo.field, 'run1DEvtPairs')) pinfo.field.run1DEvtPairs = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'run1DEvtPairMaxRate')) pinfo.field.run1DEvtPairMaxRate = cell(1, nspike); end
   if (~isfield(pinfo.field, 'run1DEvtPairMinRate')) pinfo.field.run1DEvtPairMinRate = cell(1, nspike); end
   if (~isfield(pinfo.field, 'run1DEvtPairRateDiffI')) pinfo.field.run1DEvtPairRateDiffI = cell(1, nspike); end
   if (~isfield(pinfo.field, 'run1DStablty')) pinfo.field.run1DStablty = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'run1DStabP')) pinfo.field.run1DStabP = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'run1DSlideShufStablty')) pinfo.field.run1DSlideShufStablty = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'run1DSlideShufStabltyZ')) pinfo.field.run1DSlideShufStabltyZ = cell(1, nspike); end 
   %if (~isfield(pinfo.field, 'run1DIDShufStablty')) pinfo.field.run1DIDShufStablty = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'run1DTrajPairs')) pinfo.field.run1DTrajPairs = cell(1, nspike); end  %directionality 
   if (~isfield(pinfo.field, 'run1DTrajPairMaxRate')) pinfo.field.run1DTrajPairMaxRate = cell(1, nspike); end
   if (~isfield(pinfo.field, 'run1DTrajPairMinRate')) pinfo.field.run1DTrajPairMinRate = cell(1, nspike); end
   if (~isfield(pinfo.field, 'run1DDirCrr')) pinfo.field.run1DDirCrr = cell(1, nspike); end  %directionality: simple correlation 
   if (~isfield(pinfo.field, 'run1DDirCrrP')) pinfo.field.run1DDirCrrP = cell(1, nspike); end
   if (~isfield(pinfo.field, 'run1DDirCrossCrrPcr')) pinfo.field.run1DDirCrossCrrPcr = cell(1, nspike); end
   if (~isfield(pinfo.field, 'run1DDirCrossCrrPloc')) pinfo.field.run1DDirCrossCrrPloc = cell(1, nspike); end
   if (~isfield(pinfo.field, 'run1DDirCrossCrrPpv')) pinfo.field.run1DDirCrossCrrPpv = cell(1, nspike); end
   if (~isfield(pinfo.field, 'run1DDirOverlap')) pinfo.field.run1DDirOverlap = cell(1, nspike); end  %directionality: overlap
   if (~isfield(pinfo.field, 'run1DDirLocShift')) pinfo.field.run1DDirLocShift = cell(1, nspike); end  %directionality: location shift 
   if (~isfield(pinfo.field, 'run1DDirRateDiffI')) pinfo.field.run1DDirRateDiffI = cell(1, nspike); end  %rate difference index for each field 
   
   if (~isfield(pinfo.field, 'PF1DNfield')) pinfo.field.PF1DNfield = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PFevtMeanRate')) pinfo.field.PFevtMeanRate = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1Devt')) pinfo.field.PF1Devt = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DBoundStart')) pinfo.field.PF1DBoundStart = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DBoundEnd')) pinfo.field.PF1DBoundEnd = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DLocPeakX')) pinfo.field.PF1DLocPeakX = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DLocComX')) pinfo.field.PF1DLocComX = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DLocStartX')) pinfo.field.PF1DLocStartX = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DLocEndX')) pinfo.field.PF1DLocEndX = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DLength')) pinfo.field.PF1DLength = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DSize')) pinfo.field.PF1DSize = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DSkewness')) pinfo.field.PF1DSkewness = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DKurtosis')) pinfo.field.PF1DKurtosis = cell(1, nspike); end
   if (~isfield(pinfo.field, 'PF1DInMeanrate')) pinfo.field.PF1DInMeanrate = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'PF1DInPeakrate')) pinfo.field.PF1DInPeakrate = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'PF1DBaseRate')) pinfo.field.PF1DBaseRate = cell(1, nspike); end 
   if (~isfield(pinfo.field, 'PF1DTrajMeanRate')) pinfo.field.PF1DTrajMeanRate = cell(1, nspike); end 
      
   %%%data variables to assign
   if (~isfield(data, 'field')) data.field = []; end
   if (~isfield(data.field, 'sess2DRateMaps')) data.field.sess2DRateMaps = cell(1, nspike); end
   if (~isfield(data.field, 'evt1DRateMaps')) data.field.evt1DRateMaps = cell(1, nspike); end
end

for (jjjk = 1:numel(cellind))
    i = cellind(jjjk); nshuffle = 100;
    minlap = pinfo.parm.f1DMinLapNum(i); 
    threshold1Drate = pinfo.parm.f1DThresholdRate(i); max1Dgap = pinfo.parm.f1DMaxGap(i); min1Dpeakrate = pinfo.parm.f1DMinPeakRate(i); 
    timeunit = pinfo.parm.timeunit(i); base1DrateN = pinfo.parm.f1DBaseRateGridN(i); 
    threshold2Drate = pinfo.parm.f2DThresholdRate(i); max2Dgap = pinfo.parm.f2DMaxGap(i); min2Dpeakrate = pinfo.parm.f2DMinPeakRate(i); 
    base2DrateN = pinfo.parm.f2DBaseRateGridN(i); 
    smParm.d2sigma = pinfo.parm.fSmooth2DSigma(i); smParm.d2Nsig = pinfo.parm.fSmooth2DNSigma(i); smParm.d2sm = pinfo.parm.fSmooth2D{i};
    smParm.d1sigma = pinfo.parm.fSmooth1DSigma(i); smParm.d1Nsig = pinfo.parm.fSmooth1DNSigma(i); smParm.d1sm = pinfo.parm.fSmooth1D{i};
    RC(1) = pinfo.parm.fRotate2DCenX(i); RC(2) = pinfo.parm.fRotate2DCenY(i); theta = pinfo.parm.fRotate2DAngle(i); 
    dirNpix = pinfo.parm.fDirShift(i); rotornot = pinfo.parm.fRotate{i};
    %disp(strcat('-----> compute field properties ---', pinfo.general.parmfile{i}));
    disp(['-----> compute field properties and spatial info (', num2str(jjjk), ' out of ', num2str(numel(cellind)), '): ', pinfo.general.parmfile{i}]);
    %%%%%compute session rate maps
    sessions = pinfo.general.sessionname{i}; nsess = numel(sessions); sesstype = pinfo.parm.sessType{i}; D2ratenow = []; XYbin = [];
    for (tt = 1:nsess)
        D2ratenow {tt} = []; gridoccup{tt} = []; XYbin{tt} = []; n2Dact(tt) = NaN;
        ep.start = pinfo.general.sessionstartT{i}(tt); ep.ent = pinfo.general.sessionendT{i}(tt); 
        [spiketime, epid] = SpikeEventFilter(timeunit*data.spike.spiketime{i}, ep);
        %%%session position data
        finaldirnow = pinfo.general.finaldir{i}; sessnow = sessions{tt}; 
        posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, sessnow) );
        if (numel(posid) ~= 1)
            disp(['-------------> field not computed: no or more than 1 positon files match the session: ', finaldirnow, '___', sessnow]);
        else
            postime = bhdata.pos.postimestamp{posid}*behav.parm.timeunit(posid); framerate = behav.parm.framerate(posid);
            XYbin{tt} = bhdata.sess.gridXYbin{posid}; gridoccup{tt} = bhdata.sess.gridOccuptime{posid}; 
            %%%find which color to use
            if ~isempty(behav.parm.sessPmarker{posid})
               cid = find( strcmp(behav.general.posMarker{posid}, behav.parm.sessPmarker{posid}) );    
               if ~isempty(cid)
                  XX = bhdata.pos.XX{posid}{cid}*behav.parm.pixelXSize(posid); 
                  YY = bhdata.pos.YY{posid}{cid}*behav.parm.pixelYSize(posid); 
               else
                  XX = []; YY = [];
               end
            else 
               XX = []; YY = [];
            end
            [D2ratenow{tt}, n2Dact(tt)] = get2Dratemaps(spiketime, XYbin{tt}, gridoccup{tt}, postime, XX, YY, smParm, 1/framerate);
        end
    end
    data.field.sess2DRateMaps{i} = D2ratenow;
    %%%%%compute session 2D map properties: spatial info, sparsity, active area
    for (tt = 1:nsess) 
        [spinfo, spinfotime, sparsty] = findratemapprop(D2ratenow{tt}, gridoccup{tt});
        pinfo.field.sessSptlInfo{i}(tt) = spinfo; pinfo.field.sessSptlInfotime{i}(tt) = spinfotime; 
        pinfo.field.sessSparsty{i}(tt) = sparsty; pinfo.field.sessActvArea{i}(tt) = n2Dact(tt);
    end
    %%%%compute session 2D map stability
    ncom = 0;
    for (tt = 1:nsess)
        [ny1, nx1] = size(D2ratenow{tt}); ratett = pinfo.firing.sessmeanrate{i}(tt);
        for (ss = tt+1:nsess)
            [ny2, nx2] = size(D2ratenow{ss}); ratess = pinfo.firing.sessmeanrate{i}(ss);
            if (strcmp(sesstype{tt}, sesstype{ss})) && ( (nx1>1)&&(ny1>1)&&(nx2>1)&&(ny2>1) ) ...
                    && strcmp(sesstype{tt}, 'open') && contains(sessions{tt}, 'boxA') && contains(sessions{ss}, 'boxA')
                ncom = ncom + 1; 
                [stb0, pp0] = correlateratemaps(D2ratenow{tt}, D2ratenow{ss}, XYbin{tt}, XYbin{ss}, RC, 0);
                if (strcmpi(rotornot, 'yes'))
                   [rotstb, rotpp] = correlateratemaps(D2ratenow{tt}, D2ratenow{ss}, XYbin{tt}, XYbin{ss}, RC, theta);
                else
                   rotstb = NaN; rotpp = NaN; 
                end
                pinfo.field.sessPairs{i}{ncom} = strcat(sessions{tt}, '_', sessions{ss});
                pinfo.field.sessPairMinRate{i}{ncom} = min([ratett ratess]);  pinfo.field.sessPairMaxRate{i}{ncom} = max([ratett ratess]);
                pinfo.field.sess2DStablty{i}{ncom} = stb0; pinfo.field.sess2DRotStablty{i}{ncom} = rotstb; 
                pinfo.field.sess2DStabP{i}{ncom} = pp0; pinfo.field.sess2DRotStabP{i}{ncom} = rotpp; 
                %%%%do a shuffled version
                scrr = zeros(1, nshuffle); SZ = NaN; rcrr = zeros(1, nshuffle); RZ = NaN; 
                for jjj = 1:nshuffle
                    scrr(jjj) = correlateratemaps(D2ratenow{tt}, D2slideshuffle(D2ratenow{ss}, gridoccup{ss}), XYbin{tt}, XYbin{ss}, RC, 0);
                    if (strcmpi(rotornot, 'yes'))
                       rcrr(jjj) = correlateratemaps(D2ratenow{tt}, D2slideshuffle(D2ratenow{ss}, gridoccup{ss}), XYbin{tt}, XYbin{ss}, RC, theta);
                    end
                end
                if std(scrr)~=0 SZ = (stb0-mean(scrr))/std(scrr); end
                if std(rcrr)~=0 RZ = (rotstb-mean(rcrr))/std(rcrr); end
                pinfo.field.sess2DSlideShufStablty{i}{ncom} = scrr(1); pinfo.field.sess2DSlideShufRotStablty{i}{ncom} = rcrr(1); 
                pinfo.field.sess2DSlideShufStabltyZ{i}{ncom} = SZ; pinfo.field.sess2DSlideShufRotStabltyZ{i}{ncom} = RZ; 
            end
        end
    end
    %%%%compute 2D place fields
    nfield = 0;
    for (tt = 1:nsess) 
        nttfield = 0; baserate = 0;
        [ny, nx] = size(D2ratenow{tt});
        if (nx > 1) && (ny > 1)
            baserate = computebaserate(D2ratenow{tt}, base2DrateN);
            pf = find2Dfieldprop(D2ratenow{tt}, XYbin{tt}, threshold2Drate, max2Dgap, min2Dpeakrate, baserate);
            for (j = 1:numel(pf))
             nfield = nfield + 1; nttfield = nttfield + 1;
             pinfo.field.PF2DSize{i}(nfield) = pf(j).size; pinfo.field.PF2DInMeanrate{i}(nfield) = pf(j).inmeanrate;
             pinfo.field.PF2DInPeakrate{i}(nfield) = pf(j).inpeakrate; 
             pinfo.field.PF2Dsess{i}{nfield} = sessions{tt}; pinfo.field.PF2DArea{i}(nfield) = pf(j).area; 
             pinfo.field.PF2DpeakX{i}(nfield) = pf(j).peakX; pinfo.field.PF2DpeakY{i}(nfield) = pf(j).peakY;
             pinfo.field.PF2DSessMeanRate{i}(nfield) = pinfo.firing.sessmeanrate{i}(tt); 
            end
        end
        pinfo.field.PF2DNfield{i}(tt) = nttfield; pinfo.field.PFsessMeanRate{i}(tt) = pinfo.firing.sessmeanrate{i}(tt); 
        pinfo.field.PF2DBaseRate{i}(tt) = baserate; 
    end
    disp(['-----------> number of 2D fields identified: ', num2str(nfield)]);
    
    %%%%%%%%compute event (linear track) 1D rate maps
    evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i}; 
    finaldirnow = pinfo.general.finaldir{i};
    D1ratenow = []; evSess = []; xbin = []; occutime = [];
    for (j = 1:numel(evTime))
        D1ratenow{j} = []; evSess{j} = []; xbin{j} = []; occutime{j} = []; n1Dact(j) = NaN;
        if (strcmp(evType{j}, 'run')) 
            evSess{j} = identifysession(evTime{j}, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
            posid = []; evid = [];
            if (~isempty(evSess{j}))
                posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess{j}) );
            end
            if numel(posid) == 1
                if (isfield(behav.general, 'eventname'))
                    evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
                else
                    evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
                end
            end
            if (numel(posid)~=1)||(numel(evid)~=1)
               disp(['-------------> 1D field not computed: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
            elseif contains(behav.parm.eventPosltr{posid}{evid}, '.ltr')  %if this is run event & can be linearized
               nevt = numel(evTime{j}.start);
                if (nevt<minlap)
                    disp(['-------------> 1D field not computed: not enough laps: ', finaldirnow, '___', evName{j}]);
                else
                    xbin{j} = bhdata.event.Xbin{posid}{evid}; occutime{j} = bhdata.event.Occuptime{posid}{evid};  
                    xjoint{j} = bhdata.event.Xjoint{posid}{evid}; framerate = behav.parm.framerate(posid);
                    [spiketime, ~] = SpikeEventFilter(timeunit*data.spike.spiketime{i}, evTime{j});
                    lappostime = bhdata.event.LapAllPostimestamp{posid}{evid}; 
                    lapx = bhdata.event.LapAllX{posid}{evid};
                    [D1ratenow{j}, n1Dact(j)] = getlinearmap(spiketime, evTime{j}, lappostime, lapx, xbin{j}, occutime{j}, smParm, 1/framerate);
                end
            end
        end
    end
    data.field.evt1DRateMaps{i} = D1ratenow;
    
    %%%%%%%%%%compute 1D rate map properties
    nnnev = numel(evTime); pinfo.field.runSptlInfo{i} = NaN*ones(1,nnnev); pinfo.field.runSptlInfotime{i} = NaN*ones(1,nnnev);
    pinfo.field.runSparsty{i} = NaN*ones(1,nnnev); pinfo.field.runActvArea{i} = NaN*ones(1,nnnev);
    for (j = 1:numel(evTime))
        if (strcmp(evType{j}, 'run')) 
            [spinfo, spinfotime, sparsty] = findratemapprop(D1ratenow{j}, occutime{j});
            pinfo.field.runSptlInfo{i}(j) = spinfo; pinfo.field.runSptlInfotime{i}(j) = spinfotime; 
            pinfo.field.runSparsty{i}(j) = sparsty; pinfo.field.runActvArea{i}(j) = n1Dact(j);
        end
    end
    
    %%%%compute evt 1D map stability and field directionality
   if dirNpix>=0
      [dirpairs, sessevpairs, dirpairnames, sessevpairnames] = pairupevts(evName, evSess, evType);  
      pinfo.field.run1DTrajPairs{i} = dirpairnames; 
    for (tt = 1:numel(dirpairs)) %%%directionality: shift rate maps by # bins, then compute the overlap
        pinfo.field.run1DDirRateDiffI{i}(tt) = NaN; pinfo.field.run1DDirCrr{i}(tt) = NaN;
        pinfo.field.run1DDirOverlap{i}(tt) = NaN; pinfo.field.run1DDirLocShift{i}(tt) = NaN; 
        pinfo.field.run1DDirCrrP{i}(tt) = NaN; pinfo.field.run1DDirCrossCrrPcr{i}(tt) = NaN; 
        pinfo.field.run1DDirCrossCrrPloc{i}(tt) = NaN; pinfo.field.run1DDirCrossCrrPpv{i}(tt) = NaN;
        %%%%added below for identifying active trajectories
        mmrate1 = pinfo.firing.evtmeanrate{i}(dirpairs{tt}(1)); mmrate2 = pinfo.firing.evtmeanrate{i}(dirpairs{tt}(2));
        pinfo.field.run1DTrajPairMaxRate{i}(tt) = max([mmrate1 mmrate2]); pinfo.field.run1DTrajPairMinRate{i}(tt) = min([mmrate1 mmrate2]);
        rate1 = D1ratenow{dirpairs{tt}(1)}; rate2 = D1ratenow{dirpairs{tt}(2)};
        if (~isempty(rate1)) && (~isempty(rate2))
           xx1 = xbin{dirpairs{tt}(1)}; xx2 = xbin{dirpairs{tt}(2)}; 
           %%%% revision: decide to change to a different alignment --- align at joint points
           xjoint1 = xjoint{dirpairs{tt}(1)}; xjoint2 = xjoint{dirpairs{tt}(2)}; 
           [rate1, rate2, x] = realignjointxx(xjoint1, xjoint2, rate1, rate2, xx1, xx2, 'yesflip', dirpairnames); %assume back/forth traj linearization reversed -- this is not true in many cases: like figure-8 maze
           [crr, pp] = correlateratemaps(rate1, rate2, x, x, 0, 0); %%here rate1/rate1 same direction
           [crossR, crossX, crossP] = Dircrosscrr(rate1, rate2, mean(diff(x)), dirNpix);
           pinfo.field.run1DDirCrossCrrPcr{i}(tt) = crossR; pinfo.field.run1DDirCrossCrrPloc{i}(tt) = crossX;
           pinfo.field.run1DDirCrossCrrPpv{i}(tt) = crossP;
           %%%% done revision
           %[overlap, locshift] = computedirectionoverlap(rate1, rate2, x, x, dirNpix);
           pinfo.field.run1DDirCrr{i}(tt)= crr; pinfo.field.run1DDirCrrP{i}(tt) = pp;
           %pinfo.field.run1DDirOverlap{i}(tt) = overlap; pinfo.field.run1DDirLocShift{i}(tt) = locshift; 
           if (max(rate1) + max(rate2))>0
               pinfo.field.run1DDirRateDiffI{i}(tt) = abs(max(rate1)-max(rate2))/abs(max(rate1)+max(rate2));
           else
               pinfo.field.run1DDirRateDiffI{i}(tt) = NaN;
           end
        end
    end
    pinfo.field.run1DEvtPairs{i} = sessevpairnames; %%%%stability
    for (tt = 1:numel(sessevpairs)) %%%stability
        pinfo.field.run1DStablty{i}(tt) = NaN; pinfo.field.run1DStabP{i}(tt) = NaN; pinfo.field.run1DEvtPairRateDiffI{i}(tt) = NaN;
        pinfo.field.run1DSlideShufStablty{i}(tt) = NaN; pinfo.field.run1DSlideShufStabltyZ{i}(tt) = NaN; 
        rate1 = []; xx1 = []; rate2 = []; xx2 = []; mmrate1 = []; mmrate2 = []; xjoint1 = []; xjoint2 = []; %%maa = 0; mbb = 0;
        for (k = 1:numel(sessevpairs{tt}{1}))
            if (~isempty(D1ratenow{sessevpairs{tt}{1}(k)})) && (~isempty(D1ratenow{sessevpairs{tt}{2}(k)}))
               %rate1 = [rate1; D1ratenow{sessevpairs{tt}{1}(k)}]; xx1 = [xx1 maa+xbin{sessevpairs{tt}{1}(k)}]; %maa = max(xx1);
               %rate2 = [rate2; D1ratenow{sessevpairs{tt}{2}(k)}]; xx2 = [xx2 mbb+xbin{sessevpairs{tt}{2}(k)}]; %mbb = max(xx2);
               rate1 = D1ratenow{sessevpairs{tt}{1}(k)}; xx1 = xbin{sessevpairs{tt}{1}(k)}; %maa = max(xx1);
               rate2 = D1ratenow{sessevpairs{tt}{2}(k)}; xx2 = xbin{sessevpairs{tt}{2}(k)}; %mbb = max(xx2);
               xjoint1 = xjoint{sessevpairs{tt}{1}(k)}; xjoint2 = xjoint{sessevpairs{tt}{2}(k)}; 
            end
            mmrate1(k) = pinfo.firing.evtmeanrate{i}(sessevpairs{tt}{1}(k)); mmrate2(k) = pinfo.firing.evtmeanrate{i}(sessevpairs{tt}{2}(k));
        end
        pinfo.field.run1DEvtPairMaxRate{i}(tt) = max([max(mmrate1) max(mmrate2)]); pinfo.field.run1DEvtPairMinRate{i}(tt) = min([min(mmrate1) min(mmrate2)]);
        if (~isempty(rate1)) && (~isempty(rate2))
           if (max(rate1) + max(rate2))>0
              [rate1, rate2, x] = realignjointxx(xjoint1, xjoint2, rate1, rate2, xx1, xx2, 'noflip', sessevpairnames); %%%%this revision is to align two trajs across sessions at joints instaead start
              [crr, pp] = correlateratemaps(rate1, rate2, x, x, 0, 0);
              pinfo.field.run1DStablty{i}(tt) = crr; pinfo.field.run1DStabP{i}(tt) = pp; 
              scrr = zeros(1, nshuffle); Z = NaN; 
              for jjj = 1:nshuffle
                   scrr(jjj) = correlateratemaps(rate1, slideshuffle(rate2), x, x, 0, 0);
              end
              if std(scrr)~=0 Z = (crr-mean(scrr))/std(scrr); end
              pinfo.field.run1DSlideShufStablty{i}(tt) = scrr(1); pinfo.field.run1DSlideShufStabltyZ{i}(tt) = Z; 
              pinfo.field.run1DEvtPairRateDiffI{i}(tt) = abs(max(rate1)-max(rate2))/abs(max(rate1)+max(rate2));
           end
        end
    end
   end
  
    %%%%compute 1D place fields
    nfield = 0; baserate = NaN*ones(numel(evTime),1);
    for (tt = 1:numel(evTime)) 
        nx = numel(xbin{tt}); nevtfield = 0;
        if (nx > 1)
        baserate(tt) = computebaserate(D1ratenow{tt}, base1DrateN);
        pf = find1Dfieldprop(D1ratenow{tt}, xbin{tt}, threshold1Drate, max1Dgap, min1Dpeakrate, baserate(tt));
        for (j = 1:numel(pf))
             nfield = nfield + 1; nevtfield = nevtfield + 1;
             pinfo.field.PF1DBoundStart{i}(nfield) = pf(j).boundstart; pinfo.field.PF1DBoundEnd{i}(nfield) = pf(j).boundend;
             pinfo.field.PF1DLength{i}(nfield) = pf(j).leng; pinfo.field.PF1DSize{i}(nfield) = pf(j).size; 
             pinfo.field.PF1DInMeanrate{i}(nfield) = pf(j).inmeanrate; pinfo.field.PF1DInPeakrate{i}(nfield) = pf(j).inpeakrate; 
             pinfo.field.PF1DSkewness{i}(nfield) = pf(j).skew; pinfo.field.PF1DKurtosis{i}(nfield) = pf(j).kurtosis;
             pinfo.field.PF1DLocPeakX{i}(nfield) = pf(j).peakX; pinfo.field.PF1DLocComX{i}(nfield) = pf(j).comX;
             pinfo.field.PF1DLocStartX{i}(nfield) = pf(j).startX; pinfo.field.PF1DLocEndX{i}(nfield) = pf(j).endX;
             pinfo.field.PF1Devt{i}{nfield} = evName{tt}; 
             pinfo.field.PF1DTrajMeanRate{i}(nfield) = pinfo.firing.evtmeanrate{i}(tt);         
        end
        end
        pinfo.field.PF1DNfield{i}(tt) = nevtfield; pinfo.field.PFevtMeanRate{i}(tt) = pinfo.firing.evtmeanrate{i}(tt);
        pinfo.field.PF1DBaseRate{i}(tt) = baserate(tt);
    end
    disp(['-----------> number of 1D fields identified: ', num2str(nfield)]);
end
%%%%%%%%%%%%%%%END OF THE MAIN ROUTINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pf = find2Dfieldprop(rateingrid, xybin, zerothreshold, maxgap, minpeakrate, baselinerate) %%%%do 1D first, 2D is harder to delinearate
pf = [];  %%%%%ATTENTION: rateingrid(y,x)
%%algorithm: similar to the 1D case, subtract from baseline firing rate (for ctx), then examine maximum peaks one by one.
%%           For each peak, accumulate the pixels surround the peak
%%input: rateingrid - spike firing rate in each grid on lineared track grid
%%       zerothreshold - a threshold (in percentage of average firing rate)
%%                    below which the area is not in its place field
%%       maxgap - maximum gap (rate below zero threshold) in a number of pixels within a place field
%%output: numfield - number of place fields
numfield = 0; binsize = xybin{1}(2)-xybin{1}(1); maxgap = round(maxgap/binsize); %%maxgap now in number of bins
gridnum = size(rateingrid); %%%%this is the grid
%%%%subtract from the baseline rate
rateingrid = rateingrid - baselinerate; iii = find(rateingrid<0); if (~isempty(iii)) rateingrid(iii) = zeros(1, numel(iii)); end
%%%%find the peak: rate in peakratenow, loc = [peakind1 peakind2]
[peakR, peakind1] = max(rateingrid); [peakratenow, peakind2] = max(peakR); peakind1 = peakind1(peakind2); 
while (peakratenow >= minpeakrate) %%%%find all the non-zero rate pixels surrounding the peak
    numfield = numfield + 1;
    ind = []; ind(1,1) = peakind1; ind(1,2) = peakind2; np = 1; %%this store all the pixels in the field
    for (i = 1:gridnum(1))
        for (j = 1:gridnum(2))
            if (rateingrid(i,j)>=peakratenow*zerothreshold) %%if the grid rate is high (rates in grids that belong to an already identified field are asigned to 0, see below
                np = np + 1; ind(np,:) = [i j];
            end
        end
    end
    %%%%%now need to check the connectivity (xgap+ygap <= 2*maxgap) in ind
    isneibor = zeros(np,1); isneibor(1) = 1; nzero = find(isneibor == 0); %%%%decide if a grid is connected to the peak ind(1,1)
    for (i = 1:numel(nzero))
        gridtocheck = nzero(i); connectedgrid = ind(isneibor == 1,:);   
        if (min(abs(connectedgrid(:,1)-ind(gridtocheck,1)) + abs(connectedgrid(:,2)-ind(gridtocheck,2))) <= 2*maxgap)
                                                                        %%%this is the key gap check for connectivity
            isneibor(gridtocheck) = 1;
        end
    end
    nzeronow = find(isneibor == 0); 
    while (~isempty(nzeronow)) && (numel(nzeronow) < numel(nzero)) %%%if still room to check
        nzero = nzeronow;
        for (i = 1:numel(nzero))
            gridtocheck = nzero(i); connectedgrid = ind(isneibor == 1,:);   
            if (min(abs(connectedgrid(:,1)-ind(gridtocheck,1)) + abs(connectedgrid(:,2)-ind(gridtocheck,2))) <= 2*maxgap)
               isneibor(gridtocheck) = 1;
            end
        end
        nzeronow = find(isneibor == 0);
    end
    ind = ind(isneibor==1,:); np = size(ind,1); ratenow = zeros(np,1); %%%%this stores all the grids in the field
    for (i = 1:np)
        ratenow(i) = rateingrid(ind(i,1), ind(i,2)); rateingrid(ind(i,1), ind(i,2)) = 0; %%zero rates for the grids in the current field
    end
    %%%%%%characteristics of the current field
    pf(numfield).peakX = xybin{1}(peakind2); pf(numfield).peakY = xybin{2}(peakind1); pf(numfield).inpeakrate = peakratenow + baselinerate;
    pf(numfield).area = np*binsize*binsize; pf(numfield).size = binsize*binsize*sum(ratenow); 
    pf(numfield).inmeanrate = baselinerate + mean(ratenow);
    %%%%%%leave leftovers for computinng the next field
    [peakR, peakind1] = max(rateingrid); [peakratenow, peakind2] = max(peakR); peakind1 = peakind1(peakind2); 
end

function pf = find1Dfieldprop(rateingrid, xbin, zerothreshold, maxgap, minpeakrate, baselinerate)
pf = [];
%%algorithm: subtract from baseline firing rate (for ctx), then examine maximum peak layer by layer 
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
%%%%re-align the 1D grids
xx2 = xx2 - (max(xx2)+min(xx2))/2; xx1 = xx1 - (min(xx1)+max(xx1))/2; xx2 = -xx2;
minx = max([min(xx1) min(xx2)]); maxx = min([max(xx1) max(xx2)]);
gg = minx:binsize:maxx; ind1 = zeros(1, numel(gg)); ind2 = zeros(1, numel(gg)); %%%%re-grid the space
    for (j = 1:numel(gg))
        [mmm, kkk] = min(abs(xx1-gg(j))); ind1(j) = kkk;
        [mmm, kkk] = min(abs(xx2-gg(j))); ind2(j) = kkk;
    end
rate1 = rate1(ind1); rate2 = rate2(ind2);
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
crr = 0; pp = NaN;
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

function rate = slideshuffle(rate)
[mm,nn] = size(rate); 
if mm > 1
   rate = rate((1 + mod((1:mm)+round(rand*mm), mm)), :); %%slide shuffling
end
if nn > 1
   rate = rate(:, (1 + mod((1:nn)+round(rand*nn), nn)) ); %%slide shuffling
end
function rate = D2slideshuffle(rate, occup)
%%%%2D slide shuffle is hard: this is psudo slide shuffling
[mm,nn] = size(rate); r = reshape(rate, mm*nn, 1); o = reshape(occup, mm*nn, 1);
iii = find(o>0); R = r(iii); nnow = numel(R);
R = R((1 + mod((1:nnow)+round(rand*nnow), nnow)), :); %%slide shuffling
r(iii)= R; rate =reshape(r, mm, nn);

function [dirpairs, sessevpairs, dirpairnames, sessevpairnames] = pairupevts(evName, evSess, evType)
%%%dirpairs{npair}[1,2], dirpairnames{npair}{1,2} = traj pairs' indices and names for computing directionality 
%%%sessevpairs{npair}{1,2}, sessevpairnames{npair}{1,2} = two groups of events' indices and names for computing stability 
dirpairs = []; sessevpairs = []; dirpairnames = []; sessevpairnames = []; %disp(evSess); disp(evSess(find(strcmp(evType, 'run'))));
ndirpair = 0; nevpair = 0; %disp(evSess); disp(evType); 
aa = evSess(find(strcmp(evType, 'run'))); ii = ones(1, numel(aa)); 
for i=1:numel(aa) if isempty(aa{i}) ii(i) = 0; end; end; aa = aa(find(ii));
allsess = unique(aa);
%allsess = unique(evSess(find(strcmp(evType, 'run')))); 
for (i = 1:numel(allsess)) %%%the following assumes that the two events in a session are the same trajectory in opposite directions - completely wrong in many experiments
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
            %nevpair = nevpair + 1; sessevpairs{nevpair} = pairnow; sessevpairnames{nevpair} = pairnames; %%%%this is for double trajectory stablity
            for (k = 1:numel(pairnow))
                nevpair = nevpair + 1; sessevpairs{nevpair} = pairnow{k}; sessevpairnames{nevpair} = pairnames{k}; %%%%this is for single trajectory stablity
            end
        end
    end
end

function [pairnow, pairnames, ok] = ifmatches(sessEvName1, evind1, sessEvName2, evind2) 
pairnow = [] ; pairnames = []; pairnames = []; nelm = 0; ok = 1;
%if (numel(evind1) ~= numel(evind2))
%    ok = 0;
%else
    [all2ndstr, all2ndtraj] = strtok(sessEvName2, '_');
    for (i = 1:numel(evind1))
        [str, trajnow] = strtok(sessEvName1{i}, '_'); 
        jj = find(strcmp(lower(all2ndtraj), lower(trajnow)));
        if (numel(jj) == 1)
            %nelm = nelm + 1; pairnow{nelm}{1}(1) = evind1(i); pairnow{nelm}{2}(1) = evind2(jj);
            %pairnames{nelm}{1} = sessEvName1{i}; pairnames{nelm}{2} = sessEvName2{jj};
            nelm = nelm + 1; pairnow{nelm}{1} = evind1(i); pairnow{nelm}{2} = evind2(jj);
            pairnames{nelm} = strcat(sessEvName1{i}, '-', sessEvName2{jj});
        end
    end
%end
if (nelm==0) ok = 0; end

% function [pairnow, pairnames, ok] = ifmatches(sessEvName1, evind1, sessEvName2, evind2) %%%this is for double trajectory stability
% pairnow{1} = [] ; pairnow{2} = []; pairnames{1} = []; pairnames{2} = []; nelm = 0; ok = 1;
% %if (numel(evind1) ~= numel(evind2))
% %    ok = 0;
% %else
%     [all2ndstr, all2ndtraj] = strtok(sessEvName2, '_');
%     for (i = 1:numel(evind1))
%         [str, trajnow] = strtok(sessEvName1{i}, '_'); 
%         jj = find(strcmp(all2ndtraj, trajnow));
%         if (numel(jj) == 1)
%             nelm = nelm + 1; pairnow{1}(nelm) = evind1(i); pairnow{2}(nelm) = evind2(jj);
%             pairnames{1} = strcat(pairnames{1},'+', sessEvName1{i}); pairnames{2} = strcat(pairnames{2},'+', sessEvName2{jj});
%         end
%     end
% %end
% if (isempty(pairnow{1})) ok = 0; end
    
function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end

function [D2rate, nact] = get2Dratemaps(spiketime, XYbin, gridoccup, postime, XX, YY, smParm, frametime)
nx = numel(XYbin{1}); ny = numel(XYbin{2}); D2rate = NaN*ones(ny,nx); nact = NaN;
frametime = 1.03*frametime; %%%slightly expand frametime to assure detection, due to slightly variable frame-to-frame time
if (nx > 1) && (ny > 1)
    binsize = XYbin{1}(2) - XYbin{1}(1); sigma = round(smParm.d2sigma/binsize); %%%sigma now in number of bins
    spiketime = sort(spiketime); [postime, iii] = sort(postime); XX = XX(iii); YY = YY(iii); 
    lastpoint = 1;  %this only for saving time
    x = NaN*ones(size(spiketime)); y = NaN*ones(size(spiketime));
    for (j = 1:numel(spiketime))
        for (k = lastpoint:numel(postime)-1) %find corresponding time in position data
            if abs(postime(k) - spiketime(j)) <= frametime 
            %if (postime(k) <= spiketime(j)) && (postime(k+1) > spiketime(j)) 
                x(j) = XX(k); y(j) = YY(k); lastpoint = k; break; 
            end
        end
    end
    %%%%%count spikes
   countx = cell(1, nx); county = cell(1, ny); count = zeros(ny, nx); %%%this asignment is for convenient plot using pcolor(x-columnx, y-rows)
   for (i = 1:nx) countx{i} = find((x>=XYbin{1}(i)) & (x<XYbin{1}(i)+binsize)); end
   for (i = 1:ny) county{i} = find((y>=XYbin{2}(i)) & (y<XYbin{2}(i)+binsize)); end
   %countx{nx} = find(x>=XYbin{1}(nx)); county{ny} = find(y>=XYbin{2}(ny));
   %%%%% now spike count at grid (i,j) = number of common indices in countx{i} and county{j}
   for (i = 1:nx)
   for (j = 1:ny)
       common = []; count(j,i) = 0;
       if (~isempty(countx{i}) & ~isempty(county{j}))
          common = intersect(countx{i}, county{j});
          if (~isempty(common)) count(j,i)=numel(common); end  %%assign (j,i) because pcolor requires 
       end
   end
   end
   %%%%%compute firing rate
   for (j = 1:ny)
   for (i = 1:nx)
       if (gridoccup(j,i) > 0) D2rate(j,i) = count(j,i)/gridoccup(j,i); end
   end
   end
   nact = numel(find(count>0))/numel(find(gridoccup>0));
   %%now do a 2d smoothing
   if (strcmpi(smParm.d2sm, 'yes'))
      D2rate = TwoDSmooth_new(D2rate, gridoccup, sigma, sigma, smParm.d2Nsig, smParm.d2Nsig);
   end
end

function [D1rate, nact] = getlinearmap(spiketime, evTime, lappostime, lapx, xbin, occutime, smParm, frametime)
nx = numel(xbin); D1rate = zeros(nx,1); %D1rate = NaN*ones(nx,1); 
nlap = numel(lappostime); count = zeros(nx,1); spikex = []; nact = NaN;
frametime = 1.03*frametime; %%%slightly expand frametime to assure detection, due to slightly variable frame-to-frame time
if (nx > 1)
    binsize = xbin(2) - xbin(1); sigma = round(smParm.d1sigma/binsize); %%sigma now in number of bins
%%%%%count spikes in xbin lap bu lap
for (i = 1:nlap)
    counti = zeros(nx,1);
    spikenow = sort( spiketime( (spiketime>=evTime.start(i)) & (spiketime<=evTime.ent(i)) ) );
    spikex = NaN*ones(size(spikenow)); [laptimenow, iii] = sort(lappostime{i}); allxnow = lapx{i}(iii);
    lastpoint = 1;  %this only for saving time
    for (j = 1:numel(spikenow))
         for (k = lastpoint:numel(laptimenow)) %find corresponding time in position data
             if abs(laptimenow(k) - spikenow(j)) <= frametime 
             %if (laptimenow(k) <= spikenow(j)) && (laptimenow(k+1) > spikenow(j)) 
                 spikex(j) = allxnow(k); lastpoint = k; break; 
             end
         end
    end
    for (j = 1:nx) counti(j) = numel(find((spikex>=xbin(j)) & (spikex<xbin(j)+binsize))); end
    %counti(nx) = numel(find(spikex>=xbin(nx)));
    count = count + counti;
end
%%%%%compute firing rate
for (i = 1:nx)
     if (occutime(i) ~= 0) D1rate(i) = count(i)/occutime(i); end
end
nact = numel(find(count>0))/numel(find(occutime>0));
%%now do a 1d smoothing
if (strcmpi(smParm.d1sm, 'yes'))
    otime = ones(size(occutime));
    D1rate = OneDSmooth_new(D1rate, otime, sigma, smParm.d1Nsig); %%%smooth first: smoothing does not remove NaN's (occutime = 0)
    %%%need to do a 1D interpolation of the NaN data points
    X = (1:numel(D1rate))'; iii = find(D1rate > -10); 
    if numel(iii)>1
       R = D1rate(iii); Y=X(iii);
       D1rate = interpn(Y, R, X);
       %%%%The following has the problem of reducing rate at the end of traj
       %XX = [-smParm.d1Nsig*sigma:1:smParm.d1Nsig*sigma]; 
       %wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
       %D1rate = conv(D1rate, wind); D1rate = D1rate((nw-1)/2+1:numel(D1rate)-(nw-1)/2);
    end
end
end

function [spinfo, spinfotime, sparsty] = findratemapprop(rate, occu)
%%%%%%another variable: mutual information is not computed here
spinfo = NaN; spinfotime = NaN; sparsty = NaN;  
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

function [r1, r2, x] = realignjointxx(xjoint1, xjoint2, rate1, rate2, xx1, xx2, flipflag, evtflag)
%%%%realign trajectories assuming the center landmarks are matched
r1 = []; x = []; r2 = []; rate1 = rate1'; rate2 = rate2'; binsize = mean(diff(xx1));
if numel(xjoint1) == numel(xjoint2)
   nj = numel(xjoint1); 
   if strncmpi(flipflag, 'yes' ,1) %%%%assuming rate1/rate2 xx1/xx2 are opposite trajectories
       xjoint2 = max(xjoint2) - flip(xjoint2); rate2 = flip(rate2); 
       xx2 = max(xjoint2)-flip(xx2); %%%need to be very careful here, some xx not start at 0 or end at max(xjoint)
   end
   if (nj == 2)
       if (numel(rate1) == numel(rate2))
           r1 = rate1; r2 = rate2; x = xx1;
       elseif (numel(rate1) < numel(rate2))
           r1 = interpn(xx1, rate1, xx2); r2 = rate2; x = xx2; 
       else
           r2 = interpn(xx2, rate2, xx1); r1 = rate1; x = xx1;
       end
   elseif (nj>=3)&&(nj<=4)
     xjoint = zeros(size(xjoint1));
     for (i = 1:nj) xjoint(i) = mean([xjoint1(i) xjoint2(i)]); end
     for (i = 2:numel(xjoint1)) %%
        jj1 = find( (xx1>=xjoint1(i-1)) & (xx1<xjoint1(i))); jj2 = find( (xx2>=xjoint2(i-1)) & (xx2<xjoint2(i))); 
        if (numel(jj1) == numel(jj2))
            r1 = [r1 rate1(jj1)]; r2 = [r2 rate2(jj2)]; xnow = binsize*(1:numel(jj1)); x = [x xjoint(i-1)+xnow]; 
        else
            xnow = binsize:binsize:(xjoint(i)-xjoint(i-1));
            r1now = interpn(xx1(jj1)-xjoint1(i-1), rate1(jj1), xnow); 
            r2now = interpn(xx2(jj2)-xjoint2(i-1), rate2(jj2), xnow);
            r1 = [r1 r1now]; r2 = [r2 r2now]; x = [x xjoint(i-1)+xnow];
        end
    end
   else
     xjoint = zeros(size(xjoint1));
     for (i = 1:nj) xjoint(i) = mean([xjoint1(i) xjoint2(i)]); end
     for (i = 2:numel(xjoint1)) %%
     %for (i = 3:numel(xjoint1)-1) %%%get rid the smoothing effect at the edge
        jj1 = find( (xx1>=xjoint1(i-1)) & (xx1<xjoint1(i))); jj2 = find( (xx2>=xjoint2(i-1)) & (xx2<xjoint2(i))); 
        if (numel(jj1) == numel(jj2))
            r1 = [r1 rate1(jj1)]; r2 = [r2 rate2(jj2)]; xnow = binsize*(1:numel(jj1)); x = [x xjoint(i-1)+xnow]; 
        else
            xnow = binsize:binsize:(xjoint(i)-xjoint(i-1));
            if (numel(jj1)>=2) && (numel(jj2)>=2)
               r1now = interpn(xx1(jj1)-xjoint1(i-1), rate1(jj1), xnow); 
               r2now = interpn(xx2(jj2)-xjoint2(i-1), rate2(jj2), xnow);
            else
                jjk = min([numel(jj1) numel(jj2)]); if jjk== 1
                r1now = mean(rate1(jj1)); r2now = mean(rate2(jj2)); 
                else
                    r1now = []; r2now = [];
                end
            end
            r1 = [r1 r1now]; r2 = [r2 r2now]; x = [x xjoint(i-1)+xnow];
        end
    end
   end
else
   disp(['-----------> warning: trajecotry joint points do not match: ', evtflag]); 
end
% figure; ha = axes('NextPlot', 'add'); plot(x, r1, 'Parent', ha, 'Color', [0 0 1]); plot(x, r2, 'Parent', ha, 'Color', [1 0 0]);
% for (i = 1:nj)
%     plot([xjoint(i) xjoint(i)], [0 10], 'Parent', ha, 'Color', [0.5 0.5 0.5]);
% end

function [crossR, crossX, crossP] = Dircrosscrr(rate1, rate2, binsize, maxlag)
%figure; plot(rate1); figure; plot(rate2);
crossR = NaN; crossX = NaN; crossP = NaN; maxbin = round(maxlag/binsize);
ii = find( (~isnan(rate1)) & (~isnan(rate2)) ); rate1 = rate1(ii); rate2 = rate2(ii);
[ccc, lags] = Utilities_FindXCorrCoef_self(rate1,rate2,maxbin);  %%input row vectors
[crossR, i] = max(ccc); crossX = binsize*lags(i); np = numel(rate1);
k = lags(i); %%%this is the shifted time point with max cc
%disp([numel(rate1) numel(rate2) crossX k np maxbin binsize])
if (np-abs(k))>5 %if shifting is not too large (close to edge)
  if (k>=0) 
   [C,P] = corrcoef(rate1(k+1:np), rate2(1:np-k)); crossP = P(1,2);
  else
   [C,P] = corrcoef(rate1(1:np+k), rate2(1-k:np)); crossP = P(1,2);
  end
end





