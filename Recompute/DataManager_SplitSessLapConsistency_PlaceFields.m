function DataManager_SplitSessLapConsistency_PlaceFields
%Compute the dynamics of spatial firing profile: - split laps in a session into
%subgroups (early/late) and then compute dynamics as in lap-by-lap dynamics

%For (open field) session: original PFDynam can be used to do the same;
%just need to set the parameter behav.parm.sessSegTime

hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); %tagmark = get(gcbo, 'Tag');
hgroup = getappdata(hf, 'hgroup'); %hfield = getappdata(hf, 'hfield');
ok = 1; oc = 0; plotparm = getappdata(hf, 'plotparm'); 
disp('Split track sessions to groups of laps and re-compute spatial variables...');
%get selected cellind
groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
if isempty(cellind) 
    ok = 0; disp('-----> no cells selected; aborted'); 
end 
if ok
    if plotparm.linkbehav == 1
       behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
    else
       disp('-----> no behaviroal database linked; aborted'); ok = 0;
    end
end
if ok
    input = inputdlg({'Split by time/lap number'; 'Duration/lap number'; 'Min time/lap number'}, 'Split method', 3, {'lap'; '5'; '5'}); 
    if (~isempty(input))
        parm.splitmethod = input{1}; parm.grouplength = str2num(input{2}); parm.grpminlapnum = str2num(input{3}); 
    else
        ok = 0;
    end
end 
if ok
   ow = plotparm.overwrite; %if ow=1, plot into the current database
   [writefilename, ok] = getoutputfile(hf, ow);
end
if ok
    if isfield(pinfo, 'fielddynam')
        [pinfo, data] = splitsession_fielddynamvars(pinfo, data, behav, bhdata, cellind, parm); oc = 1;
    else
        disp('-----> warning: field dynamics variables not computed; ');
    end
    if isfield(pinfo, 'PVcrr')
       [pinfo, data] = splitsession_PVcrr(pinfo, data, behav, bhdata, cellind, parm); oc = 1;
    else
       disp('-----> warning: PV correlation variables not computed');
    end
end
if ok
    if oc
        save(writefilename, 'pinfo', 'data', '-mat', '-v7.3');
        if (ow)
            iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
            fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
            newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
        else
            hmain = DataManager_DataManager_Callback; plotparm.linkspike = 1; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
        end
        setappdata(hmain, 'plotparm', plotparm);
        DataManager_PlotSpikeDatabase(hmain, pinfo, data, 'Cell', 'clname', '.spikedb');
        set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
        setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); pinfo = []; data = [];
    else
        disp('-----> database not altered; no new data to save');
    end 
end
disp('****************');

function [pinfo, data] = splitsession_fielddynamvars(pinfo, data, behav, bhdata, cellind, parm)
%%%%%% Spatial field dynamics
if (~isempty(cellind))
   nspike = numel(pinfo.general.parmfile);
   %%%%assign parm variables
   if (~isfield(pinfo.parm, 'SplitMethod')) pinfo.parm.SplitMethod = cell(1, nspike); end   
   if (~isfield(pinfo.parm, 'lapgrpLength')) pinfo.parm.lapgrpLength = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'lapgrpMinNumber')) pinfo.parm.lapgrpMinNumber = cell(1, nspike); end
   for (ttt = 1:numel(cellind))
       pinfo.parm.SplitMethod{cellind(ttt)} = parm.splitmethod;
       pinfo.parm.lapgrpLength{cellind(ttt)} = parm.grouplength; 
       pinfo.parm.lapgrpMinNumber{cellind(ttt)} = parm.grpminlapnum;
   end
   %%% define Spatial correlation
   if (~isfield(pinfo.fielddynam, 'lapgrpPairName')) pinfo.fielddynam.lapgrpPairName = cell(1, nspike); end   
   if (~isfield(pinfo.fielddynam, 'lapgrpPairMeanrate1')) pinfo.fielddynam.lapgrpPairMeanrate1 = cell(1, nspike); end  
   if (~isfield(pinfo.fielddynam, 'lapgrpPairMeanrate2')) pinfo.fielddynam.lapgrpPairMeanrate2 = cell(1, nspike); end  
   if (~isfield(pinfo.fielddynam, 'lapgrpCrr')) pinfo.fielddynam.lapgrpCrr = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpCrrP')) pinfo.fielddynam.lapgrpCrrP = cell(1, nspike); end
   %%% Rate curve parameters
   if (~isfield(pinfo.fielddynam, 'lapgrpName')) pinfo.fielddynam.lapgrpName = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpMeanrate')) pinfo.fielddynam.lapgrpMeanrate = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpPeakrate')) pinfo.fielddynam.lapgrpPeakrate = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpPeakX')) pinfo.fielddynam.lapgrpPeakX = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpSpatialInfo')) pinfo.fielddynam.lapgrpSpatialInfo = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpSpatialInfotime')) pinfo.fielddynam.lapgrpSpatialInfotime = cell(1, nspike); end
   %%% define dominant field parameters
   if (~isfield(pinfo.fielddynam, 'lapgrpfMeanRate')) pinfo.fielddynam.lapgrpfMeanRate = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpfPeakRate')) pinfo.fielddynam.lapgrpfPeakRate = cell(1, nspike); end   
   if (~isfield(pinfo.fielddynam, 'lapgrpfComX')) pinfo.fielddynam.lapgrpfComX = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpfPeakX')) pinfo.fielddynam.lapgrpfPeakX = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpfFieldSize')) pinfo.fielddynam.lapgrpfFieldSize = cell(1, nspike); end 
   if (~isfield(pinfo.fielddynam, 'grpavgFstart')) pinfo.fielddynam.grpavgFstart = cell(1, nspike); end 
   if (~isfield(pinfo.fielddynam, 'grpavgFend')) pinfo.fielddynam.grpavgFend = cell(1, nspike); end 
   %%%%%data output
   data.splitgrp.grpevnames = cell(1, nspike); data.splitgrp.grplapind = cell(1, nspike); data.splitgrp.grplapOccutime = cell(1, nspike);
   data.splitgrp.grpratemaps = cell(1, nspike); data.splitgrp.grpOccutime = cell(1, nspike);
end
for (jjjk = 1:numel(cellind))
    i = cellind(jjjk); 
    disp(strcat('-----> compute (split lap) field dynamics ---', pinfo.general.parmfile{i}));
    fprate = pinfo.field.PF1DInPeakrate{i}; 
    if (~isempty(fprate)) %%%Only do this if fields exist on linear trajectories
        smParm.d1sigma = pinfo.parm.fSmooth1DSigma(i); smParm.d1Nsig = pinfo.parm.fSmooth1DNSigma(i); smParm.d1sm = pinfo.parm.fSmooth1D{i};
        lapevName = pinfo.fielddynam.run1DEvtname{i}; withinlapNum = pinfo.fielddynam.runWithinLapNum{i}; 
        xxbin = data.fielddynam.evt1Dxbin{i}; D1ratenow = data.fielddynam.evt1DRateMapsAligned{i};
        [grpname, grplapind, grplapOccutime] = resolvegrplaps(pinfo, data, behav, bhdata, lapevName, withinlapNum, parm, i, xxbin); 
        ng = numel(grpname); grpratemaps = cell(1,ng); grpOccutime = cell(1,ng); 
        for ii=1:ng
            [grpratemaps{ii}, grpOccutime{ii}] = findavgratemaps(D1ratenow(grplapind{ii}), grplapOccutime{ii}, xxbin);
            grpratemaps{ii} = smoothratecurvenow(grpratemaps{ii}, grpOccutime{ii}, xxbin(2)-xxbin(1), smParm);
        end
        [avgmap, occutime] = findavgratemaps(grpratemaps, grpOccutime, xxbin); 
        avgmap = smoothratecurvenow(avgmap, occutime, xxbin(2)-xxbin(1), smParm); %hf = figure; plot(xxbin, avgmap);
%%%%%%%%%%%%assignd data output
        data.splitgrp.grpevnames{i} = grpname; data.splitgrp.grplapind{i} = grplapind; data.splitgrp.grplapOccutime{i} = grplapOccutime; 
        data.splitgrp.grpratemaps{i} = grpratemaps; data.splitgrp.grpOccutime{i} = grpOccutime; 
%%%%%%%%%%%%%%%%ratemap spatial variables and spatial correlations across groups
        pinfo.fielddynam.lapgrpName{i} = cell(1, ng); pinfo.fielddynam.lapgrpMeanrate{i} = NaN*ones(1, ng);
        pinfo.fielddynam.lapgrpPeakrate{i} = NaN*ones(1, ng); pinfo.fielddynam.lapgrpPeakX{i} = NaN*ones(1, ng);
        pinfo.fielddynam.lapgrpSpatialInfo{i} = NaN*ones(1, ng); pinfo.fielddynam.lapgrpSpatialInfotime{i} = NaN*ones(1, ng);         
        pinfo.fielddynam.lapgrpPairName{i} = cell(1, ng*(ng-1)/2); 
        pinfo.fielddynam.lapgrpCrr{i} = cell(1, ng*(ng-1)/2); pinfo.fielddynam.lapgrpCrrP{i} = cell(1, ng*(ng-1)/2); 
        nc = 0;
        for (m=1:ng)  
            name1 = strcat(grpname{m}, '_', num2str(2-mod(m,2)));
            [spinfo, spinfotime, ~, mmeanrate1] = findratemapprop(grpratemaps{m}, grpOccutime{m});
            [vv, II] = max(grpratemaps{m}); pinfo.fielddynam.lapgrpPeakrate{i}(m) = vv; pinfo.fielddynam.lapgrpPeakX{i}(m) = xxbin(II); 
            pinfo.fielddynam.lapgrpName{i}{m} = name1; pinfo.fielddynam.lapgrpMeanrate{i}(m) = mmeanrate1;
            pinfo.fielddynam.lapgrpSpatialInfo{i}(m) = spinfo; pinfo.fielddynam.lapgrpSpatialInfotime{i}(m) = spinfotime;
            for(n=m+1:ng)
                name2 = strcat(grpname{n}, '_', num2str(2-mod(n,2))); nc = nc + 1;
                pinfo.fielddynam.lapgrpPairName{i}{nc} = strcat(name1, '-', name2);
                [~, ~, ~, mmeanrate2] = findratemapprop(grpratemaps{n}, grpOccutime{n});
                pinfo.fielddynam.lapgrpPairMeanrate1{i}{nc} = mmeanrate1; pinfo.fielddynam.lapgrpPairMeanrate2{i}{nc} = mmeanrate2; 
                [crr, pp0] = correlateratemaps(grpratemaps{m}, grpratemaps{n}, xxbin, xxbin, 0, 0);
                pinfo.fielddynam.lapgrpCrr{i}{nc} = crr; pinfo.fielddynam.lapgrpCrrP{i}{nc} = pp0; 
            end
        end
%%%%%%%%%%%%%%place field properties
        pinfo.fielddynam.lapgrpfMeanRate{i} = NaN*ones(1, ng); pinfo.fielddynam.lapgrpfPeakRate{i} = NaN*ones(1, ng);
        pinfo.fielddynam.lapgrpfComX{i} = NaN*ones(1, ng); pinfo.fielddynam.lapgrpfPeakX{i} = NaN*ones(1, ng);        
        pinfo.fielddynam.lapgrpfFieldSize{i} = NaN*ones(1, ng);
        %%%%%%%%%%%%field computing parameters     
        threshold1Drate = pinfo.parm.f1DThresholdRate(i); max1Dgap = pinfo.parm.f1DMaxGap(i); 
        min1Dpeakrate = pinfo.parm.f1DMinPeakRate(i)/3; %%%% drop this to make sure the original peak field is recovered 
        base1DrateN = pinfo.parm.f1DBaseRateGridN(i); 
   %%%%%%%%%%redefine fields on overall averaged map: to get field start/end
        baserate = computebaserate(avgmap, base1DrateN); 
        pf = find1Dfieldprop(avgmap, xxbin, threshold1Drate, max1Dgap, min1Dpeakrate, baserate);
        if (numel(pf)>=1)
           disp('------------> group-averaged place field found');
           ppeakrate = []; for (j = 1:numel(pf)) ppeakrate(j) = pf(j).inpeakrate; end
           [~, ii] = max(ppeakrate); fstart = pf(ii).boundstart; fend = pf(ii).boundend;
   %%%%%%%%%%%compute group field variables
           for (m = 1:ng)
               [mmrate, pprate, fsize, comX, ploc, ~, ~] = find1DFFprop(grpratemaps{m}, xxbin, fstart, fend);
               pinfo.fielddynam.lapgrpfMeanRate{i}(m) = mmrate; pinfo.fielddynam.lapgrpfPeakRate{i}(m) = pprate;
               pinfo.fielddynam.lapgrpfComX{i}(m) = comX; pinfo.fielddynam.lapgrpfPeakX{i}(m) = ploc; 
               pinfo.fielddynam.lapgrpfFieldSize{i}(m) = fsize;
           end
           pinfo.fielddynam.grpavgFstart{i} = fstart; pinfo.fielddynam.grpavgFend{i} = fend;
        else
            disp('------------> no group-averaged 1D field found');
        end
    else
        disp('------------> lap-group fields not computed: no place fields found');
    end
end


function [pinfo, data] = splitsession_PVcrr(pinfo, data, behav, bhdata, cellind, parm)
%%%%%% Spatial field dynamics
if (~isempty(cellind))
   nspike = numel(pinfo.general.parmfile);
   %%%%assign parm variables
   if (~isfield(pinfo.parm, 'SplitMethod')) pinfo.parm.SplitMethod = cell(1, nspike); end   
   if (~isfield(pinfo.parm, 'lapgrpLength')) pinfo.parm.lapgrpLength = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'lapgrpMinNumber')) pinfo.parm.lapgrpMinNumber = cell(1, nspike); end
   for (ttt = 1:numel(cellind))
       pinfo.parm.SplitMethod{cellind(ttt)} = parm.splitmethod;
       pinfo.parm.lapgrpLength{cellind(ttt)} = parm.grouplength; 
       pinfo.parm.lapgrpMinNumber{cellind(ttt)} = parm.grpminlapnum;
   end
   %%% define Spatial correlation
   if (~isfield(pinfo.fielddynam, 'lapgrpPairName')) pinfo.fielddynam.lapgrpPairName = cell(1, nspike); end   
   if (~isfield(pinfo.fielddynam, 'lapgrpPairMeanrate1')) pinfo.fielddynam.lapgrpPairMeanrate1 = cell(1, nspike); end  
   if (~isfield(pinfo.fielddynam, 'lapgrpPairMeanrate2')) pinfo.fielddynam.lapgrpPairMeanrate2 = cell(1, nspike); end  
   if (~isfield(pinfo.fielddynam, 'lapgrpCrr')) pinfo.fielddynam.lapgrpCrr = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpCrrP')) pinfo.fielddynam.lapgrpCrrP = cell(1, nspike); end
   %%% Rate curve parameters
   if (~isfield(pinfo.fielddynam, 'lapgrpName')) pinfo.fielddynam.lapgrpName = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpMeanrate')) pinfo.fielddynam.lapgrpMeanrate = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpPeakrate')) pinfo.fielddynam.lapgrpPeakrate = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpPeakX')) pinfo.fielddynam.lapgrpPeakX = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpSpatialInfo')) pinfo.fielddynam.lapgrpSpatialInfo = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpSpatialInfotime')) pinfo.fielddynam.lapgrpSpatialInfotime = cell(1, nspike); end
   %%% define dominant field parameters
   if (~isfield(pinfo.fielddynam, 'lapgrpfMeanRate')) pinfo.fielddynam.lapgrpfMeanRate = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpfPeakRate')) pinfo.fielddynam.lapgrpfPeakRate = cell(1, nspike); end   
   if (~isfield(pinfo.fielddynam, 'lapgrpfComX')) pinfo.fielddynam.lapgrpfComX = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpfPeakX')) pinfo.fielddynam.lapgrpfPeakX = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'lapgrpfFieldSize')) pinfo.fielddynam.lapgrpfFieldSize = cell(1, nspike); end 
   if (~isfield(pinfo.fielddynam, 'grpavgFstart')) pinfo.fielddynam.grpavgFstart = cell(1, nspike); end 
   if (~isfield(pinfo.fielddynam, 'grpavgFend')) pinfo.fielddynam.grpavgFend = cell(1, nspike); end 
   %%%%%data output
   data.splitgrp.grpevnames = cell(1, nspike); data.splitgrp.grplapind = cell(1, nspike); data.splitgrp.grplapOccutime = cell(1, nspike);
   data.splitgrp.grpratemaps = cell(1, nspike); data.splitgrp.grpOccutime = cell(1, nspike);
end
for (jjjk = 1:numel(cellind))
    i = cellind(jjjk); 
    disp(strcat('-----> compute (split lap) field dynamics ---', pinfo.general.parmfile{i}));
    fprate = pinfo.field.PF1DInPeakrate{i}; 
    if (~isempty(fprate)) %%%Only do this if fields exist on linear trajectories
        smParm.d1sigma = pinfo.parm.fSmooth1DSigma(i); smParm.d1Nsig = pinfo.parm.fSmooth1DNSigma(i); smParm.d1sm = pinfo.parm.fSmooth1D{i};
        lapevName = pinfo.fielddynam.run1DEvtname{i}; withinlapNum = pinfo.fielddynam.runWithinLapNum{i}; 
        xxbin = data.fielddynam.evt1Dxbin{i}; D1ratenow = data.fielddynam.evt1DRateMapsAligned{i};
        [grpname, grplapind, grplapOccutime] = resolvegrplaps(pinfo, data, behav, bhdata, lapevName, withinlapNum, parm, i, xxbin); 
        ng = numel(grpname); grpratemaps = cell(1,ng); grpOccutime = cell(1,ng); 
        for ii=1:ng
            [grpratemaps{ii}, grpOccutime{ii}] = findavgratemaps(D1ratenow(grplapind{ii}), grplapOccutime{ii}, xxbin);
            grpratemaps{ii} = smoothratecurvenow(grpratemaps{ii}, grpOccutime{ii}, xxbin(2)-xxbin(1), smParm);
        end
        [avgmap, occutime] = findavgratemaps(grpratemaps, grpOccutime, xxbin); 
        avgmap = smoothratecurvenow(avgmap, occutime, xxbin(2)-xxbin(1), smParm); %hf = figure; plot(xxbin, avgmap);
%%%%%%%%%%%%assignd data output
        data.splitgrp.grpevnames{i} = grpname; data.splitgrp.grplapind{i} = grplapind; data.splitgrp.grplapOccutime{i} = grplapOccutime; 
        data.splitgrp.grpratemaps{i} = grpratemaps; data.splitgrp.grpOccutime{i} = grpOccutime; 
%%%%%%%%%%%%%%%%ratemap spatial variables and spatial correlations across groups
        pinfo.fielddynam.lapgrpName{i} = cell(1, ng); pinfo.fielddynam.lapgrpMeanrate{i} = NaN*ones(1, ng);
        pinfo.fielddynam.lapgrpPeakrate{i} = NaN*ones(1, ng); pinfo.fielddynam.lapgrpPeakX{i} = NaN*ones(1, ng);
        pinfo.fielddynam.lapgrpSpatialInfo{i} = NaN*ones(1, ng); pinfo.fielddynam.lapgrpSpatialInfotime{i} = NaN*ones(1, ng);         
        pinfo.fielddynam.lapgrpPairName{i} = cell(1, ng*(ng-1)/2); 
        pinfo.fielddynam.lapgrpCrr{i} = cell(1, ng*(ng-1)/2); pinfo.fielddynam.lapgrpCrrP{i} = cell(1, ng*(ng-1)/2); 
        nc = 0;
        for (m=1:ng)  
            name1 = strcat(grpname{m}, '_', num2str(2-mod(m,2)));
            [spinfo, spinfotime, ~, mmeanrate1] = findratemapprop(grpratemaps{m}, grpOccutime{m});
            [vv, II] = max(grpratemaps{m}); pinfo.fielddynam.lapgrpPeakrate{i}(m) = vv; pinfo.fielddynam.lapgrpPeakX{i}(m) = xxbin(II); 
            pinfo.fielddynam.lapgrpName{i}{m} = name1; pinfo.fielddynam.lapgrpMeanrate{i}(m) = mmeanrate1;
            pinfo.fielddynam.lapgrpSpatialInfo{i}(m) = spinfo; pinfo.fielddynam.lapgrpSpatialInfotime{i}(m) = spinfotime;
            for(n=m+1:ng)
                name2 = strcat(grpname{n}, '_', num2str(2-mod(n,2))); nc = nc + 1;
                pinfo.fielddynam.lapgrpPairName{i}{nc} = strcat(name1, '-', name2);
                [~, ~, ~, mmeanrate2] = findratemapprop(grpratemaps{n}, grpOccutime{n});
                pinfo.fielddynam.lapgrpPairMeanrate1{i}{nc} = mmeanrate1; pinfo.fielddynam.lapgrpPairMeanrate2{i}{nc} = mmeanrate2; 
                [crr, pp0] = correlateratemaps(grpratemaps{m}, grpratemaps{n}, xxbin, xxbin, 0, 0);
                pinfo.fielddynam.lapgrpCrr{i}{nc} = crr; pinfo.fielddynam.lapgrpCrrP{i}{nc} = pp0; 
            end
        end
%%%%%%%%%%%%%%place field properties
        pinfo.fielddynam.lapgrpfMeanRate{i} = NaN*ones(1, ng); pinfo.fielddynam.lapgrpfPeakRate{i} = NaN*ones(1, ng);
        pinfo.fielddynam.lapgrpfComX{i} = NaN*ones(1, ng); pinfo.fielddynam.lapgrpfPeakX{i} = NaN*ones(1, ng);        
        pinfo.fielddynam.lapgrpfFieldSize{i} = NaN*ones(1, ng);
        %%%%%%%%%%%%field computing parameters     
        threshold1Drate = pinfo.parm.f1DThresholdRate(i); max1Dgap = pinfo.parm.f1DMaxGap(i); 
        min1Dpeakrate = pinfo.parm.f1DMinPeakRate(i)/3; %%%% drop this to make sure the original peak field is recovered 
        base1DrateN = pinfo.parm.f1DBaseRateGridN(i); 
   %%%%%%%%%%redefine fields on overall averaged map: to get field start/end
        baserate = computebaserate(avgmap, base1DrateN); 
        pf = find1Dfieldprop(avgmap, xxbin, threshold1Drate, max1Dgap, min1Dpeakrate, baserate);
        if (numel(pf)>=1)
           disp('------------> group-averaged place field found');
           ppeakrate = []; for (j = 1:numel(pf)) ppeakrate(j) = pf(j).inpeakrate; end
           [~, ii] = max(ppeakrate); fstart = pf(ii).boundstart; fend = pf(ii).boundend;
   %%%%%%%%%%%compute group field variables
           for (m = 1:ng)
               [mmrate, pprate, fsize, comX, ploc, ~, ~] = find1DFFprop(grpratemaps{m}, xxbin, fstart, fend);
               pinfo.fielddynam.lapgrpfMeanRate{i}(m) = mmrate; pinfo.fielddynam.lapgrpfPeakRate{i}(m) = pprate;
               pinfo.fielddynam.lapgrpfComX{i}(m) = comX; pinfo.fielddynam.lapgrpfPeakX{i}(m) = ploc; 
               pinfo.fielddynam.lapgrpfFieldSize{i}(m) = fsize;
           end
           pinfo.fielddynam.grpavgFstart{i} = fstart; pinfo.fielddynam.grpavgFend{i} = fend;
        else
            disp('------------> no group-averaged 1D field found');
        end
    else
        disp('------------> lap-group fields not computed: no place fields found');
    end
end

function [grpname, grplapind, grplapOccutime] = resolvegrplaps(pinfo, data, behav, bhdata, lapevName, withinlapNum, parm, cellind, xxbin)
%%%%% get lap numbers within each session
ind = find(withinlapNum == 1); nn = numel(ind); %%%assuming within event lap number always starts at 1
sessevname = cell(1, nn); sesslapind = cell(1, nn); sesswithinlap = cell(1, nn);
for (i = 1:nn-1)
    sessevname{i} = lapevName{ind(i)}; sesslapind{i} = ind(i):ind(i+1)-1; sesswithinlap{i} = withinlapNum(ind(i):ind(i+1)-1);
end
sessevname{nn} = lapevName{ind(nn)}; sesslapind{nn} = ind(nn):numel(withinlapNum); sesswithinlap{nn} = withinlapNum(ind(nn):numel(withinlapNum));
%%%%initial assignments
ng= 2*nn; grpname = cell(1,ng); grplapind = cell(1,ng); grplapOccutime = cell(1,ng); 
%%%%%get event times and occutimes
nlap = numel(lapevName); lapT.start = NaN*ones(1, nlap); lapT.ent = NaN*ones(1, nlap); lapoccutime = cell(1, nlap);
finaldirnow = pinfo.general.finaldir{cellind}; evName = pinfo.general.eventname{cellind}; evTime = data.events.eventtimes{cellind}; evType = pinfo.parm.eventtype{cellind}; 
if isfield(data.fielddynam, 'evt1DOccutAligned') %%%if aligned occutime already computed in fielddynam
   for (i=1:nlap)
       lapoccutime{i} = zeros(1, numel(xxbin));
       iii = find(strcmp(evName, lapevName{i})); 
       if (numel(iii) == 1)
          j = iii(1); lapT.start(i) = evTime{j}.start(withinlapNum(i)); lapT.ent(i) = evTime{j}.ent(withinlapNum(i));
          lapoccutime{i} = data.fielddynam.evt1DOccutAligned{cellind}{i};
       end
   end
else %%% need to figure out aligned lap-occutime on site
   refjoint = getrefjoint(pinfo, data, behav, bhdata, cellind, lapevName);
   for (i=1:nlap)
       lapoccutime{i} = zeros(1, numel(xxbin));
       iii = find(strcmp(evName, lapevName{i})); 
       if (numel(iii) == 1)
          j = iii(1); lapT.start(i) = evTime{j}.start(withinlapNum(i)); lapT.ent(i) = evTime{j}.ent(withinlapNum(i));
          %%%%locate event position data
          evSess = identifysession(evTime{j}, pinfo.general.sessionname{cellind}, pinfo.general.sessionstartT{cellind}, pinfo.general.sessionendT{cellind});
          posid = []; evid = [];
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
          if (numel(posid)==1)||(numel(evid)==1)
             occutnow = bhdata.event.LapOccuptime{posid}{evid}; xnow = bhdata.event.Xbin{posid}{evid};
             xjoint = bhdata.event.Xjoint{posid}{evid}; 
             lapoccutime{i} = realignjointtoref(xjoint, refjoint, occutnow{withinlapNum(i)}, xnow, xxbin, 'noflip', evName{j}); 
          end
       end
   end
end
%%%%%%now identify early/late groups within each session
if strncmpi(parm.splitmethod, 'lap', 2)
   for i = 1:nn
       iii = find(sesswithinlap{i} <= parm.grouplength);
       if numel(iii)>= parm.grpminlapnum
          earlyind = sesslapind{i}(iii); 
          grpname{(i-1)*2+1} = sessevname{i}; grplapind{(i-1)*2+1} = earlyind; grplapOccutime{(i-1)*2+1} = lapoccutime(earlyind); 
       end
       iii = find(max(sesswithinlap{i}) - sesswithinlap{i} <= parm.grouplength);
       if numel(iii)>= parm.grpminlapnum
          lateind = sesslapind{i}(iii); 
          grpname{(i-1)*2+2} = sessevname{i}; grplapind{(i-1)*2+2} = lateind; grplapOccutime{(i-1)*2+2} = lapoccutime(lateind); 
       end
   end
else
   for i = 1:nn
       tstartnow = lapT.start(sesslapind{i}); tendnow = lapT.ent(sesslapind{i});
       ts = min(tstartnow); te = max(tendnow);
       iii = find(tendnow-ts <= parm.grouplength); 
       if numel(iii)>= parm.grpminlapnum
          earlyind = sesslapind{i}(iii); 
          grpname{(i-1)*2+1} = sessevname{i}; grplapind{(i-1)*2+1} = earlyind; grplapOccutime{(i-1)*2+1} = lapoccutime(earlyind); 
       end
       iii = find(te - tstartnow <= parm.grouplength); 
       if numel(iii)>= parm.grpminlapnum
          lateind = sesslapind{i}(iii); 
          grpname{(i-1)*2+2} = sessevname{i}; grplapind{(i-1)*2+2} = lateind; grplapOccutime{(i-1)*2+2} = lapoccutime(lateind); 
       end
   end 
end

function occut = realignjointtoref(xjoint, refjoint, occutime, xx, refxx, flipflag, trajflag)
%%%%realign trajectories assuming the center landmarks are matched
occut = []; 
if numel(xjoint) == numel(refjoint)
   nj = numel(xjoint); 
   if strncmpi(flipflag, 'yes' ,1) %%%%assuming rate/refrate xx/refxx are opposite trajectories
       xjoint = max(xjoint) - flip(xjoint); occutime = flip(occutime);
       xx = max(xjoint)-flip(xx); %%%need to be very careful here, some xx not start at 0 or end at max(xjoint)
   end
   if (nj == 2)
       if (numel(xx) == numel(refxx))
           occut = occutime;
       elseif (numel(xx) < numel(refxx))
           occut = interpn(xx, occutime, refxx); 
       end
   else
       for (i = 2:numel(xjoint)) %%
           jj1 = find( (xx>=xjoint(i-1)) & (xx<xjoint(i))); jj2 = find( (refxx>=refjoint(i-1)) & (refxx<refjoint(i))); 
           if (numel(jj1) == numel(jj2))
               occut = [occut occutime(jj1)]; 
           else
               if (numel(jj1)>=2)
                   xnow = refxx(jj2) - refpoint(i-1);                
                   o1now = interpn(xx(jj1)-xjoint(i-1), occutime(jj1), xnow); o1now = o1now*sum(occutime(jj1))/sum(o1now); %%%to keep overall time the same
               else
                   if numel(jj1) == 1
                      o1now = occutime(jj1)*ones(1, numel(jj2))/numel(jj2); 
                   else
                       o1now = [];
                   end
               end
               occut = [occut o1now];
           end
       end
   end
else
    disp(['----------------> warning: trajectory joint points do not match for ', trajflag]);
end

function refjoint = getrefjoint(pinfo, data, behav, bhdata, cellind, lapevName)
refjoint = []; finaldirnow = pinfo.general.finaldir{cellind}; 
evName = pinfo.general.eventname{cellind}; evTime = data.events.eventtimes{cellind}; evType = pinfo.parm.eventtype{cellind};    
reflapnum = pinfo.fielddynam.runRefLapNum{cellind}; i=reflapnum(1); iii = find(strcmp(evName, lapevName{i})); 
if (numel(iii) == 1)
       j = iii(1); 
       %%%%locate event position data
       evSess = identifysession(evTime{j}, pinfo.general.sessionname{cellind}, pinfo.general.sessionstartT{cellind}, pinfo.general.sessionendT{cellind});
       posid = []; evid = [];
       if (~isempty(evSess))
           posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess) );
       else
           disp(['--------> warning: no session identified for the event: ', finaldirnow, '___', evName{j}]);  
       end
       if numel(posid) == 1
            if (isfield(behav.general, 'eventname'))
                evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
            else
                evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
            end
       else
            disp(['--------> warning: no or more than 1 positon files match the session: ', finaldirnow, '___', evSess]); 
       end
       if (numel(posid)~=1)||(numel(evid)~=1)
            disp(['--------> warning: no or more than 1 position event files match the event: ', finaldirnow, '___', evName{j}]);
       else
            refjoint = bhdata.event.Xjoint{posid}{evid}; 
       end
else
       disp(['--------> warning: no or multiple event files found match with ', lapevName{i}]);
end

function [avgratemap, occu] = findavgratemaps(lapratenow, lapOccutime, xxbin) %%%output as column vectors; xxbin - row vector
numspikes = zeros(size(xxbin))'; occu = zeros(size(xxbin))'; 
for (i = 1:numel(lapratenow))
    aa = lapratenow{i}; bb = lapOccutime{i}; if (size(aa,1)~=size(bb,1)) || (size(aa,2)~=size(bb,2)) bb = bb'; end 
    numspikes = numspikes + aa.*bb; occu = occu + bb;
end
avgratemap = numspikes ./ occu; 
%disp(['******* ' num2str([numel(numspikes) numel(occu)]) ' ******']);
function D1rate = smoothratecurvenow(D1rate, occutime, binsize, smParm)
%%now do a 1d smoothing
sigma = round(smParm.d1sigma/binsize);
if (strcmpi(smParm.d1sm, 'yes'))
    otime = ones(size(occutime));
    D1rate = OneDSmooth_new(D1rate, otime, sigma, smParm.d1Nsig); %%%smooth first - NaNs (occutime = 0) are not removed by smoothing
    %%%need to do a 1D interpolation of the NaN data points
    X = (1:numel(D1rate))'; iii = find(D1rate > -10); 
    if numel(iii)>1
       R = D1rate(iii); Y=X(iii); D1rate = interpn(Y, R, X);
    end
end
    
function pf = find1Dfieldprop(rateingrid, xbin, zerothreshold, maxgap, minpeakrate, baselinerate)
pf = []; %%%xbin - row vector; rateingrid -column vector
%%algorithm: subtract from baseline firing rate (for ctx), then examine maximum peak layer by layer 
%%break spike firing area into a number of place fields
%%work only with lineared track
%%input: rateingrid - spike firing rate in each grid on lineared track grid
%%%%%%       zerothreshold - a threshold (in percentage of average firing rate)
%%%%%%                    below which the area is not in its place field
%%%%%%       maxgap - maximum gap (rate below zero threshold) in a number of pixels within a place field
%%output: numfield - number of place fields
numfield = 0; binsize = xbin(2)-xbin(1); maxgap = round(maxgap/binsize); %%maxgap now in number of bins
gridnum = numel(rateingrid); gridleftindex = 1:gridnum; %%%this is the leftover indices in grid for new field identification
rateingrid = rateingrid - baselinerate; %iii = find(rateingrid<0); if (~isempty(iii)) rateingrid(iii) = zeros(1, numel(iii)); end
[peakratenow, peakindexnow] = max(rateingrid);
while (peakratenow >= minpeakrate) 
    numfield = numfield + 1;
    startindex = []; endindex = []; threstartindex = []; threendindex = []; %%%field locations in rateingrid
    %%%%%determine the left bound of the current field
    zeropoint = 0;
    if (peakindexnow == 1)
        startindex = 1; threstartindex = 1;
    else
        for (k = peakindexnow-1:-1:1)
            if (rateingrid(k) < peakratenow * zerothreshold) || (isnan(rateingrid(k))) 
               zeropoint = zeropoint + 1; 
            else
               zeropoint = 0;
            end
            if (zeropoint > maxgap)  
               startindex = k + floor(zeropoint/2); threstartindex = k + zeropoint; break;
            elseif ((~ismember(k, gridleftindex))||(k==1)) && (zeropoint >= floor(maxgap/2))
               startindex = k + zeropoint-floor(maxgap/2); threstartindex = k + zeropoint; break;
            elseif ((~ismember(k, gridleftindex))||(k==1)) && (zeropoint < floor(maxgap/2))
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
            if (rateingrid(k) < peakratenow * zerothreshold) || (isnan(rateingrid(k))) 
               zeropoint = zeropoint + 1; 
            else
               zeropoint = 0;
            end
            if (zeropoint >= maxgap)  
               endindex = k - floor(zeropoint/2); threendindex = k - zeropoint; break;
            elseif ((~ismember(k, gridleftindex))||(k==gridnum)) && (zeropoint >= floor(maxgap/2))
               endindex = k - zeropoint+floor(maxgap/2); threendindex = k - zeropoint; break;
            elseif ((~ismember(k, gridleftindex))||(k==gridnum)) && (zeropoint < floor(maxgap/2))
               endindex = k; threendindex = k - zeropoint; break;   
            end
        end
    end
    %%%%%%characteristics of the current field
    pf(numfield).peakX = xbin(peakindexnow); pf(numfield).startX = xbin(threstartindex); 
    pf(numfield).endX = xbin(threendindex);
    pf(numfield).boundstart = xbin(startindex); pf(numfield).boundend = xbin(endindex);
    pf(numfield).leng = (threendindex-threstartindex+1)*binsize; pf(numfield).size = binsize*sum(rateingrid(threstartindex:threendindex));
    pf(numfield).inmeanrate = baselinerate + mean(rateingrid(threstartindex:threendindex));
    pf(numfield).inpeakrate = peakratenow + baselinerate;
           %%%the center location index and skewness: ATTENTION: rateingrid is a (unnormalized) density function, but not samples
    locmean = NaN; skew = NaN; kurt = NaN;
    if (endindex-startindex > 2)
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
    end
    pf(numfield).comX = locmean; pf(numfield).skew = skew; pf(numfield).kurtosis = kurt;
    %%%%%%leave leftovers for computinng the next field
    gridleftindex = setdiff(gridleftindex, startindex:endindex); %%%leftover grid indices
    ratenow = rateingrid(gridleftindex); [peakratenow, peaknewindex] = max(ratenow); peakindexnow = gridleftindex(peaknewindex);
end

function [mrate, prate, fsize, comX, ploc, skew, kurt] = find1DFFprop(D1map, xxbin, fstart, fend)
mrate = NaN; prate = NaN; fsize=NaN; comX=NaN; ploc=NaN; skew=NaN; kurt = NaN; 
iii = find(~isnan(D1map)); 
if ~isempty(iii)
    D1map=D1map(iii); xxbin = xxbin(iii); 
    if ~isempty(find(D1map>0))
       [prate,ii] = max(D1map); ploc = xxbin(ii);
    else
        prate = 0;
    end
    indstart = 1; indend = numel(xxbin);    
for (i = 2:numel(xxbin))
    if ((xxbin(i)>fstart) && (xxbin(i-1)<=fstart)) 
       indstart = i; break
    end
end
for (i = 2:numel(xxbin))
    if ((xxbin(i)>fend) && (xxbin(i-1)<=fend)) 
       indend = i; break
    end
end
if (indend > indstart)
   ratenow = D1map(indstart:indend); xnow = xxbin(indstart:indend);
   iii = find(ratenow > -10); 
   if (numel(iii) > 2)
       ratenow = ratenow(iii); xnow = xnow (iii); mrate = mean(ratenow); %%%[prate,ii] = max(ratenow); ploc = xnow(ii); 
       fsize = sum(ratenow)*min(diff(xxbin));
       densfun = ratenow/sum(ratenow); comX = xnow * densfun; locvar = sqrt( ((xnow-comX).^2) * densfun );
       if (locvar > 0)
           X = (xnow-comX)/locvar; skew = (X.^3) * densfun; kurt = (X.^4)*densfun - 3;
       end
   end
end
end

function [crr,pp] = correlateratemaps(rate1, rate2, XY1, XY2, RC, theta)
%%%%here RC and theta are the rotation center and angle: RC[x y](coordinates in pixels), 
%%%% theta (angles in degree): 0, 90, -90, 180, 1 (x mirror), -1 (y mirror)  
%%%% assume XYgrid1 and XYgrid2 (could be 1D or 2D) have the same bin sizes
%crr = NaN; pp = NaN;
crr = 0; pp = NaN;
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

function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end

function [spinfo, spinfotime, sparsty, mrate] = findratemapprop(rate, occu)
%%%%%%another variable: mutual information is not computed here
spinfo = NaN; spinfotime = NaN; sparsty = NaN; mrate = NaN;  
[ny, nx] = size(rate); rate = reshape(rate, nx*ny, 1); occu = reshape(occu, nx*ny, 1); mr = max(rate);
if (sum(occu) > 0)
   occu = occu / sum(occu); jjj = find(rate>-1); meanrate = sum(occu(jjj).*rate(jjj)); mrate = meanrate;
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
if (nn > 0) baserate = mean(ratenow(1:nn)); end

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

